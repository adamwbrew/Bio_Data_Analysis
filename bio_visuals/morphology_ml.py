"""Exploratory, animal-balanced PLA morphology phenotypes.

These are learned measurement bins, not cell types or a running classifier.
Only log10(PLA area) and minimum/maximum diameter enter fitting. ChAT flags,
condition, intensity and location are attached AFTER fitting for description.
Every arm is fitted separately; spatial maps remain in native coordinates.
"""
import json
import warnings
import hashlib

import numpy as np
import pandas as pd
from scipy.cluster.vq import kmeans2
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

HOME = 'Home cage'
RUNNING = 'Running'
SEED = 1729
PARAMETERS = dict(features=['log10(area)', 'min_diameter / max_diameter'],
                  k_primary=3, k_sensitivity=[2, 3, 4], sample_per_animal=3000,
                  n_init=8, iterations=60, random_seed=SEED,
                  min_cluster_per_animal=10, spatial_bin=500., min_spatial_n=10,
                  scaling='Median and IQR fitted to the equal-animal sample; if IQR=0 use SD, then 1 for a constant feature',
                  invalid='Exclude nonfinite/nonpositive area, nonfinite diameters, max<=0, min<0 or min>max from ML only; zero min with positive max is retained',
                  weighting='Same training count per animal within each arm; the smallest available count caps every animal',
                  sampling='Fixed seed per animal, reused across leave-one-animal-out fits when the limiting count is unchanged; no condition balancing or outcome-aware sampling',
                  selection='K=3 is prespecified for interpretation; K=2 and K=4 are sensitivity checks, not selected by condition differences',
                  stability='Leave one animal out, refit scaler and centroids on the other five; held-out agreement to the full-data reference is descriptive stability, NOT predictive accuracy',
                  normalization='Equal-animal common phenotype mixture, with at least 10 objects and complete association in every animal per retained phenotype',
                  spatial='500-unit native bins; phenotype percentage among valid PLA objects, masked below 10 valid PLA; no cross-animal registration or tissue-area inference')


def _mean_complete(values):
    values = np.asarray(values, float)
    return float(values.mean()) if len(values) and np.isfinite(values).all() else np.nan


def _difference(frame, column):
    return (_mean_complete(frame.loc[frame.running_group.eq(HOME), column]) -
            _mean_complete(frame.loc[frame.running_group.eq(RUNNING), column]))


def _validate(objects, metadata):
    need = ['job', 'object_type', 'object_id', 'area', 'min_diameter', 'max_diameter',
            'chat_present', 'x', 'y']
    missing = set(need) - set(objects.columns)
    if missing:
        raise ValueError('Objects missing ML inputs: ' + ', '.join(sorted(missing)))
    keys = ['job', 'animal_id', 'running_group', 'pla_arm']
    if set(keys) - set(metadata.columns):
        raise ValueError('Metadata needs job, animal_id, running_group and pla_arm')
    m = metadata[keys].copy()
    if m.isna().any().any():
        raise ValueError('Missing section identity')
    m = m.astype(str)
    if m.job.duplicated().any() or m[['animal_id', 'pla_arm']].duplicated().any():
        raise ValueError('Require one primary section per animal and arm')
    if set(objects.job.astype(str)) != set(m.job):
        raise ValueError('Object jobs must exactly match primary metadata')
    for arm, group in m.groupby('pla_arm', observed=True):
        if len(group) != 6 or group.animal_id.nunique() != 6:
            raise ValueError(f'{arm}: six animals required')
        if group.running_group.value_counts().to_dict() != {HOME: 3, RUNNING: 3}:
            raise ValueError(f'{arm}: require three animals in each condition')
    if m.groupby('animal_id').running_group.nunique().gt(1).any():
        raise ValueError('Condition assignment must agree across each animal\'s arms')
    extra = [c for c in ['source_file'] if c in objects]
    p = objects.loc[objects.object_type.eq('PLA'), need + extra].copy()
    p['job'] = p.job.astype(str)
    p = p.merge(m, on='job', validate='many_to_one')
    if p[['job', 'object_id']].duplicated().any():
        raise ValueError('Object identifiers must be unique within a section')
    if 'source_file' not in p:
        source = next((c for c in ('source_file', 'filename', 'file') if c in metadata), None)
        p['source_file'] = p.job.map(metadata.assign(job=metadata.job.astype(str)).set_index('job')[source]) if source else ''
    for col in ['area', 'min_diameter', 'max_diameter', 'chat_present', 'x', 'y']:
        try:
            p[col] = p[col].to_numpy(dtype=float, na_value=np.nan)
        except (ValueError, TypeError) as exc:
            raise ValueError(f'{col} must be numeric or missing') from exc
    if not p.chat_present.dropna().isin([0, 1]).all():
        raise ValueError('ChAT flags must be binary or missing')
    finite = np.isfinite(p[['area', 'min_diameter', 'max_diameter']]).all(axis=1)
    p['valid'] = finite & p.area.gt(0) & p.max_diameter.gt(0) & p.min_diameter.ge(0) & p.min_diameter.le(p.max_diameter)
    p['shape'] = np.nan
    p.loc[p.valid, 'shape'] = p.loc[p.valid, 'min_diameter'] / p.loc[p.valid, 'max_diameter']
    p['log10_area'] = np.nan
    p.loc[p.valid, 'log10_area'] = np.log10(p.loc[p.valid, 'area'])
    # Stable row order makes sampling reproducible when source rows are shuffled.
    return p.sort_values(['job', 'object_id'], kind='mergesort').reset_index(drop=True), m


def _features(frame):
    return frame[['log10_area', 'shape']].to_numpy(float)


def _balanced_sample(frame, cap, seed):
    counts = frame.groupby('animal_id', observed=True).size()
    take = min(int(cap), int(counts.min()))
    if take < 4:
        raise ValueError('Need at least four valid morphology objects per training animal')
    pieces = []
    for animal in sorted(counts.index):
        # Fixed per-animal draws isolate animal omission from a different draw
        # of the remaining animals (unless the limiting count itself changes).
        offset = int(hashlib.sha256(str(animal).encode()).hexdigest()[:8], 16)
        rng = np.random.RandomState((int(seed) + offset) % (2 ** 32))
        part = frame.loc[frame.animal_id.eq(animal)]
        pieces.append(part.iloc[np.sort(rng.choice(len(part), take, replace=False))])
    return pd.concat(pieces), take


def _fit(frame, k, cap, seed, n_init):
    sample, take = _balanced_sample(frame, cap, seed)
    x = _features(sample)
    center = np.median(x, axis=0)
    scale = np.percentile(x, 75, axis=0) - np.percentile(x, 25, axis=0)
    sd = np.std(x, axis=0)
    scale = np.where(scale > 0, scale, np.where(sd > 0, sd, 1.))
    z = (x - center) / scale
    if len(np.unique(z, axis=0)) < k:
        raise ValueError(f'Fewer than {k} distinct valid morphology measurements')
    best = None
    for start in range(n_init):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            centroids, _ = kmeans2(z, k, iter=60, minit='++', missing='warn',
                                   seed=seed + 1009 * start)
        distances = cdist(z, centroids, metric='sqeuclidean')
        labels = np.argmin(distances, axis=1)
        if len(np.unique(labels)) != k or not np.isfinite(centroids).all():
            continue
        inertia = float(distances[np.arange(len(z)), labels].sum())
        if best is None or inertia < best[0]:
            best = inertia, centroids
    if best is None:
        raise ValueError('Clustering produced empty classes in every start')
    raw_centroids = best[1] * scale + center
    order = np.lexsort((raw_centroids[:, 1], raw_centroids[:, 0]))
    return dict(center=center, scale=scale, centroids=best[1][order],
                raw_centroids=raw_centroids[order], inertia=best[0],
                n_per_animal=take, fit_animals=sorted(frame.animal_id.unique()),
                sample=sample)


def _assign(frame, model):
    distances = cdist((_features(frame) - model['center']) / model['scale'],
                      model['centroids'], metric='sqeuclidean')
    labels = distances.argmin(axis=1)
    return labels + 1, distances[np.arange(len(frame)), labels]


def _adjusted_rand(a, b):
    """Permutation-invariant agreement, using the Hubert–Arabie correction."""
    a, b = np.asarray(a), np.asarray(b)
    if len(a) != len(b):
        raise ValueError('Partition lengths differ')
    if len(a) < 2:
        return 1.
    _, ai = np.unique(a, return_inverse=True)
    _, bi = np.unique(b, return_inverse=True)
    table = np.bincount(ai * (bi.max() + 1) + bi,
                        minlength=(ai.max() + 1) * (bi.max() + 1)).reshape(ai.max() + 1, bi.max() + 1)
    choose2 = lambda x: np.sum(np.asarray(x, float) * (np.asarray(x, float) - 1) / 2)
    index = choose2(table)
    a_pairs, b_pairs = choose2(table.sum(axis=1)), choose2(table.sum(axis=0))
    expected = a_pairs * b_pairs / (len(a) * (len(a) - 1) / 2)
    denominator = (a_pairs + b_pairs) / 2 - expected
    return float((index - expected) / denominator) if denominator else 1.


def _animal_rows(valid, labels, k, originals):
    labeled = valid.assign(cluster=labels)
    rows = []
    for job, group in labeled.groupby('job', observed=True, sort=True):
        first = group.iloc[0]
        for cluster in range(1, k + 1):
            part = group.loc[group.cluster.eq(cluster)]
            rows.append(dict(job=str(job), animal_id=first.animal_id,
                             running_group=first.running_group, pla_arm=first.pla_arm,
                             cluster=cluster, n=len(part), n_positive=float(part.chat_present.sum()) if part.chat_present.notna().all() else np.nan,
                             n_total_valid=len(group), proportion_pct=100 * len(part) / len(group),
                             association_pct=100 * _mean_complete(part.chat_present),
                             retained_pct=100 * len(group) / originals[str(job)]))
    return pd.DataFrame(rows)


def _standardize(animals, raw, k, minimum):
    """Common equal-animal mixture; fail unsupported strata, never fill zeros."""
    records, weight_rows, standardized_animals = [], [], []
    for arm, arm_rows in animals.groupby('pla_arm', observed=True, sort=True):
        arm_raw = raw.loc[raw.pla_arm.eq(arm)]
        eligible, weights = [], []
        for cluster in range(1, k + 1):
            part = arm_rows.loc[arm_rows.cluster.eq(cluster)]
            ok = len(part) == 6 and part.n.ge(minimum).all() and np.isfinite(part.association_pct).all()
            eligible.append(ok)
            weights.append(float(part.proportion_pct.mean()) / 100)
        eligible, weights = np.array(eligible), np.array(weights)
        common = np.zeros(k, float)
        if eligible.any():
            common[eligible] = weights[eligible] / weights[eligible].sum()
        for cluster in range(1, k + 1):
            weight_rows.append(dict(pla_arm=arm, k=k, cluster=cluster,
                                    eligible=bool(eligible[cluster - 1]),
                                    reference_weight=float(weights[cluster - 1]),
                                    common_weight=float(common[cluster - 1]) if eligible[cluster - 1] else np.nan))
        animal_rows = []
        for animal, group in arm_rows.groupby('animal_id', observed=True, sort=True):
            g = group.sort_values('cluster')
            first = g.iloc[0]
            kept = g.loc[eligible]
            r = dict(pla_arm=arm, k=k, animal_id=animal, job=first.job,
                     running_group=first.running_group,
                     standardized_pct=float(np.dot(common[eligible], kept.association_pct)) if eligible.any() else np.nan,
                     valid_raw_pct=float(np.dot(g.n / g.n.sum(), g.association_pct)) if np.isfinite(g.association_pct.loc[g.n.gt(0)]).all() else np.nan,
                     retained_pct=100 * kept.n.sum() / g.n.sum())
            r['eligible_raw_pct'] = (float(np.dot(kept.n / kept.n.sum(), kept.association_pct))
                                     if len(kept) and kept.n.sum() > 0 and np.isfinite(kept.association_pct).all() else np.nan)
            # Empty clusters have no rate; zero population contributes zero.
            nonempty = g.loc[g.n.gt(0)]
            r['valid_raw_pct'] = float(np.dot(nonempty.n / nonempty.n.sum(), nonempty.association_pct)) if np.isfinite(nonempty.association_pct).all() else np.nan
            animal_rows.append(r)
        animal_frame = pd.DataFrame(animal_rows)
        standardized_animals.extend(animal_rows)
        records.append(dict(pla_arm=arm, k=k,
                            raw_difference_pp=_difference(arm_raw, 'raw_association_pct'),
                            valid_difference_pp=_difference(animal_frame, 'valid_raw_pct'),
                            eligible_difference_pp=_difference(animal_frame, 'eligible_raw_pct'),
                            standardized_difference_pp=_difference(animal_frame, 'standardized_pct'),
                            home_standardized_pct=_mean_complete(animal_frame.loc[animal_frame.running_group.eq(HOME), 'standardized_pct']),
                            running_standardized_pct=_mean_complete(animal_frame.loc[animal_frame.running_group.eq(RUNNING), 'standardized_pct']),
                            retained_pct=float(animal_frame.retained_pct.mean()),
                            min_retained_pct=float(animal_frame.retained_pct.min()),
                            min_cluster_n=int(arm_rows.n.min()), eligible_clusters=int(eligible.sum())))
    return records, weight_rows, standardized_animals


def _spatial(group, width, minimum, full_group=None):
    full_group = group if full_group is None else full_group
    finite = np.isfinite(group[['x', 'y']]).all(axis=1)
    g = group.loc[finite]
    all_finite = np.isfinite(full_group[['x', 'y']]).all(axis=1)
    all_g = full_group.loc[all_finite]
    if not len(all_g):
        return None
    x0, y0 = np.floor(all_g[['x', 'y']].min().to_numpy() / width) * width
    xi = np.floor((g.x.to_numpy() - x0) / width).astype(int)
    yi = np.floor((g.y.to_numpy() - y0) / width).astype(int)
    all_xi = np.floor((all_g.x.to_numpy() - x0) / width).astype(int)
    all_yi = np.floor((all_g.y.to_numpy() - y0) / width).astype(int)
    nx, ny = int(all_xi.max() + 1), int(all_yi.max() + 1)
    counts = np.bincount(yi * nx + xi, minlength=nx * ny).reshape(ny, nx)
    all_counts = np.bincount(all_yi * nx + all_xi, minlength=nx * ny).reshape(ny, nx)
    excluded_counts = all_counts - counts
    excluded_pct = np.divide(100. * excluded_counts, all_counts, out=np.zeros_like(counts, float), where=all_counts > 0)
    fractions = {}
    for cluster in (1, 2, 3):
        keep = g.cluster.to_numpy() == cluster
        n = np.bincount((yi * nx + xi)[keep], minlength=nx * ny).reshape(ny, nx)
        pct = np.divide(100. * n, counts, out=np.zeros_like(n, float), where=counts > 0)
        fractions[str(cluster)] = [[float(v) if counts[j, i] >= minimum else None for i, v in enumerate(row)] for j, row in enumerate(pct)]
    first = group.iloc[0]
    return dict(job=str(first.job), animal_id=first.animal_id, running_group=first.running_group,
                pla_arm=first.pla_arm, width=float(width), min_count=int(minimum),
                extent=[float(x0), float(x0 + nx * width), float(y0), float(y0 + ny * width)],
                x=(x0 + (np.arange(nx) + .5) * width).tolist(),
                y=(y0 + (np.arange(ny) + .5) * width).tolist(), n=counts.tolist(),
                n_all=all_counts.tolist(), n_excluded=excluded_counts.tolist(),
                excluded_fraction_pct=[[float(v) if all_counts[j, i] >= minimum else None for i, v in enumerate(row)] for j, row in enumerate(excluded_pct)],
                cluster_fractions=fractions,
                missing_coordinate_n=int((~finite).sum()),
                all_pla_missing_coordinate_n=int((~all_finite).sum()))


def run_analysis(objects, metadata, sample_per_animal=3000, random_seed=SEED,
                 n_init=8, min_cluster_n=10, spatial_width=500., min_spatial_n=10):
    """Fit and audit exploratory morphology clusters, returning tables + JSON.

    Raw primary summaries are never modified. Labels, ChAT and coordinates have
    no path into fitting or K selection. A held-out animal is absent from BOTH
    scaling and centroid fitting; agreement with the full-data reference still
    is stability of descriptive classes, not unbiased generalization accuracy.
    """
    if any(isinstance(v, bool) or int(v) != v or v < 1 for v in
           [sample_per_animal, n_init, min_cluster_n, min_spatial_n]):
        raise ValueError('Counts and starts must be positive integers')
    if not np.isfinite(spatial_width) or spatial_width <= 0:
        raise ValueError('spatial_width must be positive')
    sample_per_animal, n_init, min_cluster_n, min_spatial_n = map(int, (sample_per_animal, n_init, min_cluster_n, min_spatial_n))
    p, m = _validate(objects, metadata)
    original_counts = p.groupby('job', observed=True).size().to_dict()
    if set(original_counts) != set(m.job):
        raise ValueError('Every primary section needs PLA objects')
    valid = p.loc[p.valid].copy()
    if set(valid.job) != set(m.job):
        raise ValueError('Every primary section needs valid morphology objects')
    audit, raw = [], []
    for job, group in p.groupby('job', observed=True, sort=True):
        first = group.iloc[0]
        identity = dict(job=str(job), animal_id=first.animal_id, running_group=first.running_group, pla_arm=first.pla_arm)
        audit.append(dict(**identity, n_pla=len(group), n_valid=int(group.valid.sum()),
                          n_excluded=int((~group.valid).sum()), retained_pct=100 * group.valid.mean(),
                          nonfinite_morphology_n=int((~np.isfinite(group[['area', 'min_diameter', 'max_diameter']]).all(axis=1)).sum()),
                          nonpositive_area_n=int(group.area.le(0).sum()),
                          nonpositive_max_diameter_n=int(group.max_diameter.le(0).sum()),
                          min_greater_than_max_n=int(group.min_diameter.gt(group.max_diameter).sum()),
                          negative_min_n=int(group.min_diameter.lt(0).sum()),
                          retained_zero_min_n=int((group.valid & group.min_diameter.eq(0)).sum()),
                          both_zero_diameters_n=int((group.min_diameter.eq(0) & group.max_diameter.eq(0)).sum()),
                          missing_association_n=int(group.chat_present.isna().sum()),
                          missing_valid_coordinates_n=int((group.valid & ~np.isfinite(group[['x', 'y']]).all(axis=1)).sum())))
        raw.append(dict(**identity, raw_association_pct=100 * _mean_complete(group.chat_present)))
    raw = pd.DataFrame(raw)
    result = dict(profiles=[], animals=[], contrasts=[], standardized=[], stability=[],
                  stability_clusters=[], audit=pd.DataFrame(audit), review_targets=[],
                  scatter=[], spatial=[], weights=[], standardized_animals=[], fit_audit=[],
                  area_only_standardized=[])
    models = []
    for arm_index, (arm, frame) in enumerate(valid.groupby('pla_arm', observed=True, sort=True)):
        base_seed = int(random_seed) + arm_index * 100000
        # Diameter-degenerate detections remain in this independent area-only
        # sensitivity. Constant shape=0 removes shape entirely from distances;
        # it does not impute missing diameters in the morphology model.
        area_frame = p.loc[p.pla_arm.eq(arm) & np.isfinite(p.area) & p.area.gt(0)].copy()
        area_frame['shape'] = 0.
        area_frame['log10_area'] = np.log10(area_frame.area)
        area_model = _fit(area_frame, 3, sample_per_animal, base_seed, n_init)
        area_labels, _ = _assign(area_frame, area_model)
        area_animals = _animal_rows(area_frame, area_labels, 3, original_counts)
        area_standard, _, _ = _standardize(area_animals, raw, 3, min_cluster_n)
        for record in area_standard:
            record['feature'] = 'log10(area) only; no diameter requirement'
            record['valid_pla_n'] = len(area_frame)
            record['morphology_pla_n'] = len(frame)
            record['input_coverage_pct'] = 100 * len(area_frame) / int(p.pla_arm.eq(arm).sum())
            record['morphology_coverage_pct'] = 100 * len(frame) / int(p.pla_arm.eq(arm).sum())
        result['area_only_standardized'].extend(area_standard)
        models.append(dict(pla_arm=arm, k=3, feature_set='area_only',
                           center=area_model['center'].tolist(), scale=area_model['scale'].tolist(),
                           centroids_log10_area_shape=area_model['raw_centroids'].tolist(),
                           training_animals=area_model['fit_animals'], n_per_animal=area_model['n_per_animal'],
                           seed=base_seed, n_init=n_init, inertia=area_model['inertia']))
        primary, labels, distance = None, None, None
        for k in (2, 3, 4):
            model = _fit(frame, k, sample_per_animal, base_seed, n_init)
            assigned, d = _assign(frame, model)
            animals = _animal_rows(frame, assigned, k, original_counts)
            standard, weights, standardized_animals = _standardize(animals, raw, k, min_cluster_n)
            result['standardized'].extend(standard)
            result['weights'].extend(weights)
            result['standardized_animals'].extend(standardized_animals)
            if k == 3:
                primary, labels, distance = model, assigned, d
                result['animals'].extend(animals.to_dict('records'))
            models.append(dict(pla_arm=arm, k=k, feature_set='area_and_shape', center=model['center'].tolist(), scale=model['scale'].tolist(),
                               centroids_log10_area_shape=model['raw_centroids'].tolist(),
                               training_animals=model['fit_animals'], n_per_animal=model['n_per_animal'],
                               seed=base_seed, n_init=n_init, inertia=model['inertia']))
        fitted = frame.assign(cluster=labels, centroid_distance=distance)
        for cluster, group in fitted.groupby('cluster', observed=True, sort=True):
            c = primary['raw_centroids'][cluster - 1]
            result['profiles'].append(dict(pla_arm=arm, cluster=int(cluster), median_area=float(group.area.median()),
                                          median_shape=float(group['shape'].median()), n=len(group),
                                          zero_min_pct=100 * float(group.min_diameter.eq(0).mean()),
                                          center_log10_area=float(c[0]), center_shape=float(c[1])))
        for animal, heldout in frame.groupby('animal_id', observed=True, sort=True):
            train = frame.loc[frame.animal_id.ne(animal)]
            fold = _fit(train, 3, sample_per_animal, base_seed, n_init)
            held_labels, _ = _assign(heldout, fold)
            reference, _ = _assign(heldout, primary)
            cost = cdist((fold['raw_centroids'] - primary['center']) / primary['scale'],
                         primary['centroids'], metric='sqeuclidean')
            row, col = linear_sum_assignment(cost)
            remap = {int(a + 1): int(b + 1) for a, b in zip(row, col)}
            matched = np.array([remap[int(label)] for label in held_labels])
            result['stability'].append(dict(pla_arm=arm, animal_id=animal,
                                            ari=_adjusted_rand(reference, held_labels),
                                            agreement_pct=100 * float(np.mean(reference == matched)),
                                            fit_animals=','.join(fold['fit_animals']), heldout_n=len(heldout),
                                            n_per_training_animal=fold['n_per_animal']))
            for cluster in (1, 2, 3):
                a, b = reference == cluster, matched == cluster
                union = int((a | b).sum())
                result['stability_clusters'].append(dict(pla_arm=arm, animal_id=animal, cluster=cluster,
                                                        jaccard=float((a & b).sum() / union) if union else np.nan,
                                                        reference_n=int(a.sum()), refit_n=int(b.sum())))
            result['fit_audit'].append(dict(pla_arm=arm, omitted_animal=animal,
                                            training_animals=','.join(fold['fit_animals']),
                                            n_training_animals=len(fold['fit_animals']),
                                            n_per_animal=fold['n_per_animal'], training_n=len(fold['sample']),
                                            heldout_n=len(heldout),
                                            scaler_log_area_center=float(fold['center'][0]),
                                            scaler_shape_center=float(fold['center'][1])))
        for job, group in fitted.groupby('job', observed=True, sort=True):
            rng = np.random.RandomState(base_seed + sum(map(ord, str(job))))
            sample = group.iloc[np.sort(rng.choice(len(group), min(500, len(group)), replace=False))]
            result['scatter'].extend(sample[['pla_arm', 'animal_id', 'cluster', 'area', 'shape', 'chat_present']].to_dict('records'))
            grid = _spatial(group, spatial_width, min_spatial_n, p.loc[p.job.eq(job)])
            if grid is not None:
                result['spatial'].append(grid)
            for cluster, part in group.groupby('cluster', observed=True, sort=True):
                for role, target in [('Typical measurement', float(part.centroid_distance.min())),
                                     ('Far-from-center measurement (99th percentile)', float(part.centroid_distance.quantile(.99)))]:
                    item = part.iloc[int(np.argmin(abs(part.centroid_distance.to_numpy() - target)))]
                    result['review_targets'].append(dict(job=str(job), animal_id=item.animal_id, running_group=item.running_group,
                                                          pla_arm=arm, cluster=int(cluster), selection=role,
                                                          object_id=item.object_id, source_file=str(item.source_file),
                                                          x=float(item.x), y=float(item.y), area=float(item.area), shape=float(item['shape']),
                                                          chat_present=float(item.chat_present),
                                                          window_xmin=float(item.x - 500), window_xmax=float(item.x + 500),
                                                          window_ymin=float(item.y - 500), window_ymax=float(item.y + 500),
                                                          reason='Inspect segmentation and ChAT overlap in the original image; selection is not grounds for automatic exclusion'))
    animals = pd.DataFrame(result['animals'])
    for (arm, cluster), group in animals.groupby(['pla_arm', 'cluster'], observed=True, sort=True):
        for metric in ('proportion_pct', 'association_pct'):
            home = _mean_complete(group.loc[group.running_group.eq(HOME), metric])
            running = _mean_complete(group.loc[group.running_group.eq(RUNNING), metric])
            result['contrasts'].append(dict(pla_arm=arm, cluster=int(cluster), metric=metric,
                                            home_pct=home, running_pct=running, difference_pp=home - running))
    for key in list(result):
        if key != 'spatial' and not isinstance(result[key], pd.DataFrame):
            result[key] = pd.DataFrame(result[key])
    result['definitions'] = dict(parameters=dict(PARAMETERS, sample_per_animal=int(sample_per_animal),
                                                random_seed=int(random_seed), n_init=int(n_init),
                                                min_cluster_per_animal=int(min_cluster_n), spatial_bin=float(spatial_width),
                                                min_spatial_n=int(min_spatial_n)),
                                 models=models, difference='Home cage minus Running; reverse is exact negation',
                                 exclusions_scope='ML only. Primary unfiltered association remains unchanged; raw_difference_pp retains its full denominator.',
                                 uncertainty='Six animals, three per condition. No inferential p-values, predictive claims, cell identities or causal mechanisms are inferred.',
                                 interpretation='Morphology clusters are descriptive measurement partitions. They may reflect segmentation or section/acquisition effects.',
                                 profiles='Medians use all assigned valid objects; they are descriptive pooled profiles. Fit sampling is balanced by animal.',
                                 source_review=f"{len(result['review_targets'])} native-coordinate review targets: one typical and one 99th-distance-percentile detection per animal/arm/phenotype. Original microscopy images are unavailable.",
                                 alignment='Native maps only; unapproved image alignment never enters fitting, normalization or stability checks.')
    # Enforce serializable, finite model definitions independently of the views.
    json.dumps(result['definitions'], allow_nan=False)
    return result

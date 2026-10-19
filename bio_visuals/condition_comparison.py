"""Animal-level Home-cage minus Running comparisons for the HALO handoff.

Object rows are subsamples, not independent experimental replicates. The two
PLA arms remain separate; no spatial alignment or anatomical subtraction is
performed here. All weighting and small-sample diagnostics are explicit.
"""
import itertools
import json

import numpy as np
import pandas as pd

from .handoff import ANIMAL_ORDER, ARMS, GROUPS, SIZE_EDGES, SIZE_LABELS

HOME = 'Home cage'
RUNNING = 'Running'
PRIMARY = 'association_pct'
METRICS = {
    'association_pct': dict(label='PLA objects with ChAT signal', unit='percentage points',
                            caution='Object association, not neuron counts, binding strength or metabolic flux.'),
    'pla_count': dict(label='PLA objects per section', unit='objects',
                     caution='Unnormalized export count; sampled area, section level and detection can differ.'),
    'chat_count': dict(label='ChAT objects per section', unit='objects',
                      caution='Segmented detections, not neurons; sampled tissue area is unavailable.'),
    'mean_pla_area_um2': dict(label='Mean PLA object area', unit='µm²',
                            caution='Sensitive to large detections and segmentation.'),
    'median_pla_area_um2': dict(label='Median PLA object area', unit='µm²',
                              caution='A per-animal median; group summaries average those medians.'),
    'large_pla_pct': dict(label='PLA objects with area ≥1 µm²', unit='percentage points',
                         caution='Prespecified size sensitivity measure, not a validated biological class.'),
    'overlap_area_ratio_pct': dict(label='Actual ChAT-overlap area / total PLA area', unit='percentage points',
                                  caution='Sum of exported overlap areas / sum of PLA areas within each animal; not area of all positive objects.'),
}
ESTIMATORS = {
    'raw': 'All PLA objects (primary)',
    'restricted': 'Eligible size strata only, original within-animal mix',
    'standardized': 'Eligible size strata, same pooled size mix for every animal',
}


def _complete_mean(values):
    values = np.asarray(values, float)
    return float(values.mean()) if len(values) and np.isfinite(values).all() else np.nan


def _complete_sum(values):
    values = np.asarray(values, float)
    return float(values.sum()) if np.isfinite(values).all() else np.nan


def _ratio(num, den):
    return 100 * num / den if np.isfinite(num) and np.isfinite(den) and den > 0 else np.nan


def _validated_inputs(objects, metadata):
    required_meta = ['job', 'animal_id', 'running_group', 'pla_arm']
    required_objects = ['job', 'object_type', 'area', 'chat_present', 'overlap_area']
    for frame, required, name in ((metadata, required_meta, 'metadata'),
                                  (objects, required_objects, 'objects')):
        missing = set(required) - set(frame.columns)
        if missing:
            raise ValueError(f'{name} is missing required columns: {sorted(missing)}')
    m = metadata[required_meta].copy()
    if m.isna().any().any():
        raise ValueError('Missing section identity or condition assignment')
    for column in required_meta:
        m[column] = m[column].astype(str)
    if len(m) != 12 or m.job.duplicated().any() or m[['animal_id', 'pla_arm']].duplicated().any():
        raise ValueError('Expected one section per animal per arm: twelve unique section jobs')
    if set(m.pla_arm) != set(ARMS.values()):
        raise ValueError('Expected the two distinct handoff PLA arms')
    for arm, rows in m.groupby('pla_arm', observed=True):
        if set(rows.animal_id) != set(GROUPS) or len(rows) != 6:
            raise ValueError(f'{arm}: expected the six handoff animals exactly once')
        if not rows.running_group.eq(rows.animal_id.map(GROUPS)).all():
            raise ValueError(f'{arm}: group assignments differ from the scientist handoff')
    if objects.job.isna().any() or set(objects.job.astype(str)) != set(m.job):
        raise ValueError('Object jobs must match all twelve primary sections; exclude replacement duplicates first')
    if objects.object_type.isna().any() or not objects.object_type.isin(['PLA', 'ChAT']).all():
        raise ValueError('Unknown or missing object type')
    # Missing measurements remain unknown. Strings and nonbinary flags must not
    # be silently converted to truth values or interpreted as missing negatives.
    known_flags = objects.chat_present.dropna()
    if not known_flags.isin([True, False, 0, 1]).all():
        raise ValueError('ChAT-presence values must be binary or missing')
    for column in ('area', 'overlap_area'):
        try:
            values = objects[column].to_numpy(dtype=float, na_value=np.nan)
        except (TypeError, ValueError) as exc:
            raise ValueError(f'{column} must contain numeric measurements or missing values') from exc
        if np.isinf(values).any() or (values[np.isfinite(values)] < 0).any():
            raise ValueError(f'{column} must be nonnegative, finite measurements or missing values')
    return m


def _contrast(frame, metric, arm):
    h = frame.loc[frame.running_group.eq(HOME)].sort_values('animal_id')
    r = frame.loc[frame.running_group.eq(RUNNING)].sort_values('animal_id')
    hv, rv = h[metric].to_numpy(float), r[metric].to_numpy(float)
    complete = np.isfinite(np.r_[hv, rv]).all()
    loo = []
    for hi, ri in itertools.product(range(3), repeat=2):
        value = float(np.delete(hv, hi).mean() - np.delete(rv, ri).mean()) if complete else np.nan
        loo.append(dict(pla_arm=arm, metric=metric, omitted_home=h.iloc[hi].animal_id,
                        omitted_running=r.iloc[ri].animal_id, home_minus_running=value))
    row = dict(pla_arm=arm, metric=metric, unit=METRICS.get(metric, METRICS[PRIMARY])['unit'],
               home_mean=_complete_mean(hv), running_mean=_complete_mean(rv),
               home_minus_running=_complete_mean(hv) - _complete_mean(rv),
               n_home=len(h), n_running=len(r), n_finite_home=int(np.isfinite(hv).sum()),
               n_finite_running=int(np.isfinite(rv).sum()),
               leave_pair_out_min=min(x['home_minus_running'] for x in loo) if complete else np.nan,
               leave_pair_out_max=max(x['home_minus_running'] for x in loo) if complete else np.nan)
    return row, loo


def comparison_results(objects, metadata, min_per_stratum=10):
    """Return named data tables and definitions for a descriptive comparison.

    ``min_per_stratum`` is the prespecified minimum PLA count in EACH of the
    six animals for a common size stratum (default 10). Size weights are the
    equal-animal mean of each animal's size proportions, renormalized across
    eligible strata. Standardization is a sensitivity analysis; it neither
    controls confounding nor establishes a mediator or causal explanation.

    Group means require all three animal measurements. Missing values are
    preserved, never zero-imputed or silently dropped. The exact enumeration
    is only for the primary unfiltered association metric and only meaningful
    under exchangeable animal labels; the handoff does not establish that.
    """
    if isinstance(min_per_stratum, bool) or int(min_per_stratum) != min_per_stratum or min_per_stratum < 1:
        raise ValueError('min_per_stratum must be a positive integer')
    min_per_stratum = int(min_per_stratum)
    m = _validated_inputs(objects, metadata)
    animals, strata = [], []
    keys = m.set_index('job')
    for job, frame in objects.groupby('job', observed=True, sort=True):
        job = str(job)
        identity = dict(job=job, **keys.loc[job].to_dict())
        p = frame.loc[frame.object_type.eq('PLA')]
        area = p.area.to_numpy(dtype=float, na_value=np.nan)
        positive = p.chat_present.to_numpy(dtype=float, na_value=np.nan)
        overlap = p.overlap_area.to_numpy(dtype=float, na_value=np.nan)
        full_area = len(area) > 0 and np.isfinite(area).all()
        total_area = _complete_sum(area)
        animals.append(dict(**identity, association_pct=100 * _complete_mean(positive),
                            pla_count=len(p), chat_count=int(frame.object_type.eq('ChAT').sum()),
                            mean_pla_area_um2=_complete_mean(area),
                            median_pla_area_um2=float(np.median(area)) if full_area else np.nan,
                            large_pla_pct=100 * float(np.mean(area >= 1)) if full_area else np.nan,
                            overlap_area_ratio_pct=_ratio(_complete_sum(overlap), total_area),
                            missing_area_count=int(np.isnan(area).sum()),
                            missing_association_count=int(np.isnan(positive).sum())))
        for index, label in enumerate(SIZE_LABELS):
            use = (area >= SIZE_EDGES[index]) & (area < SIZE_EDGES[index + 1])
            pos = positive[use]
            strata.append(dict(**identity, size_bin=label, size_index=index, n_pla=int(use.sum()),
                               n_positive=_complete_sum(pos), fraction_pct=100 * _complete_mean(pos),
                               size_share_pct=100 * use.sum() / len(p) if len(p) else np.nan))
    animal_metrics, size_animals = pd.DataFrame(animals), pd.DataFrame(strata)
    contrasts, leave_pair_out, cross_animal, allocations, exact_summary = [], [], [], [], []
    size_groups, size_contrasts, weights, standardization, standardized_animals = [], [], [], [], []
    for arm, animal_frame in animal_metrics.groupby('pla_arm', observed=True, sort=True):
        for metric in METRICS:
            row, loo = _contrast(animal_frame, metric, arm)
            contrasts.append(row)
            leave_pair_out.extend(loo)
        home = animal_frame.loc[animal_frame.running_group.eq(HOME)].sort_values('animal_id')
        running = animal_frame.loc[animal_frame.running_group.eq(RUNNING)].sort_values('animal_id')
        for h, r in itertools.product(home.to_dict('records'), running.to_dict('records')):
            cross_animal.append(dict(pla_arm=arm, home_animal=h['animal_id'], running_animal=r['animal_id'],
                                     home_value=h[PRIMARY], running_value=r[PRIMARY],
                                     home_minus_running=h[PRIMARY] - r[PRIMARY]))
        ordered = animal_frame.set_index('animal_id').loc[ANIMAL_ORDER]
        primary_values = ordered[PRIMARY].to_numpy(float)
        observed = _complete_mean(home[PRIMARY]) - _complete_mean(running[PRIMARY])
        finite = np.isfinite(primary_values).all()
        null_differences = []
        for choice in itertools.combinations(range(6), 3):
            complement = [i for i in range(6) if i not in choice]
            value = float(primary_values[list(choice)].mean() - primary_values[complement].mean()) if finite else np.nan
            null_differences.append(value)
            allocations.append(dict(pla_arm=arm, home_animals=','.join(ANIMAL_ORDER[i] for i in choice),
                                    running_animals=','.join(ANIMAL_ORDER[i] for i in complement),
                                    home_minus_running=value,
                                    is_observed=set(ANIMAL_ORDER[i] for i in choice) == set(home.animal_id)))
        # Both the observed allocation and its complete reversal are enumerated;
        # there is no Monte Carlo +1 correction and minimum two-sided p is 2/20.
        p_value = float(np.mean(np.abs(null_differences) >= abs(observed) - 1e-12)) if finite else np.nan
        exact_summary.append(dict(pla_arm=arm, observed_difference=observed, two_sided_p=p_value,
                                  allocations=20, min_attainable_p=0.1,
                                  status='Exploratory label-enumeration diagnostic; exchangeability unverified' if finite
                                         else 'Unavailable: at least one animal association is missing'))
        s = size_animals.loc[size_animals.pla_arm.eq(arm)]
        eligible = []
        pooled_weights = []
        reasons = []
        for index, label in enumerate(SIZE_LABELS):
            part = s.loc[s.size_index.eq(index)]
            unsupported = part.loc[(part.n_pla < min_per_stratum) | ~np.isfinite(part.fraction_pct), 'animal_id'].tolist()
            ok = not unsupported
            eligible.append(ok)
            pooled_weights.append(_complete_mean(part.size_share_pct) / 100)
            reasons.append('All six animals have sufficient observed PLA support' if ok else
                           'Insufficient count or missing association in: ' + ', '.join(unsupported))
            group_rows = []
            for group in (HOME, RUNNING):
                g = part.loc[part.running_group.eq(group)]
                record = dict(pla_arm=arm, running_group=group, size_bin=label, size_index=index,
                              mean_fraction_pct=_complete_mean(g.fraction_pct),
                              mean_size_share_pct=_complete_mean(g.size_share_pct),
                              min_animal_n=int(g.n_pla.min()), max_animal_n=int(g.n_pla.max()))
                size_groups.append(record)
                group_rows.append(record)
            size_contrasts.append(dict(pla_arm=arm, size_bin=label, size_index=index,
                                       home_mean=group_rows[0]['mean_fraction_pct'],
                                       running_mean=group_rows[1]['mean_fraction_pct'],
                                       home_minus_running=group_rows[0]['mean_fraction_pct'] - group_rows[1]['mean_fraction_pct'],
                                       eligible=ok))
        eligible = np.asarray(eligible, bool)
        pooled_weights = np.asarray(pooled_weights, float)
        weight_total = pooled_weights[eligible].sum()
        common_weights = np.full(len(SIZE_LABELS), np.nan)
        if eligible.any() and np.isfinite(weight_total) and weight_total > 0:
            common_weights[eligible] = pooled_weights[eligible] / weight_total
        else:
            eligible[:] = False
        for index, label in enumerate(SIZE_LABELS):
            weights.append(dict(pla_arm=arm, size_bin=label, size_index=index, included=bool(eligible[index]),
                                reason=reasons[index], pooled_weight=pooled_weights[index], common_weight=common_weights[index],
                                min_animal_n=int(s.loc[s.size_index.eq(index), 'n_pla'].min())))
        arm_standard = []
        for animal in animal_frame.to_dict('records'):
            part = s.loc[s.animal_id.eq(animal['animal_id'])].sort_values('size_index')
            keep = part.loc[eligible]
            retained_count = int(keep.n_pla.sum())
            retained_pct = _ratio(retained_count, animal['pla_count'])
            restricted = _ratio(_complete_sum(keep.n_positive), retained_count)
            standardized = float(np.dot(common_weights[eligible], keep.fraction_pct.to_numpy(float))) if eligible.any() else np.nan
            for estimator, value in (('raw', animal[PRIMARY]), ('restricted', restricted), ('standardized', standardized)):
                record = dict(pla_arm=arm, animal_id=animal['animal_id'], running_group=animal['running_group'],
                              estimator=estimator, association_pct=value,
                              retained_pct=100. if estimator == 'raw' and animal['pla_count'] else retained_pct)
                arm_standard.append(record)
                standardized_animals.append(record)
        arm_standard = pd.DataFrame(arm_standard)
        for estimator in ESTIMATORS:
            data = arm_standard.loc[arm_standard.estimator.eq(estimator)]
            contrast, _ = _contrast(data, PRIMARY, arm)
            standardization.append(dict(pla_arm=arm, estimator=estimator, home_mean=contrast['home_mean'],
                                        running_mean=contrast['running_mean'], home_minus_running=contrast['home_minus_running'],
                                        leave_pair_out_min=contrast['leave_pair_out_min'], leave_pair_out_max=contrast['leave_pair_out_max'],
                                        eligible_strata='; '.join(np.asarray(SIZE_LABELS)[eligible]),
                                        excluded_strata='; '.join(np.asarray(SIZE_LABELS)[~eligible]),
                                        min_retained_pct=float(data.retained_pct.min()), max_retained_pct=float(data.retained_pct.max())))
    definitions = dict(direction='Home cage minus Running: positive means higher in non-running animals',
                       groups=[HOME, RUNNING], primary_metric=PRIMARY, metrics=METRICS, estimators=ESTIMATORS,
                       biological_unit='Animal; three animals per condition and PLA arm',
                       weighting='Every animal has equal weight; never pool object rows between animals',
                       uncertainty='Leave-one-animal-out from EACH group gives nine two-versus-two differences. The range is a sensitivity check, not a confidence interval.',
                       cross_animal='All nine Home-cage–Running differences share animals; they are not nine independent replicates or matched pairs.',
                       exact_test='All 20 allocations of six animals to groups of three, two-sided absolute mean difference. Meaningful only if animal labels are exchangeable; randomization and batch independence are not established here. Minimum attainable two-sided p is 0.10. Not a confirmatory significance test.',
                       size_edges_um2=[0, .03, .1, .3, 1, None], size_labels=SIZE_LABELS,
                       min_per_stratum=min_per_stratum,
                       size_weights='Equal-animal mean of size proportions (count in stratum / all PLA objects), renormalized over strata with at least the minimum count and observed association in all six animals.',
                       size_interpretation='Raw is primary. Restricted isolates changes from excluding unsupported sizes. Standardized uses one common size mixture; this sensitivity analysis does not establish why groups differ or adjust all confounding. Its leave-pair-out range holds the full-cohort eligibility and common weights fixed.',
                       spatial_status='No pixelwise or anatomical subtraction: landmarks, orientation, scale and tissue support require review.',
                       missing='Missing measurements stay missing; group means require all three animals. Unknown area is outside known size strata and reduces retained coverage.')
    return dict(definitions=definitions, animal_metrics=animal_metrics, contrasts=pd.DataFrame(contrasts),
                cross_animal=pd.DataFrame(cross_animal), leave_pair_out=pd.DataFrame(leave_pair_out),
                exact_allocations=pd.DataFrame(allocations), exact_summary=pd.DataFrame(exact_summary),
                size_animals=size_animals, size_groups=pd.DataFrame(size_groups), size_contrasts=pd.DataFrame(size_contrasts),
                standard_weights=pd.DataFrame(weights), standardization=pd.DataFrame(standardization),
                standardized_animals=pd.DataFrame(standardized_animals))


def comparison_payload(results):
    """Convert results to strict JSON-ready records, with unknown values null."""
    output = {}
    for name, value in results.items():
        output[name] = json.loads(value.to_json(orient='records', force_ascii=False)) if isinstance(value, pd.DataFrame) else value
    # Do not let an added non-serializable field or NaN enter browser payloads.
    json.dumps(output, allow_nan=False)
    return output


def condition_registration(objects, reviews, width=1000., min_pla=10):
    """Expose reviewed group maps, or an explicit incomplete-review status.

    This adapter keeps the comparison direction consistent with all new
    condition results. It never runs registration on partial, duplicated or
    unapproved section sets. Unobserved/unsupported bins serialize as null.
    """
    from . import handoff as h

    status = h.review_status(reviews)  # Review before any coordinate computation.
    total = len(reviews)
    approved = int(status.ready.sum()) if len(status) else 0
    needs_review = status.loc[~status.ready].to_dict('records') if len(status) else []
    result = dict(ready=False, approved=approved, total=total, needs_review=needs_review,
                  message='', maps={})
    if total != 12:
        result['message'] = (f'Group spatial subtraction is unavailable: received {total} of the required '
                             '12 section reviews. All six animals in both PLA arms are required.')
        return result
    identities = [(r['animal_id'], r['pla_arm']) for r in reviews]
    expected = set(itertools.product(GROUPS, ARMS.values()))
    jobs = [str(r['job']) for r in reviews]
    correct_groups = all(r.get('running_group') == GROUPS.get(r['animal_id']) for r in reviews)
    if (len(set(jobs)) != 12 or set(identities) != expected or not correct_groups or
            set(objects.job.astype(str)) != set(jobs)):
        result['message'] = ('Group spatial subtraction is unavailable: section identities, animal assignments '
                             'or raw-data jobs do not match the twelve distinct handoff sections.')
        return result
    if approved != 12:
        result['message'] = (f'Group spatial subtraction is unavailable: {approved}/12 sections pass landmark, '
                             'orientation, common-scale and support review. Whole-section condition comparisons remain available.')
        return result
    if not np.isfinite(width) or width <= 0:
        raise ValueError('width must be positive and finite')
    if isinstance(min_pla, bool) or not np.isfinite(min_pla) or int(min_pla) != min_pla or min_pla < 1:
        raise ValueError('min_pla must be a positive integer')
    registered = h.registered_fraction_maps(objects, reviews, width=width, min_pla=int(min_pla))

    def nullable(value):
        array = np.asarray(value)
        return np.where(np.isfinite(array), array, None).tolist()

    maps = {}
    for arm, data in registered['arms'].items():
        maps[arm] = dict(x_edges=nullable(registered['x_edges']), y_edges=nullable(registered['y_edges']),
                         home_mean=nullable(data[HOME]['mean']), running_mean=nullable(data[RUNNING]['mean']),
                         home_n=nullable(data[HOME]['n']), running_n=nullable(data[RUNNING]['n']),
                         home_minus_running=nullable(-np.asarray(data['difference'], float)))
    result.update(ready=True, maps=maps,
                  message='Reviewed spatial association: Home cage minus Running, with equal-animal means and three animals supporting each condition bin. Blank bins are unknown, not zero.')
    json.dumps(result, allow_nan=False)
    return result

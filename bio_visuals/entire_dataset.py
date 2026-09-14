"""Validated multi-image HALO exports; repeat analyses remain separate.

Source-image identity and literal section labels are provenance, not treatment
assignments or independent animal identifiers. No coordinate registration is
inferred from the CSVs. Explicit replacement decisions are auditable.
"""
from pathlib import Path, PureWindowsPath
import csv
import hashlib
import re
import numpy as np
import pandas as pd
from .data import MEASUREMENTS, OBJECT_TYPES, SOURCE_COLUMNS, quality_control, summarize

JOB_RE = re.compile(r'job(?P<job>\d+)PLA')
IMAGE_RE = re.compile(r'(?P<image>20260515_\d+_\d+_[^_]+_Ellie_PLA_05152026_\d+)')
REGION_RE = re.compile(r'(?P<condition>\d+[MF])_(?:TDP43-(?P<pla_a>Pfkp|Hk1)|(?P<pla_b>Pfkp|Hk1)-TDP43)_tight')
TRUNCATED_RERUN = {'4978': '5403'}
SIGNATURE_FIELDS = ['object_type', 'xmin', 'xmax', 'ymin', 'ymax', 'area',
                    'inner_area', 'outer_area', 'overlap_area', 'overlap_pct',
                    'intensity', 'min_diameter', 'max_diameter', 'median_diameter',
                    'pla_present', 'chat_present']


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _metadata(path):
    path = Path(path)
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != SOURCE_COLUMNS:
            raise ValueError(f'Unexpected CSV schema in {path.name}')
        first = next(reader, None)
    if first is None:
        raise ValueError(f'Empty export: {path.name}')
    job_match = JOB_RE.search(path.name)
    image_match = IMAGE_RE.search(PureWindowsPath(first['Image Location']).name)
    region_match = REGION_RE.fullmatch(first['Analysis Region'])
    if not all((job_match, image_match, region_match)):
        raise ValueError(f'Cannot parse image/job/region metadata: {path.name}')
    if not path.name.startswith(image_match.group('image')):
        raise ValueError(f'Filename and image metadata disagree: {path.name}')
    return dict(job=job_match.group('job'), source_file=path.name,
                image_location=first['Image Location'], image_id=image_match.group('image'),
                region=first['Analysis Region'], condition=region_match.group('condition'),
                pla_type=region_match.group('pla_a') or region_match.group('pla_b'),
                algorithm=first['Algorithm Name'])


def discover_exports(root='.', verify_hashes=True, replacements=None):
    """Inventory all exports; only reviewed replacements may be excluded.

    Unlike a largest-file heuristic, an unknown repeated region stops selection
    until its provenance has been reviewed. Missing replacements also stop it.
    """
    root = Path(root).expanduser().resolve()
    replacements = TRUNCATED_RERUN if replacements is None else replacements
    paths = sorted(root.glob('*.csv'))
    if not paths:
        raise FileNotFoundError(f'No HALO CSV exports found in {root}')
    records = []
    for path in paths:
        record = _metadata(path)
        with path.open(encoding='utf-8-sig', newline='') as handle:
            reader = csv.reader(handle)
            next(reader)
            rows = sum(1 for _ in reader)
        record.update(rows=rows, bytes=path.stat().st_size,
                      sha256=_sha256(path) if verify_hashes else None)
        records.append(record)
    manifest = pd.DataFrame(records).sort_values(['image_id', 'condition', 'pla_type', 'job']).reset_index(drop=True)
    if manifest.job.duplicated().any():
        raise ValueError('A job appears more than once in the complete folder')
    manifest['duplicate_group'] = manifest.image_id + '|' + manifest.condition + '|' + manifest.pla_type
    manifest['is_primary'] = True
    manifest['duplicate_status'] = 'primary'
    manifest['replacement_job'] = ''
    by_job = manifest.set_index('job')
    for short, full in replacements.items():
        if short not in by_job.index:
            continue
        if full not in by_job.index:
            raise ValueError(f'Interrupted job {short} requires replacement {full}')
        if by_job.loc[short, 'duplicate_group'] != by_job.loc[full, 'duplicate_group']:
            raise ValueError('Replacement must reference the same image and region')
        if by_job.loc[full, 'rows'] <= by_job.loc[short, 'rows']:
            raise ValueError('Reviewed replacement no longer has the expected larger row count')
        manifest.loc[manifest.job == short, ['is_primary', 'duplicate_status', 'replacement_job']] = [False, 'interrupted_superseded', full]
        manifest.loc[manifest.job == full, 'duplicate_status'] = 'primary_replacement'
    if manifest.loc[manifest.is_primary, 'duplicate_group'].duplicated().any():
        raise ValueError('Unreviewed repeated region: specify a validated replacement decision')
    image_ids = sorted(manifest.image_id.unique())
    manifest['image_label'] = manifest.image_id.map({v: f'Image {i+1}' for i, v in enumerate(image_ids)})
    return manifest


def load_export(path, record=None):
    """Validate every row, including constant metadata, exact IDs and flags."""
    path = Path(path)
    record = _metadata(path) if record is None else record
    dtypes = {name: 'float64' for name in ['Object Id', *MEASUREMENTS]}
    frame = pd.read_csv(path, encoding='utf-8-sig', dtype=dtypes)
    if frame.columns.tolist() != SOURCE_COLUMNS:
        raise ValueError(f'Unexpected schema in job {record["job"]}')
    if frame.empty or frame.isna().any().any() or not np.isfinite(frame[list(dtypes)].to_numpy()).all():
        raise ValueError(f'Missing/nonfinite measurements in job {record["job"]}')
    for column, key in [('Analysis Region', 'region'), ('Image Location', 'image_location'), ('Algorithm Name', 'algorithm')]:
        if frame[column].nunique() != 1 or str(frame[column].iloc[0]) != record[key]:
            raise ValueError(f'Inconsistent {column} in job {record["job"]}')
    ids = frame['Object Id'].to_numpy()
    if ((ids < 0) | (ids >= 2**53) | (ids != np.floor(ids))).any() or frame['Object Id'].duplicated().any():
        raise ValueError(f'Invalid or duplicate Object Id in job {record["job"]}')
    for column in ['PLA - TRITC_E1 present', 'ChAt-50 - FITC_E1 present']:
        if not frame[column].isin([0, 1]).all():
            raise ValueError(f'Nonbinary flag in job {record["job"]}: {column}')
    if not frame['Object Type'].isin(OBJECT_TYPES).all():
        raise ValueError(f'Unknown object type in job {record["job"]}')
    frame = frame.drop(columns=['Image Location', 'Algorithm Name']).rename(columns={
        **MEASUREMENTS, 'Analysis Region': 'region', 'Object Id': 'object_id', 'Object Type': 'object_type'})
    frame['object_type'] = frame.object_type.map(OBJECT_TYPES)
    frame['object_id'] = frame.object_id.astype('int64')
    for field in ['pla_present', 'chat_present']:
        frame[field] = frame[field].astype(bool)
    for field in ['job', 'image_id', 'condition', 'pla_type']:
        frame[field] = record[field]
    frame['image_label'] = record.get('image_label', record['image_id'])
    frame['roi_mode'] = 'outline'
    frame['pair_id'] = frame.condition + '_' + frame.pla_type + '-TDP43'
    frame['x'] = (frame.xmin + frame.xmax) / 2
    frame['y'] = (frame.ymin + frame.ymax) / 2
    for field in ['job', 'image_id', 'condition', 'pla_type', 'image_label', 'roi_mode', 'pair_id', 'region', 'object_type']:
        frame[field] = frame[field].astype('category')
    return frame


def load_complete_dataset(root='.', verify_hashes=True):
    """Return primary objects and a manifest auditing all 13 CSVs."""
    root = Path(root).expanduser().resolve()
    manifest = discover_exports(root, verify_hashes)
    frames, audits = [], []
    for record in manifest.to_dict('records'):
        frame = load_export(root / record['source_file'], record)
        if len(frame) != record['rows']:
            raise ValueError(f'Row count changed while loading {record["job"]}')
        audit = quality_control(frame).iloc[0].to_dict()
        audits.append({k: v for k, v in audit.items() if k not in ('region', 'n_objects')})
        if record['is_primary']:
            frames.append(frame)
    objects = pd.concat(frames, ignore_index=True)
    for field in ['job', 'image_id', 'condition', 'pla_type', 'image_label', 'roi_mode', 'pair_id', 'region', 'object_type']:
        objects[field] = objects[field].astype('category')
    manifest = manifest.merge(pd.DataFrame(audits), on='job', validate='one_to_one')
    return objects, manifest


def compare_initial_files(manifest, initial_root):
    """Check byte identity of the copied outlines; never concatenate copies."""
    records = []
    for row in manifest.to_dict('records'):
        path = Path(initial_root) / row['source_file']
        if path.exists():
            if not row['sha256']:
                raise ValueError('Hash verification is required to establish copied-file identity')
            digest = _sha256(path)
            records.append(dict(job=row['job'], source_file=row['source_file'],
                                initial_sha256=digest, complete_sha256=row['sha256'],
                                identical_file=digest == row['sha256'], rows=row['rows']))
    return pd.DataFrame(records)


def signature_overlap(short, full):
    """Exact 16-field multiset comparison; no hash-only or Object Id match."""
    a = short.groupby(SIGNATURE_FIELDS, observed=True, dropna=False).size().rename('short_count')
    b = full.groupby(SIGNATURE_FIELDS, observed=True, dropna=False).size().rename('full_count')
    counts = pd.concat([a, b], axis=1).fillna(0)
    matched = int(counts.min(axis=1).sum())
    return dict(short_rows=len(short), full_rows=len(full), identical_measurement_rows=matched,
                short_unmatched_rows=len(short)-matched, full_additional_rows=len(full)-matched)


def run_summary(objects):
    """One row per export, including counts, overlap and distribution tails."""
    stats = summarize(objects)
    pla = stats.loc[stats.object_type == 'PLA'].copy()
    chat = stats.loc[stats.object_type == 'ChAT', ['job', 'n_objects']].rename(columns={'n_objects': 'chat_objects'})
    pla = pla.rename(columns={'n_objects': 'pla_objects', 'other_marker_count': 'pla_with_chat',
                              'other_marker_pct': 'pla_association_pct', 'mean_object_overlap_pct': 'mean_overlap_pct'})
    meta = objects[['job', 'image_id', 'image_label', 'condition', 'pla_type']].drop_duplicates()
    if meta.job.duplicated().any():
        raise ValueError('Each job must identify one source image, literal label and PLA type')
    result = pla.merge(chat, on='job', validate='one_to_one').merge(meta, on='job', validate='one_to_one')
    extents = objects.groupby('job', observed=True).agg(x_min=('xmin','min'), x_max=('xmax','max'),
                                                     y_min=('ymin','min'), y_max=('ymax','max')).reset_index()
    result = result.merge(extents, on='job', validate='one_to_one')
    result['all_objects'] = result.pla_objects + result.chat_objects
    return result.sort_values(['image_id','condition','pla_type']).reset_index(drop=True)


def type_summary(summary):
    """Distinguish equal-export averaging from pooled detection weighting."""
    rows = []
    for typ, group in summary.groupby('pla_type', observed=True):
        rows.append(dict(pla_type=typ, exports=len(group), pla_objects=int(group.pla_objects.sum()),
                         mean_export_association_pct=group.pla_association_pct.mean(),
                         pooled_object_association_pct=100*group.pla_with_chat.sum()/group.pla_objects.sum(),
                         min_export_association_pct=group.pla_association_pct.min(),
                         max_export_association_pct=group.pla_association_pct.max(),
                         mean_export_weighted_overlap_pct=group.area_weighted_overlap_pct.mean()))
    return pd.DataFrame(rows)


def metadata_template(summary):
    frame = summary[['job','image_id','condition','pla_type']].astype(str).rename(columns={'condition':'sample_label'})
    for field in ['animal_id','running_group','section_id','anatomical_region','anatomical_level',
                  'species','tissue','section_relationship','acquisition_batch','notes']:
        frame[field] = ''
    return frame


def biological_summaries(summary, metadata):
    """Average sections within animal/anatomy/type before averaging animals.

    Missing design assignments return a review table instead of inventing
    treatment contrasts. No p-values are fitted automatically.
    """
    if metadata.job.duplicated().any() or set(metadata.job.astype(str)) != set(summary.job.astype(str)):
        raise ValueError('Metadata must contain exactly one row for every primary job')
    design = summary[['job','condition','pla_type','image_id']].astype(str).rename(columns={'condition':'sample_label'})
    checked = design.merge(metadata.fillna('').astype(str), on='job', suffixes=('_expected',''), validate='one_to_one')
    for field in ['sample_label','pla_type','image_id']:
        if (checked[field] != checked[field+'_expected']).any():
            raise ValueError(f'Metadata identity does not match the export: {field}')
    required = ['animal_id','running_group','section_id','anatomical_region']
    missing = checked[required].apply(lambda column: column.str.strip().eq(''))
    if missing.any().any():
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame({'field':required, 'missing_exports':missing.sum().values})
    assignments = metadata.fillna('').copy()
    assignments[required] = assignments[required].apply(lambda c: c.astype(str).str.strip())
    if assignments.groupby('animal_id').running_group.nunique().gt(1).any():
        raise ValueError('An animal spans running groups; a longitudinal design needs explicit timepoints')
    keys = ['animal_id','section_id','anatomical_region','pla_type']
    if assignments.duplicated(keys).any():
        raise ValueError('Repeated animal/section/type analysis requires an explicit repeat policy')
    joined = summary.merge(assignments[['job',*required]], on='job', validate='one_to_one')
    animal = joined.groupby(['animal_id','running_group','anatomical_region','pla_type'], observed=True).agg(
        association_pct=('pla_association_pct','mean'), weighted_overlap_pct=('area_weighted_overlap_pct','mean'),
        sections=('section_id','nunique')).reset_index()
    groups = animal.groupby(['running_group','anatomical_region','pla_type'], observed=True).agg(
        animals=('animal_id','nunique'), mean_association_pct=('association_pct','mean'),
        min_association_pct=('association_pct','min'), max_association_pct=('association_pct','max')).reset_index()
    return animal, groups, pd.DataFrame()


# Preserve existing plotting entry points; implementation lives with the figures.
def run_overview_figure(summary):
    from .complete_plots import overview_figure
    return overview_figure(summary)


def normalized_coordinate_figure(objects):
    """Compatibility view: use explicit native extents instead of false alignment."""
    from .complete_plots import layout_figure
    return layout_figure(run_summary(objects))


def object_review_queue(objects, seed=7):
    """Reproducible review coordinates, combining random and targeted PLA rows."""
    frames=[]
    for job,g in objects.loc[objects.object_type=='PLA'].groupby('job',observed=True,sort=True):
        selections={
            'random baseline':g.sample(n=min(3,len(g)),random_state=seed),
            'largest reported area':g.nlargest(min(3,len(g)),'area'),
            'zero diameter or collapsed box':g.loc[(g.min_diameter==0)|(g.median_diameter==0)|(g.max_diameter==0)|(g.xmin==g.xmax)|(g.ymin==g.ymax)].head(2),
        }
        for reason,part in selections.items():
            frames.append(part[['job','image_id','region','object_id','area','overlap_pct','xmin','xmax','ymin','ymax','x','y']].assign(review_reason=reason))
    return pd.concat(frames,ignore_index=True)


def area_sensitivity_table(objects, thresholds=(0,.03,.1,.3,1)):
    rows=[]
    for job,g in objects.loc[objects.object_type=='PLA'].groupby('job',observed=True):
        for threshold in thresholds:
            kept=g.loc[g.area>=threshold]
            rows.append(dict(job=str(job),condition=str(g.condition.iloc[0]),pla_type=str(g.pla_type.iloc[0]),
                minimum_area_um2=threshold,retained_objects=len(kept),retained_pct=100*len(kept)/len(g),
                association_pct=100*kept.chat_present.mean(),
                weighted_overlap_pct=100*kept.overlap_area.sum()/kept.area.sum() if len(kept) else np.nan))
    return pd.DataFrame(rows)

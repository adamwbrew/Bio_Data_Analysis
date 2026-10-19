"""Audit condition versus exported image layout without anatomical registration.

Bounding-box envelopes describe where detections were stored in source-image
coordinates. They are neither tissue masks nor registered anatomical regions.
Within-image contrasts are descriptive provenance checks, not matched trials.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import StrMethodFormatter

from .handoff import ARMS, GROUPS

HOME = 'Home cage'
RUNNING = 'Running'
COLORS = {HOME: '#087E8B', RUNNING: '#D06738'}
INK = '#173447'
EXTENT_DEFINITION = 'Envelope of exported object bounding boxes; not a tissue ROI or anatomical registration'


def _validated_inputs(objects, metadata):
    meta_columns = ['job', 'animal_id', 'running_group', 'pla_arm', 'image_id']
    object_columns = ['job', 'xmin', 'xmax', 'ymin', 'ymax', 'object_type', 'chat_present']
    for frame, columns, label in ((metadata, meta_columns, 'metadata'),
                                   (objects, object_columns, 'objects')):
        missing = set(columns) - set(frame.columns)
        if missing:
            raise ValueError(f'{label} is missing required columns: {sorted(missing)}')
    m = metadata[meta_columns + (['image_label'] if 'image_label' in metadata else [])].copy()
    if m.empty or m[meta_columns].isna().any().any():
        raise ValueError('Missing section identity or condition assignment')
    for column in meta_columns:
        m[column] = m[column].astype(str)
        if m[column].str.strip().eq('').any():
            raise ValueError('Section identities and condition assignments must not be blank')
    if m.job.duplicated().any() or m[['animal_id', 'pla_arm']].duplicated().any():
        raise ValueError('Expected one unique section job per animal and PLA arm')
    if not m.pla_arm.isin(ARMS.values()).all():
        raise ValueError('Unknown PLA arm; the handoff arms must remain separate')
    if not m.animal_id.isin(GROUPS).all() or not m.running_group.eq(m.animal_id.map(GROUPS)).all():
        raise ValueError('Condition assignments differ from the scientist handoff')
    if set(m.running_group) != {HOME, RUNNING}:
        raise ValueError('Both Home cage and Running conditions are required for the layout audit')
    if 'image_label' not in m:
        labels = {key: f'Image {i + 1}' for i, key in enumerate(sorted(m.image_id.unique()))}
        m['image_label'] = m.image_id.map(labels)
    if m.image_label.isna().any() or m.groupby('image_id').image_label.nunique().gt(1).any():
        raise ValueError('Each source-image identity needs one nonmissing display label')
    if objects.job.isna().any() or set(objects.job.astype(str)) != set(m.job):
        raise ValueError('Object jobs must match the primary section metadata exactly')
    if objects.object_type.isna().any() or not objects.object_type.isin(['PLA', 'ChAT']).all():
        raise ValueError('Unknown or missing object type')
    if not objects.chat_present.dropna().isin([True, False, 0, 1]).all():
        raise ValueError('ChAT-presence values must be binary or missing')
    bounds = objects[['xmin', 'xmax', 'ymin', 'ymax']].to_numpy(dtype=float)
    if not np.isfinite(bounds).all():
        raise ValueError('Object bounding boxes must have finite numeric coordinates')
    if (bounds[:, 0] > bounds[:, 1]).any() or (bounds[:, 2] > bounds[:, 3]).any():
        raise ValueError('Object bounding-box lower bounds exceed upper bounds')
    return m


def _x_order(home, running):
    if len(home) == 0 or len(running) == 0:
        return 'incomplete conditions', np.nan
    if running.x_max.max() < home.x_min.min():
        return 'Running lower-X; Home cage higher-X', float(home.x_min.min() - running.x_max.max())
    if home.x_max.max() < running.x_min.min():
        return 'Home cage lower-X; Running higher-X', float(running.x_min.min() - home.x_max.max())
    return 'X extents overlap or touch', 0.


def layout_audit(objects, metadata):
    """Return section extents, within-image contrasts, and JSON-safe audit facts.

    Extents use the minima/maxima of the ORIGINAL object bounding boxes,
    including both object classes. Association is ChAT-positive PLA / all PLA.
    Missing flags make a section association unknown, rather than negative.
    No pair assignment, registration, density estimate, or test is inferred.
    """
    m = _validated_inputs(objects, metadata)
    keys = m.set_index('job')
    rows = []
    for job, frame in objects.groupby('job', observed=True, sort=True):
        job = str(job)
        pla = frame.loc[frame.object_type.eq('PLA')]
        valid = len(pla) > 0 and pla.chat_present.notna().all()
        positives = int(pla.chat_present.sum()) if valid else np.nan
        rows.append(dict(job=job, **keys.loc[job].to_dict(),
                         x_min=float(frame.xmin.min()), x_max=float(frame.xmax.max()),
                         y_min=float(frame.ymin.min()), y_max=float(frame.ymax.max()),
                         n_objects=len(frame), n_pla=len(pla), n_positive=positives,
                         association_pct=100 * positives / len(pla) if valid else np.nan))
    extent = pd.DataFrame(rows).sort_values(['pla_arm', 'image_id', 'running_group']).reset_index(drop=True)
    extent['x_mid'] = (extent.x_min + extent.x_max) / 2
    extent['y_mid'] = (extent.y_min + extent.y_max) / 2
    extent['x_width'] = extent.x_max - extent.x_min
    extent['y_height'] = extent.y_max - extent.y_min
    extent['extent_definition'] = EXTENT_DEFINITION
    contrasts = []
    for (arm, image), frame in extent.groupby(['pla_arm', 'image_id'], observed=True, sort=True):
        home, running = [frame.loc[frame.running_group.eq(group)] for group in (HOME, RUNNING)]
        h = home.association_pct.to_numpy(float)
        r = running.association_pct.to_numpy(float)
        complete = len(h) > 0 and len(r) > 0 and np.isfinite(np.r_[h, r]).all()
        order, gap = _x_order(home, running)
        contrasts.append(dict(pla_arm=arm, image_id=image, image_label=frame.image_label.iloc[0],
                              home_animals=', '.join(home.animal_id), running_animals=', '.join(running.animal_id),
                              n_home=len(home), n_running=len(running),
                              home_mean=float(h.mean()) if complete else np.nan,
                              running_mean=float(r.mean()) if complete else np.nan,
                              home_minus_running=float(h.mean() - r.mean()) if complete else np.nan,
                              x_order=order, x_gap_export_units=gap,
                              status='Descriptive within-source-image contrast; not a matched experimental pair'
                                     if complete else 'Incomplete condition or missing association; contrast unavailable'))
    contrasts = pd.DataFrame(contrasts)
    separated = contrasts.x_order.isin(['Running lower-X; Home cage higher-X',
                                       'Home cage lower-X; Running higher-X'])
    consistent_order = bool(separated.all() and contrasts.x_order.nunique() == 1)
    deltas = contrasts.home_minus_running.to_numpy(float)
    finite = np.isfinite(deltas)
    ranges = {}
    for group, frame in extent.groupby('running_group', observed=True):
        ranges[group] = dict(x_min=float(frame.x_min.min()), x_max=float(frame.x_max.max()))
    summary = dict(section_count=len(extent), source_image_count=int(extent.image_id.nunique()),
                   coordinate_units='HALO export units; calibration unverified',
                   extent_definition=EXTENT_DEFINITION,
                   complete_condition_x_separation=consistent_order,
                   condition_x_order=str(contrasts.x_order.iloc[0]) if consistent_order else 'Not consistently separated',
                   group_x_ranges=ranges,
                   within_image_comparisons=len(contrasts),
                   finite_within_image_contrasts=int(finite.sum()),
                   positive_within_image_contrasts=int((deltas[finite] > 0).sum()),
                   negative_within_image_contrasts=int((deltas[finite] < 0).sum()),
                   all_within_image_contrasts_positive=bool(finite.all() and (deltas > 0).all()),
                   interpretation='Condition follows the same left/right exported layout in every source-image/arm panel. '
                                  'This flags a technical-review priority; it does not establish an imaging artifact.'
                                  if consistent_order else
                                  'Condition does not follow one fully separated X layout across every source-image/arm panel.',
                   action='Review source-image/composite provenance, illumination and background correction, '
                          'HALO threshold uniformity, and equivalent blinded fields from both sides. '
                          'Shared image IDs do not establish matched experimental pairs or remove layout confounding.')
    return extent, contrasts, summary


def layout_figure(section_extent, metadata=None):
    """Show original X envelopes by source image and arm; no spatial alignment.

    ``metadata`` is accepted for call-site convenience; ``section_extent`` is
    self-contained and is the only source of plotted measurements.
    """
    required = {'job', 'animal_id', 'running_group', 'pla_arm', 'image_id', 'image_label', 'x_min', 'x_max'}
    if section_extent.empty or not required.issubset(section_extent):
        raise ValueError('A nonempty layout_audit section-extent table is required')
    arms = [arm for arm in ARMS.values() if arm in set(section_extent.pla_arm)]
    images = sorted(section_extent.image_id.unique())
    fig, axes = plt.subplots(len(images), len(arms), figsize=(14, 2.55 * len(images) + 2.05),
                             squeeze=False, sharex=True, facecolor='white')
    fig.subplots_adjust(left=.09, right=.98, top=.80, bottom=.16, hspace=.60, wspace=.25)
    lo, hi = float(section_extent.x_min.min()), float(section_extent.x_max.max())
    span = max(hi - lo, 1.)
    for i, image in enumerate(images):
        for j, arm in enumerate(arms):
            ax = axes[i, j]
            panel = section_extent.loc[section_extent.image_id.eq(image) & section_extent.pla_arm.eq(arm)]
            title = f'{panel.image_label.iloc[0] if len(panel) else image} · {arm}'
            ax.set_title(title, loc='left', color=INK, fontsize=12, fontweight='bold', pad=12)
            row_positions, row_labels = [], []
            ordered = pd.concat([panel.loc[panel.running_group.eq(RUNNING)],
                                 panel.loc[panel.running_group.eq(HOME)]], ignore_index=True)
            for k, row in enumerate(ordered.itertuples()):
                y = len(ordered) - 1 - k
                color = COLORS[row.running_group]
                ax.plot([row.x_min, row.x_max], [y, y], color=color, linewidth=12,
                        solid_capstyle='butt', alpha=.88)
                ax.scatter([row.x_min, row.x_max], [y, y], marker='|', s=170, color=color, linewidths=2)
                row_positions.append(y)
                row_labels.append(f'{row.animal_id}\n{row.running_group}')
                if np.isfinite(getattr(row, 'association_pct', np.nan)):
                    ax.text((row.x_min + row.x_max) / 2, y + .22,
                            f'{row.association_pct:.1f}% ChAT-positive PLA', color=color,
                            fontsize=9, ha='center', va='bottom')
            ax.set_yticks(row_positions, row_labels, fontsize=9)
            ax.set_ylim(-.42, max(len(ordered) - 1, 0) + .70)
            ax.set_xlim(lo - .025 * span, hi + .025 * span)
            ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
            ax.tick_params(axis='x', labelsize=9, colors=INK)
            ax.tick_params(axis='y', length=0, colors=INK)
            ax.grid(axis='x', color='#e5ecef', linewidth=.7)
            ax.set_axisbelow(True)
            for spine in ('top', 'right', 'left'):
                ax.spines[spine].set_visible(False)
            ax.spines['bottom'].set_color('#becdd3')
            if i == len(images) - 1:
                ax.set_xlabel('Original source-image X coordinate · export units', color=INK, labelpad=9)
    fig.text(.035, .965, 'Check image layout before interpreting the condition difference',
             fontsize=20, weight='bold', color=INK, ha='left', va='top')
    fig.text(.035, .918, 'Each bar spans the exported detections of one section. Coordinates have not been aligned or rescaled.',
             fontsize=11, color=INK, ha='left', va='top')
    fig.legend(handles=[Line2D([0], [0], color=COLORS[group], linewidth=7, label=group)
                        for group in (RUNNING, HOME)], loc='upper left', bbox_to_anchor=(.03, .886),
               ncol=2, frameon=False, fontsize=10)
    order_flags = []
    for _, frame in section_extent.groupby(['pla_arm', 'image_id'], observed=True):
        order_flags.append(_x_order(frame.loc[frame.running_group.eq(HOME)],
                                   frame.loc[frame.running_group.eq(RUNNING)])[0])
    consistent = len(set(order_flags)) == 1 and order_flags[0] not in ('X extents overlap or touch', 'incomplete conditions')
    conclusion = ('Condition and left/right layout coincide in every panel. This is a review priority, not proof of an artifact.'
                  if consistent else 'Compare condition with storage position in every panel; incomplete or overlapping layouts remain visible.')
    fig.text(.035, .085, conclusion, fontsize=11, weight='bold', color=INK, va='top')
    fig.text(.035, .048, 'Review image assembly, illumination/background correction and detection thresholds. '
             'Bar envelopes are not tissue masks; shared image IDs are not matched trials.',
             fontsize=9.5, color=INK, va='top')
    return fig

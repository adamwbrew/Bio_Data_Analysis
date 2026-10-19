"""Scientific displays for exploratory PLA morphology clustering.

Cluster labels are arm-specific learned partitions, never cell identities. Spatial
panels use native section coordinates and cannot establish anatomical group maps.
"""
from pathlib import Path
import json
import uuid

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.lines import Line2D

from .handoff import ANIMAL_ORDER, ARMS
from .handoff_plots import style, title, INK, TEAL, ORANGE, GROUP_COLORS

CLUSTER_COLORS = {1: '#3379AD', 2: '#AA477D', 3: '#CC982A', 4: '#529682'}
GROUPS = ('Home cage', 'Running')


def _frame(result, key):
    value = result.get(key)
    return value if isinstance(value, pd.DataFrame) else pd.DataFrame(value)


def _for_arm(frame, arm):
    return frame[frame.pla_arm == arm] if 'pla_arm' in frame else frame.iloc[:0]


def _arms(result):
    frame = _frame(result, 'animals')
    present = set(frame.get('pla_arm', []))
    return [a for a in ARMS.values() if a in present]


def _short(arm):
    return arm.replace('TDP43–', '')


def _footer(fig, text):
    fig.text(.03, -.025, text, color=INK, fontsize=9, ha='left', va='top')


def _empty(ax, message='No supported measurements'):
    ax.text(.5, .5, message, transform=ax.transAxes, ha='center', va='center', color=INK)
    ax.set_axis_off()


def _group_dots(ax, data, metric, ylabel):
    clusters = sorted(data.cluster.unique())
    order = {animal: i for i, animal in enumerate(ANIMAL_ORDER)}
    for i, cluster in enumerate(clusters):
        for gi, group in enumerate(GROUPS):
            rows = data[(data.cluster == cluster) & (data.running_group == group)].copy()
            rows['_order'] = rows.animal_id.map(order)
            rows = rows.sort_values('_order')
            finite = np.isfinite(rows[metric].to_numpy(float))
            vals = rows[metric].to_numpy(float)[finite]
            if not len(vals):
                continue
            base = i + (-.18 if gi == 0 else .18)
            offsets = np.linspace(-.065, .065, len(rows))[finite]
            ax.scatter(base + offsets, vals, s=38, color=GROUP_COLORS[group],
                       edgecolors='white', linewidths=.6, zorder=4)
            ax.plot([base-.11, base+.11], [vals.mean()]*2, color=GROUP_COLORS[group], lw=2.5)
    ax.set_xticks(range(len(clusters)), [f'C{c}' for c in clusters])
    ax.set_ylim(0, 105)
    ax.set_xlim(-.5, max(len(clusters)-.5, .5))
    ax.set_ylabel(ylabel)
    ax.set_xlabel('Learned morphology cluster · separate model per arm')
    ax.grid(axis='y', color='#e8eef1')


def overview_figure(result):
    """Feature space, per-animal mixture, and held-out-from-fit association."""
    style()
    arms = _arms(result)
    fig, axes = plt.subplots(max(len(arms), 1), 3, figsize=(17, max(1, len(arms))*4.1),
                             layout='constrained', squeeze=False)
    scatter, profiles, animals = (_frame(result, key) for key in ('scatter', 'profiles', 'animals'))
    if not arms:
        for ax in axes.flat:
            _empty(ax)
    for i, arm in enumerate(arms):
        sampled = _for_arm(scatter, arm)
        data = animals[animals.pla_arm == arm]
        ax = axes[i, 0]
        for c in sorted(data.cluster.unique()):
            points = sampled[sampled.cluster == c] if len(sampled) else sampled
            if len(points):
                ax.scatter(points.area, points['shape'], s=5, alpha=.2,
                           color=CLUSTER_COLORS.get(c, '#7d8b92'), rasterized=True)
            arm_profiles = _for_arm(profiles, arm)
            profile = arm_profiles[arm_profiles.cluster == c] if len(arm_profiles) else arm_profiles
            if len(profile):
                row = profile.iloc[0]
                if row.median_area > 0 and np.isfinite(row.median_shape):
                    ax.scatter(row.median_area, row.median_shape, s=80,
                               color=CLUSTER_COLORS.get(c, '#7d8b92'), edgecolors='white', linewidths=1.5)
                    offset = (5, 7) if c == 1 else ((5, 22) if c == 2 else (5, -32))
                    ax.annotate(f'C{c}\n{row.median_area:.3g} µm² · {row.median_shape:.2f}',
                                (row.median_area, row.median_shape), xytext=offset,
                                textcoords='offset points', fontsize=8, ha='left',
                                bbox=dict(boxstyle='round,pad=.2', facecolor='white', edgecolor='none', alpha=.85))
        ax.set_xscale('log')
        ax.set_ylim(-.03, 1.07)
        ax.set_ylabel('Shape ratio · minimum / maximum diameter')
        ax.set_xlabel('PLA area (µm²) · log scale')
        ax.set_title(f'{_short(arm)} | two measured features', loc='left')
        _group_dots(axes[i, 1], data, 'proportion_pct', 'Share of eligible PLA objects (%)')
        axes[i, 1].set_title('Does the morphology mix differ?', loc='left')
        _group_dots(axes[i, 2], data, 'association_pct', 'ChAT-positive PLA within cluster (%)')
        axes[i, 2].set_title('Does association differ within a cluster?', loc='left')
    handles = [Line2D([0], [0], marker='o', color=GROUP_COLORS[g], lw=2,
                      markersize=5, label=g) for g in GROUPS]
    axes[0, 1].legend(handles=handles, loc='upper right', ncol=2, frameon=False, fontsize=9)
    title(fig, 'Machine learning describes PLA morphology; animals remain the replicates',
          'Unsupervised k=3 models use area and shape only. ChAT association and running labels are inspected afterward; neither is used to fit the clusters.')
    audit = _frame(result, 'audit')
    excluded = f"{int(audit.n_excluded.sum()):,}/{int(audit.n_pla.sum()):,} PLA objects excluded from morphology modeling. " if len(audit) else ''
    _footer(fig, 'Left: sampled objects; labels show full eligible-object medians. C1–C3 are area-ranked within each arm, not cell identities. A zero diameter can reflect measurement precision.\n'
                + excluded + 'Middle/right: animal dots and equal-animal means (n=3/group); primary unfiltered results remain unchanged.')
    return fig


def robustness_figure(result):
    """Show cluster-count sensitivity and omitted-animal refit agreement."""
    style()
    arms = _arms(result)
    fig, axes = plt.subplots(max(len(arms), 1), 2, figsize=(14, max(len(arms), 1)*4),
                             layout='constrained', squeeze=False)
    standardized, stability, area_only = (_frame(result, key) for key in ('standardized', 'stability', 'area_only_standardized'))
    for i, arm in enumerate(arms):
        data = _for_arm(standardized, arm)
        if len(data):
            data = data.sort_values('k')
        ax = axes[i, 0]
        if len(data):
            ax.axhline(0, color='#92a3ab', lw=1)
            a = area_only[area_only.pla_arm == arm] if len(area_only) else area_only
            values = [data.raw_difference_pp.iloc[0], data.valid_difference_pp.iloc[0],
                      a.standardized_difference_pp.iloc[0] if len(a) else np.nan] + list(data.standardized_difference_pp)
            labels = ['All PLA\nobserved', 'Eligible PLA\nobserved', 'Same area\nmix · k=3'] + [f'Morphology\nmix · k={k}' for k in data.k]
            colors = ['#445c6a', '#8799a3', '#805AA4'] + [TEAL]*len(data)
            ax.plot(range(3, len(values)), values[3:], color=TEAL, lw=2)
            ax.scatter(range(len(values)), values, s=65, c=colors, zorder=3)
            for j, value in enumerate(values):
                if np.isfinite(value):
                    ax.annotate(f'{value:+.2f}', (j, value), xytext=(0, 8),
                                textcoords='offset points', ha='center', fontsize=9)
            values = np.asarray(values)
            finite = values[np.isfinite(values)]
            if len(finite):
                ax.set_ylim(min(0, finite.min())-3, max(0, finite.max())+7)
            ax.set_xticks(range(len(values)), labels, fontsize=8)
            ax.set_xlim(-.5, len(values)-.5)
        else:
            _empty(ax)
        ax.set_title(f'{_short(arm)} | does the morphology mix explain the gap?', loc='left')
        ax.set_ylabel('Home cage − Running · ChAT association (pp)')
        ax.set_xlabel('Area-only analysis retains detections with unusable diameters', fontsize=9)
        ax.grid(axis='y', color='#e8eef1')
        ax = axes[i, 1]
        data = _for_arm(stability, arm)
        if len(data):
            data = data.assign(_order=data.animal_id.map({a: i for i, a in enumerate(ANIMAL_ORDER)})).sort_values('_order')
            for j, row in enumerate(data.itertuples()):
                group = 'Home cage' if row.animal_id in ANIMAL_ORDER[:3] else 'Running'
                ax.barh(j, row.agreement_pct, height=.45, color=GROUP_COLORS[group], alpha=.7)
                ax.text(101, j, f'{row.agreement_pct:.1f}%', va='center', fontsize=9)
            ax.set_yticks(range(len(data)), [f'{r.animal_id} · ARI {r.ari:.3f}' for r in data.itertuples()])
            ax.invert_yaxis()
            ax.set_xlim(0, 112)
            ax.set_xlabel('Matched cluster assignments agreeing (%)')
            ax.set_ylabel('Omitted animal · adjusted Rand index (ARI)')
            ax.grid(axis='x', color='#e8eef1')
        else:
            _empty(ax)
        ax.set_title('Refit without the named animal · k=3', loc='left')
    title(fig, 'Two checks on the exploratory result',
          'Left: compare all PLA, morphology-eligible PLA, and common mixtures. Right: refit preprocessing and clustering on five animals, then compare the omitted animal.')
    _footer(fig, 'Common weights remove differences in the observed morphology mixture only; they do not correct tissue coverage, ChAT availability, imaging or segmentation.\n'
                'Agreement compares an omitted-animal refit with the full-data reference; the reference includes that animal. This is stability, not predictive accuracy or independent validation.')
    return fig


def spatial_figure(result, arm):
    """Native per-animal atlas for all three morphology fractions, unsmoothed."""
    style()
    sections = [section for section in result.get('spatial', []) if section['pla_arm'] == arm]
    order = {a: i for i, a in enumerate(ANIMAL_ORDER)}
    sections.sort(key=lambda section: order.get(section['animal_id'], 100))
    clusters = sorted({int(c) for section in sections for c in section['cluster_fractions']})
    nr, nc = max(len(sections), 1), max(len(clusters), 1)
    fig, axes = plt.subplots(nr, nc, figsize=(max(5.2*nc, 6), max(3.8*nr, 5)),
                             layout='constrained', squeeze=False)
    image = None
    for row, section in enumerate(sections):
        xmin, xmax, ymin, ymax = section['extent']
        extent = [xmin, xmax, ymax, ymin]
        counts = np.asarray(section['n'], dtype=float)
        all_counts = np.asarray(section.get('n_all', section['n']), dtype=float)
        for col, c in enumerate(clusters):
            ax = axes[row, col]
            ax.imshow(np.where(all_counts > 0, 1., np.nan), extent=extent, origin='upper',
                      cmap=ListedColormap(['#e5ebee']), interpolation='nearest', aspect='equal')
            values = np.asarray(section['cluster_fractions'].get(str(c), np.full(counts.shape, np.nan)), dtype=float)
            # Enforce the display denominator even if a caller supplies unmasked values.
            values = np.where(counts >= section.get('min_count', 10), values, np.nan)
            image = ax.imshow(np.ma.masked_invalid(values), extent=extent, origin='upper',
                              cmap='viridis', norm=Normalize(0, 100), interpolation='nearest', aspect='equal')
            ax.set_title(f"{section['animal_id']} · {section['running_group']} | C{c}",
                         color=GROUP_COLORS.get(section['running_group'], INK), loc='left', fontsize=11)
            ax.set_xlabel('x · native export units', fontsize=9)
            ax.set_ylabel('y · native export units', fontsize=9)
            ax.ticklabel_format(style='sci', scilimits=(4, 4), axis='both')
            ax.tick_params(labelsize=8)
    if image is not None:
        fig.colorbar(image, ax=axes.ravel().tolist(), shrink=.45, pad=.012,
                     label='Share of eligible PLA objects belonging to the cluster (%)')
    else:
        for ax in axes.flat:
            _empty(ax)
    title(fig, f'{_short(arm)} | where are the learned morphology patterns?',
          'Each row is a separate animal/section. Columns are arm-specific C1–C3, ordered by fitted area; identical 0–100% color scales throughout.')
    _footer(fig, 'Native 500 × 500 export-unit bins; ≥10 eligible PLA objects per bin. No smoothing, cross-animal overlay or anatomical registration.\n'
                'Pale gray: PLA observed but too few eligible objects for a fraction; white: no PLA support. These are object-derived maps, not microscopy images or tissue masks.')
    return fig


def _json_clean(value):
    """Strict JSON: absent/nonfinite estimates stay missing, never become zeros."""
    if isinstance(value, pd.DataFrame):
        return _json_clean(value.to_dict(orient='records'))
    if isinstance(value, dict):
        return {str(key): _json_clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_json_clean(item) for item in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if value is pd.NA or value is None:
        return None
    return value


def explorer_payload(result):
    """Expose only interpretable summaries and sampled native coordinates."""
    keys = ('profiles', 'animals', 'contrasts', 'standardized', 'area_only_standardized', 'stability',
            'scatter', 'spatial', 'audit', 'review_targets', 'definitions')
    payload = {key: result.get(key, []) for key in keys}
    payload.update(arms=_arms(result), animal_order=ANIMAL_ORDER,
                   cluster_colors=CLUSTER_COLORS, group_colors=GROUP_COLORS)
    return _json_clean(payload)


def explorer_html(result, include_plotly=True):
    """Self-contained, offline Plotly view for notebooks and the linked app."""
    from plotly.offline import get_plotlyjs
    payload = json.dumps(explorer_payload(result), allow_nan=False, ensure_ascii=False,
                         separators=(',', ':')).replace('<', '\\u003c')
    template = Path(__file__).with_name('morphology_explorer.html').read_text()
    return (template.replace('__ROOT_ID__', 'morphology-' + uuid.uuid4().hex)
            .replace('__DATA__', payload)
            .replace('__PLOTLY__', get_plotlyjs() if include_plotly else ''))

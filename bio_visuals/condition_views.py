"""Scientific condition contrasts: equal-animal averages and explicit subtraction.

These displays compare measured section summaries. They never align anatomy or
turn unreviewed coordinates into a condition-effect map.
"""
from pathlib import Path
import json
import uuid

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from .handoff import ARMS, ANIMAL_ORDER, SIZE_LABELS
from .handoff_plots import style, title, INK, TEAL, ORANGE, GROUP_COLORS

DIFF_CMAP = 'BrBG'  # orange/brown = Running higher; teal = Home cage higher
GROUPS = ('Home cage', 'Running')


def _definitions(results):
    return results['definitions']['metrics']


def _short_arm(arm):
    return arm.replace('TDP43–', '')


def _footer(fig, text):
    fig.text(.03, -.025, text, fontsize=9, ha='left', va='top', color=INK)


def _format(value, metric=None):
    if not np.isfinite(value):
        return 'not estimable'
    if metric in ('pla_count', 'chat_count'):
        return f'{value:,.0f}'
    return f'{value:.2f}' if abs(value) < 1 else f'{value:.1f}'


def _dot_means(ax, data, metric, label):
    for x, group in enumerate(GROUPS):
        rows = data[data.running_group == group].copy()
        rows['_order'] = rows.animal_id.map({a:i for i,a in enumerate(ANIMAL_ORDER)})
        rows = rows.sort_values('_order')
        offsets = np.linspace(-.14, .14, len(rows))
        values = rows[metric].to_numpy(float)
        ax.scatter(x + offsets, values, s=65, color=GROUP_COLORS[group],
                   edgecolors='white', linewidths=.8, zorder=3)
        mean = values.mean()
        ax.plot([x-.23, x+.23], [mean, mean], color=GROUP_COLORS[group], lw=3, zorder=2)
        for dx, row in zip(offsets, rows.itertuples()):
            ax.annotate(row.animal_id, (x+dx, getattr(row, metric)), xytext=(4, 4),
                        textcoords='offset points', fontsize=8)
        ax.text(x, 1.02, f'Mean {_format(mean, metric)}', transform=ax.get_xaxis_transform(),
                color=GROUP_COLORS[group], ha='center', fontsize=9, weight='bold')
    ax.set_xlim(-.5, 1.55)
    ax.set_xticks([0, 1], ['Home cage\nnon-running', 'Running'])
    ax.set_ylabel(label)
    ax.grid(axis='y', color='#e8eef1', zorder=0)
    if metric.endswith('_pct'):
        ax.set_ylim(0, 104)
    else:
        ax.margins(y=.2)


def _contrast(ax, row, unit, limit=None):
    delta = float(row.home_minus_running)
    low, high = float(row.leave_pair_out_min), float(row.leave_pair_out_max)
    scale = max(abs(v) for v in [delta, low, high, 1e-6]) if limit is None else limit
    ax.axvline(0, color='#859ba5', lw=1.2)
    ax.barh(0, delta, height=.28, color=TEAL if delta >= 0 else ORANGE, alpha=.85)
    ax.plot([low, high], [0, 0], color=INK, lw=2, zorder=3)
    ax.plot([low, high], [0, 0], '|', color=INK, ms=16, mew=1.5, zorder=3)
    ax.plot(delta, 0, 'o', color='white', markeredgecolor=INK, ms=7, zorder=4)
    ax.set_xlim(-scale*1.45, scale*1.45)
    ax.set_ylim(-.5, .65)
    ax.set_yticks([])
    ax.set_xlabel('Home cage − Running\n' + unit)
    ax.set_title(f'{delta:+.2f} {unit}', loc='left', color=TEAL if delta >= 0 else ORANGE)
    ax.text(.5, .03, f'Delete-one-per-group range\n{low:+.2f} to {high:+.2f} {unit}',
            transform=ax.transAxes, ha='center', fontsize=8)
    ax.grid(axis='x', color='#e8eef1', zorder=0)


def average_difference_figure(results):
    """Animal points and means, beside signed contrasts and sensitivity ranges."""
    style()
    defs = _definitions(results)
    metrics = ['association_pct', 'overlap_area_ratio_pct']
    fig, axes = plt.subplots(2, 4, figsize=(17, 8.3), layout='constrained',
                             gridspec_kw={'width_ratios':[1.3, 1, 1.3, 1]})
    animals, contrasts = results['animal_metrics'], results['contrasts']
    for i, arm in enumerate(ARMS.values()):
        data = animals[animals.pla_arm == arm]
        for j, metric in enumerate(metrics):
            meta = defs[metric]
            _dot_means(axes[i, j*2], data, metric, meta['label'] + (' (%)' if metric.endswith('_pct') else ''))
            axes[i, j*2].set_title(f'{_short_arm(arm)} | animal values', loc='left', pad=30)
            row = contrasts[(contrasts.pla_arm == arm) & (contrasts.metric == metric)].iloc[0]
            same = contrasts[contrasts.metric == metric]
            limit = max(abs(same.leave_pair_out_min).max(), abs(same.leave_pair_out_max).max(),
                        abs(same.home_minus_running).max(), 1e-6)
            _contrast(axes[i, j*2+1], row, 'pp' if meta['unit'] == 'percentage points' else meta['unit'], limit)
    title(fig, 'What changes when Running is subtracted from non-running?',
          'Every dot is one animal. Horizontal colored bars are equal-animal means (n=3/group/arm). Positive differences mean Home cage is higher.')
    _footer(fig, 'Black intervals delete one animal from each group (2 versus 2, nine combinations); they are sensitivity ranges, not confidence intervals.\n'
                'ChAT-positive PLA is a binary object association. Overlap-area ratio uses exported geometric overlap area / total PLA area; neither establishes molecular binding.')
    return fig


def _matrix(table, row_key, col_key, value='home_minus_running'):
    rows = [a for a in ANIMAL_ORDER if a in set(table[row_key])]
    cols = [a for a in ANIMAL_ORDER if a in set(table[col_key])]
    return table.pivot(index=row_key, columns=col_key, values=value).reindex(index=rows, columns=cols), rows, cols


def robustness_figure(results):
    """Cross-animal contrasts and delete-one-per-group checks for association."""
    style()
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout='constrained')
    primary = results['definitions']['primary_metric']
    cross, leave = results['cross_animal'], results['leave_pair_out']
    leave = leave[leave.metric == primary]
    maximum = max(abs(cross.home_minus_running).max(), abs(leave.home_minus_running).max(), 1)
    norm = TwoSlopeNorm(vmin=-maximum, vcenter=0, vmax=maximum)
    for i, arm in enumerate(ARMS.values()):
        panels = [(cross[cross.pla_arm == arm], 'home_animal', 'running_animal',
                   'All 9 between-animal differences', 'Running animal', 'Home cage animal'),
                  (leave[leave.pla_arm == arm], 'omitted_home', 'omitted_running',
                   'Omit one per group', 'Omitted Running animal', 'Omitted Home cage animal')]
        for ax, (table, rk, ck, heading, xlabel, ylabel) in zip(axes[i], panels):
            matrix, rows, cols = _matrix(table, rk, ck)
            im = ax.imshow(matrix, cmap=DIFF_CMAP, norm=norm, aspect='equal')
            for y in range(len(rows)):
                for x in range(len(cols)):
                    value = matrix.iloc[y, x]
                    ax.text(x, y, f'{value:+.1f}', ha='center', va='center', fontsize=15,
                            color='white' if abs(value) > maximum*.55 else INK, weight='bold')
            ax.set_xticks(range(len(cols)), cols)
            ax.set_yticks(range(len(rows)), rows)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(f'{_short_arm(arm)} | {heading}', loc='left', pad=12)
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=.75, pad=.02,
                 label='Home cage − Running · ChAT-positive PLA (percentage points)')
    title(fig, 'Is the average difference driven by one animal?',
          'Left: each Home cage animal minus each Running animal. Right: group-mean difference after omitting the named animal from each group.')
    _footer(fig, 'The 9 cells are dependent descriptive comparisons, not 9 independent replicates or matched pairs. Actual replication remains n=3 animals/group/arm.\n'
                'Shared color limits allow comparison across panels. These checks do not remove confounding or establish a causal running effect.')
    return fig


def size_standardization_figure(results):
    """Show size-stratum associations and the controlled composition comparison."""
    style()
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.8), layout='constrained',
                             gridspec_kw={'width_ratios':[1.25, 1, 1.2]})
    sd, sg = results['size_animals'], results['size_groups']
    sc, standard = results['size_contrasts'], results['standardization']
    labels = {'raw':'All exported PLA', 'restricted':'Supported size bins', 'standardized':'Same size mix'}
    for i, arm in enumerate(ARMS.values()):
        ax, dx, sx = axes[i]
        for (_, animal), data in sd[sd.pla_arm == arm].groupby(['running_group', 'animal_id'], observed=True):
            data = data.sort_values('size_index')
            group = data.running_group.iloc[0]
            ax.plot(data.size_index, data.fraction_pct, color=GROUP_COLORS[group], alpha=.24, lw=1)
        for group in GROUPS:
            data = sg[(sg.pla_arm == arm) & (sg.running_group == group)].sort_values('size_index')
            ax.plot(data.size_index, data.mean_fraction_pct, 'o-', color=GROUP_COLORS[group],
                    lw=2.7, ms=5, label=group)
        ax.set_ylim(0, 100)
        ax.set_xticks(range(len(SIZE_LABELS)), SIZE_LABELS, rotation=20)
        ax.set_xlabel('PLA area (µm²) · fixed descriptive bins')
        ax.set_ylabel('ChAT-positive PLA (%)')
        ax.set_title(f'{_short_arm(arm)} | size-stratum averages', loc='left')
        ax.grid(axis='y', color='#e8eef1')
        if i == 0:
            ax.legend(frameon=False, loc='lower right', fontsize=9)
        rows = sc[sc.pla_arm == arm].sort_values('size_index')
        colors = [TEAL if v >= 0 else ORANGE for v in rows.home_minus_running]
        bars = dx.barh(rows.size_index, rows.home_minus_running, color=colors, alpha=.85)
        for bar, row in zip(bars, rows.itertuples()):
            if not row.eligible:
                bar.set_facecolor('#b8c5cc');bar.set_hatch('//')
            value = row.home_minus_running
            if np.isfinite(value):
                dx.text(value + (1 if value >= 0 else -1), row.size_index, f'{value:+.1f}',
                        va='center', ha='left' if value >= 0 else 'right', fontsize=9)
        dx.axvline(0, color='#859ba5', lw=1)
        dx.set_yticks(range(len(SIZE_LABELS)), SIZE_LABELS)
        dx.invert_yaxis()
        dx.margins(x=.25)
        dx.set_xlabel('Home cage − Running (pp)')
        dx.set_title('Difference within each size bin', loc='left')
        data = standard[standard.pla_arm == arm].set_index('estimator').reindex(labels)
        values = data.home_minus_running.to_numpy(float)
        sx.axvline(0, color='#859ba5', lw=1)
        sx.barh(range(3), values, color=[TEAL if v >= 0 else ORANGE for v in values], alpha=.85, height=.52)
        for j, (_, row) in enumerate(data.iterrows()):
            lo, hi = row.leave_pair_out_min, row.leave_pair_out_max
            if np.isfinite([lo, hi]).all():
                sx.plot([lo, hi], [j, j], color=INK, lw=1.5)
                sx.plot([lo, hi], [j, j], '|', color=INK, ms=12)
            if np.isfinite(row.home_minus_running):
                sx.text(row.home_minus_running, j-.32, f'{row.home_minus_running:+.2f} pp',
                        ha='center', fontsize=9, weight='bold')
        sx.set_yticks(range(3), list(labels.values()))
        sx.invert_yaxis()
        sx.margins(x=.3, y=.25)
        sx.set_xlabel('Home cage − Running (pp)')
        sx.set_title('Does matching the size mix change it?', loc='left')
    title(fig, 'Does a different mix of small and large PLA objects explain the contrast?',
          'Thin lines: individual animals. Thick lines: equal-animal means. Common size weights are applied to every animal, within each separate PLA arm.')
    minimum = results['definitions'].get('min_per_stratum', 10)
    _footer(fig, f'Only size bins with ≥{minimum} PLA objects in every animal enter the common-size comparison; hatched gray bars are unsupported bins.\n'
                'Intervals are delete-one-per-group sensitivity ranges, not confidence intervals. Standardization is descriptive; delete-one checks hold common size weights fixed.')
    return fig


def comparison_html(results, registration=None, include_plotly=True):
    """Return a portable, offline interactive condition-comparison HTML fragment.

    ``registration`` is returned by ``condition_registration``. Supplied maps
    render only when review is complete; uncovered cells remain unknown.
    ``include_plotly=False`` assumes the embedding page already loads Plotly.
    """
    from plotly.offline import get_plotlyjs
    from .condition_comparison import comparison_payload
    payload = comparison_payload(results)
    payload['registration'] = registration or {'ready':False, 'approved':0, 'total':12,
        'message':'Anatomical subtraction is pending reviewed alignment, common coordinate scale and sampling support.'}
    serialized = json.dumps(payload, allow_nan=False, ensure_ascii=False, separators=(',', ':')).replace('<', '\\u003c')
    template = Path(__file__).with_name('condition_comparison.html').read_text(encoding='utf-8')
    return template.replace('__ROOT_ID__', 'conditions-'+uuid.uuid4().hex).replace('__DATA__', serialized).replace('__PLOTLY__', get_plotlyjs() if include_plotly else '')

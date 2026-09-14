"""Distribution, denominator, and sensitivity figures using every exported row.

All tables returned by these functions contain the actual plotted values.
No p-values or animal-level uncertainty are computed from object rows.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import PercentFormatter, FuncFormatter
from .data import require_single_roi


COLORS = {"4594": "#0072B2", "4595": "#D55E00", "4596": "#009E73", "4597": "#CC79A7"}
STYLES = {"4594": "-", "4595": "--", "4596": "-", "4597": "--"}
COLORS.update({"4999": COLORS["4594"], "5000": COLORS["4595"], "5001": COLORS["4596"], "5002": COLORS["4597"]})
STYLES.update({"4999": "-", "5000": "--", "5001": "-", "5002": "--"})


def setup_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 11, "axes.labelsize": 10,
        "figure.titlesize": 16, "figure.dpi": 140,
        "savefig.dpi": 300, "axes.spines.top": False,
        "axes.spines.right": False, "axes.axisbelow": True,
        "pdf.fonttype": 42, "svg.fonttype": "none",
        "figure.facecolor": "white", "axes.facecolor": "white",
    })


def _groups(objects, object_type=None):
    require_single_roi(objects)
    data = objects if object_type is None else objects[objects.object_type == object_type]
    return data.groupby("job", observed=True, sort=True)


def _label(group):
    column = "pair_id" if "pair_id" in group else "region"
    return str(group[column].iloc[0]).replace("_", " · ").replace(" PLA", "")


def _footer(fig, text):
    fig.text(0.01, 0.015, text, ha="left", va="bottom", fontsize=9, color="#444444")


def overview(objects):
    """Compare counts and explicitly different overlap denominators."""
    rows = []
    for job, group in _groups(objects):
        pla = group[group.object_type == "PLA"]
        chat = group[group.object_type == "ChAT"]
        rows.append({"job": str(job), "region": str(group.region.iloc[0]),
                     "label": _label(group), "pla_objects": len(pla), "chat_objects": len(chat),
                     "pla_with_chat": int(pla.chat_present.sum()),
                     "pla_without_chat": int((~pla.chat_present).sum()),
                     "pla_object_fraction_pct": 100 * pla.chat_present.mean(),
                     "pla_mean_overlap_pct": pla.overlap_pct.mean(),
                     "pla_area_weighted_overlap_pct": 100 * pla.overlap_area.sum() / pla.area.sum()})
    table = pd.DataFrame(rows)
    y = np.arange(len(table))
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.8), gridspec_kw={"width_ratios": [1.2, 1, 1.2]})
    ax = axes[0]
    ax.barh(y, table.chat_objects, color="#c5c8cc", label="ChAT objects")
    ax.barh(y, table.pla_objects, left=table.chat_objects, color="#0072B2", label="PLA objects")
    for i, row in table.iterrows():
        total = row.pla_objects + row.chat_objects
        ax.text(total + 4500, i, f"{total:,}", va="center", fontsize=9)
    ax.set_xlim(0, 440000)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x/1000:g}k"))
    ax.set_yticks(y)
    ax.set_yticklabels(table.label)
    ax.set_title("A  What was detected?")
    ax.set_xlabel("Exported objects")
    ax.legend(loc="lower right", fontsize=8)

    ax = axes[1]
    ax.barh(y, table.pla_without_chat, color="#c5c8cc", label="ChAT absent")
    ax.barh(y, table.pla_with_chat, left=table.pla_without_chat, color="#D55E00", label="ChAT present")
    for i, row in table.iterrows():
        ax.text(row.pla_objects + 450, i, f"{row.pla_objects:,}", va="center", fontsize=9)
    ax.set_xlim(0, max(table.pla_objects) * 1.22)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x/1000:g}k"))
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_title("B  PLA objects only")
    ax.set_xlabel("PLA objects")
    ax.legend(loc="lower right", fontsize=8)

    ax = axes[2]
    specs = [("pla_object_fraction_pct", -0.20, "o", "#0072B2", "Objects with any overlap"),
             ("pla_mean_overlap_pct", 0, "s", "#D55E00", "Mean object area overlap"),
             ("pla_area_weighted_overlap_pct", 0.20, "^", "#009E73", "Pooled object area overlap")]
    for field, offset, marker, color, label in specs:
        ax.scatter(table[field], y + offset, marker=marker, color=color, s=45, label=label, zorder=3)
    ax.set_xlim(0, 100)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.xaxis.set_major_formatter(PercentFormatter())
    ax.set_title("C  The denominator changes the result")
    ax.set_xlabel("Overlap measure among PLA objects")
    ax.legend(loc="lower right", fontsize=8)
    for ax in axes:
        ax.invert_yaxis()
        ax.set_ylim(3.9, -0.6)
        ax.grid(axis="x", alpha=0.17)
    fig.suptitle("Four exports: object abundance and PLA–ChAT overlap", x=0.01, ha="left")
    fig.subplots_adjust(left=0.14, right=0.98, top=0.83, bottom=0.23, wspace=0.20)
    _footer(fig, "Every exported object retained. Counts are not normalized by tissue area. These are four region exports; biological replication is unknown.")
    return fig, table


def ecdf_figures(objects):
    """Exact tie-aware ECDFs; return all unique plotted values and counts."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    frames = []
    for row_i, (metric, axis_label) in enumerate([
        ("area", "Object area (µm²; log scale)"),
        ("intensity", "Average intensity (export units; log scale)")]):
        for col_i, object_type in enumerate(["PLA", "ChAT"]):
            ax = axes[row_i, col_i]
            for job, group in _groups(objects, object_type):
                values, counts = np.unique(group[metric].to_numpy(), return_counts=True)
                cumulative = counts.cumsum() / counts.sum() * 100
                ax.step(values, cumulative, where="post", color=COLORS.get(str(job)),
                        linestyle=STYLES.get(str(job), "-"), label=_label(group), linewidth=1.6)
                frames.append(pd.DataFrame({"job": str(job), "region": str(group.region.iloc[0]),
                                            "object_type": object_type, "metric": metric,
                                            "value": values, "object_count_at_value": counts,
                                            "cumulative_pct": cumulative}))
            ax.set_xscale("log")
            ax.set_ylim(0, 100)
            ax.yaxis.set_major_formatter(PercentFormatter())
            ax.set_xlabel(axis_label)
            ax.set_ylabel("Objects at or below this value")
            ax.set_title(f"{object_type} objects — {'area' if metric == 'area' else 'intensity'}")
            ax.grid(alpha=0.18)
            if row_i == 0 and col_i == 0:
                ax.legend(loc="lower right", fontsize=8)
    fig.suptitle("Size and brightness distributions: read the full range", x=0.01, ha="left")
    fig.subplots_adjust(left=0.09, right=0.98, top=0.91, bottom=0.14, wspace=0.25, hspace=0.34)
    _footer(fig, "Exact cumulative distributions of all objects; no sampling or smoothing. Log X axes reveal long tails. Intensity scale and channel comparability are unconfirmed.")
    return fig, pd.concat(frames, ignore_index=True)


def overlap_composition(objects):
    """Zero/partial/full masses and full-data overlap distribution, PLA only."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
    states = []
    for i, (job, group) in enumerate(_groups(objects, "PLA")):
        vals = group.overlap_pct.to_numpy()
        counts = [(vals == 0).sum(), ((vals > 0) & (vals < 100)).sum(), (vals == 100).sum()]
        left = 0
        for state, count, color in zip(["No overlap", "Partial overlap", "Full overlap"], counts,
                                       ["#c5c8cc", "#56B4E9", "#0072B2"]):
            pct = count / len(vals) * 100
            axes[0].barh(i, pct, left=left, color=color, label=state if i == 0 else None)
            if pct >= 8:
                axes[0].text(left + pct / 2, i, f"{pct:.1f}%", ha="center", va="center",
                             color="white" if state == "Full overlap" else "#111111", fontsize=9)
            states.append({"job": str(job), "region": str(group.region.iloc[0]), "state": state,
                           "count": int(count), "pct_of_pla_objects": pct})
            left += pct
        values, cnt = np.unique(vals, return_counts=True)
        axes[1].step(values, 100 * cnt.cumsum() / cnt.sum(), where="post", linewidth=1.8,
                     color=COLORS.get(str(job)), linestyle=STYLES.get(str(job), "-"), label=_label(group))
    axes[0].set_yticks(range(4))
    axes[0].set_yticklabels([_label(g) for _, g in _groups(objects, "PLA")])
    axes[0].invert_yaxis()
    axes[0].set_xlim(0, 100)
    axes[0].xaxis.set_major_formatter(PercentFormatter())
    axes[0].set_title("A  Three distinct overlap states")
    axes[0].set_xlabel("Percentage of PLA objects")
    axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.21), ncol=3, fontsize=8)
    axes[1].set(xlim=(-1, 101), ylim=(0, 101), xlabel="PLA object area overlapping ChAT (%)",
                ylabel="PLA objects at or below this overlap (%)", title="B  Exact overlap distribution")
    axes[1].legend(loc="lower right", fontsize=8)
    axes[1].grid(alpha=0.18)
    fig.suptitle("Some PLA objects have no overlap; others are fully overlapped", x=0.01, ha="left")
    fig.subplots_adjust(left=0.17, right=0.98, top=0.83, bottom=0.25, wspace=0.30)
    _footer(fig, "PLA objects only. Endpoint masses at 0% and 100% are preserved; a smooth density curve could obscure them.")
    return fig, pd.DataFrame(states)


def size_overlap(objects, bins=36):
    """Full-count 2D histograms; one common area grid and color scale."""
    pla = objects[objects.object_type == "PLA"]
    xedges = np.geomspace(pla.area.min(), pla.area.max() * (1 + 1e-10), bins + 1)
    yedges = np.linspace(0, 100, 26)
    histograms = []
    rows = []
    for job, group in _groups(objects, "PLA"):
        hist, _, _ = np.histogram2d(group.area, group.overlap_pct, bins=(xedges, yedges))
        histograms.append((str(job), group, hist))
        for i in range(hist.shape[0]):
            for j in range(hist.shape[1]):
                rows.append({"job": str(job), "area_left_um2": xedges[i], "area_right_um2": xedges[i+1],
                             "overlap_left_pct": yedges[j], "overlap_right_pct": yedges[j+1],
                             "pla_count": int(hist[i, j])})
    norm = LogNorm(vmin=1, vmax=max(2, max(h.max() for _, _, h in histograms)))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#f2f2f2")
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True, sharey=True)
    for ax, (job, group, hist) in zip(axes.flat, histograms):
        mesh = ax.pcolormesh(xedges, yedges, np.ma.masked_equal(hist.T, 0), cmap=cmap, norm=norm, shading="flat", rasterized=True)
        ax.set_xscale("log")
        ax.set_ylim(0, 100)
        ax.set_title(f"{_label(group)}  |  {len(group):,} PLA objects")
        ax.set_xlabel("PLA object area (µm²; log scale)")
        ax.set_ylabel("Object area overlapping ChAT (%)")
    fig.subplots_adjust(left=0.08, right=0.87, top=0.90, bottom=0.15, wspace=0.18, hspace=0.32)
    cax = fig.add_axes([0.90, 0.23, 0.016, 0.59])
    fig.colorbar(mesh, cax=cax, label="PLA objects per area × overlap bin (log color)")
    fig.suptitle("Object size and overlap: where the measurements accumulate", x=0.01, ha="left")
    _footer(fig, "All PLA objects; identical bins and color limits. Grey = zero exported objects in a measurement bin. Associations may reflect segmentation geometry.")
    return fig, pd.DataFrame(rows)


def area_sensitivity(objects, thresholds=None):
    """Exploratory minimum-area thresholds; does not modify the input data."""
    pla = objects[objects.object_type == "PLA"]
    if thresholds is None:
        thresholds = np.r_[0.0, np.geomspace(pla.area.min(), max(1.0, pla.area.quantile(0.995)), 35)]
    rows = []
    for job, group in _groups(objects, "PLA"):
        for threshold in thresholds:
            selected = group[group.area >= threshold]
            rows.append({"job": str(job), "region": str(group.region.iloc[0]),
                         "min_area_um2": float(threshold), "retained_pla_count": len(selected),
                         "retained_pla_pct": 100 * len(selected) / len(group),
                         "pla_with_chat_count": int(selected.chat_present.sum()),
                         "pla_with_chat_pct": 100 * selected.chat_present.mean() if len(selected) else np.nan})
    table = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
    for job, group in table.groupby("job", sort=True):
        label = str(group.region.iloc[0]).replace("_", " · ").replace(" PLA", "")
        for ax, metric in zip(axes, ["retained_pla_pct", "pla_with_chat_pct"]):
            ax.plot(group.min_area_um2, group[metric], color=COLORS.get(str(job)),
                    linestyle=STYLES.get(str(job), "-"), linewidth=1.8, label=label)
    for ax in axes:
        ax.set_xscale("symlog", linthresh=0.003136)
        ax.set_xlabel("Hypothetical minimum PLA area (µm²; symlog scale)")
        ax.set_ylim(0, 100)
        ax.yaxis.set_major_formatter(PercentFormatter())
        ax.grid(alpha=0.18)
    axes[0].set_title("A  How much data would be retained?")
    axes[0].set_ylabel("All PLA objects retained")
    axes[1].set_title("B  How would ChAT association change?")
    axes[1].set_ylabel("Retained PLA objects with ChAT present")
    axes[0].legend(loc="lower left", fontsize=8)
    fig.suptitle("Sensitivity to excluding small PLA objects", x=0.01, ha="left")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.83, bottom=0.23, wspace=0.23)
    _footer(fig, "Exploratory filters only: no cutoff is validated or applied to the other figures. Baseline at area ≥ 0 retains all PLA objects.")
    return fig, table


def measurement_qc(objects):
    """Plot export features to review, never automatic exclusion flags."""
    rows = []
    metrics = [("Zero minimum diameter", lambda d: d.min_diameter == 0),
               ("Zero-width or zero-height box", lambda d: (d.xmin == d.xmax) | (d.ymin == d.ymax)),
               ("Zero inner area", lambda d: d.inner_area == 0)]
    for object_type in ["PLA", "ChAT"]:
        for job, group in _groups(objects, object_type):
            for name, fn in metrics:
                count = int(fn(group).sum())
                rows.append({"job": str(job), "region": str(group.region.iloc[0]), "object_type": object_type,
                             "feature": name, "count": count, "total_objects": len(group), "pct": 100 * count / len(group)})
    table = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), sharey=True)
    names = [m[0] for m in metrics]
    for ax, typ in zip(axes, ["PLA", "ChAT"]):
        subset = table[table.object_type == typ]
        arr = subset.pivot(index="job", columns="feature", values="pct").reindex(columns=names).to_numpy()
        im = ax.imshow(arr, cmap="cividis", vmin=0, vmax=100, aspect="auto")
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                ax.text(j, i, f"{arr[i,j]:.1f}%", ha="center", va="center", color="white" if arr[i,j] < 45 else "#111111")
        ax.set_xticks(range(3))
        ax.set_xticklabels(["Zero minimum\ndiameter", "Zero box width\nor height", "Zero inner\narea"])
        ax.set_yticks(range(4))
        ax.set_yticklabels([_label(g) for _, g in _groups(objects, typ)])
        ax.set_title(f"{typ} objects")
    fig.subplots_adjust(left=0.17, right=0.88, top=0.82, bottom=0.24, wspace=0.14)
    cax = fig.add_axes([0.91, 0.25, 0.018, 0.55])
    fig.colorbar(im, cax=cax, label="Objects with this exported feature (%)")
    fig.suptitle("Measurement features to inspect against the original image", x=0.01, ha="left")
    _footer(fig, "Features are retained, not classified as errors. Coordinate rounding and HALO settings may explain them; image review is needed.")
    return fig, table


def area_concentration(objects):
    """Show the exact contribution of the largest PLA shapes to summed area."""
    rows = []
    summaries = []
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), gridspec_kw={"width_ratios": [1.2, 1]})
    for job, group in _groups(objects, "PLA"):
        areas = np.sort(group.area.to_numpy())[::-1]
        ranks = np.arange(1, len(areas) + 1)
        cumulative = 100 * areas.cumsum() / areas.sum()
        fraction = 100 * ranks / len(areas)
        top_n = max(1, int(np.ceil(0.01 * len(areas))))
        share = float(cumulative[top_n - 1])
        axes[0].plot(fraction, cumulative, color=COLORS.get(str(job)),
                     linestyle=STYLES.get(str(job), "-"), linewidth=1.8, label=_label(group))
        rows.append(pd.DataFrame({"job": str(job), "rank_largest_first": ranks,
                                  "area_um2": areas, "objects_included_pct": fraction,
                                  "summed_area_included_pct": cumulative}))
        summaries.append({"job": str(job), "label": _label(group), "largest_object_count": top_n,
                          "largest_1pct_area_share_pct": share})
    table = pd.DataFrame(summaries)
    axes[0].set(xscale="log", xlim=(0.01, 100), ylim=(0, 100),
                xlabel="Largest PLA objects included (% of objects; log scale)",
                ylabel="Summed PLA object area included (%)",
                title="A  Area accumulates quickly in the upper tail")
    axes[0].axvline(1, color="#666666", linestyle=":", linewidth=1)
    axes[0].legend(loc="lower right", fontsize=8)
    axes[0].grid(alpha=0.18)
    for i, row in table.iterrows():
        axes[1].barh(i, row.largest_1pct_area_share_pct, color=COLORS.get(row.job))
        axes[1].text(row.largest_1pct_area_share_pct + 1, i,
                     f"{row.largest_1pct_area_share_pct:.1f}%", va="center", fontsize=10)
    axes[1].set_yticks(range(len(table)))
    axes[1].set_yticklabels(table.label)
    axes[1].invert_yaxis()
    axes[1].set(xlim=(0, 100), xlabel="Summed PLA object area (%)",
                title="B  Area contributed by the largest 1%")
    axes[1].grid(axis="x", alpha=0.18)
    fig.suptitle("A small number of large PLA shapes dominate area-based summaries", x=0.01, ha="left")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.83, bottom=0.23, wspace=0.53)
    _footer(fig, "Largest 1% uses ceil(0.01 × object count). Summed object area is not tissue coverage. Large shapes are review candidates, not established artifacts.")
    return fig, {"area_concentration": pd.concat(rows, ignore_index=True), "largest_one_percent": table}

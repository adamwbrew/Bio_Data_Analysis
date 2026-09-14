"""Compare box and tissue-outline analyses without pooling repeat detections.

Jobs are paired by the dataset owner's description and explicit region labels.
No Object Id join is performed: HALO assigns new IDs in each export. Signed
grid differences compare counts in aligned bins, not identified lost cells.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from .data import ROI_PAIRS, REGIONS, summarize
from .spatial import _edges, _hist

ROI_COLORS = {"box": "#737b85", "outline": "#0072B2"}
SIGNATURE_COLUMNS = ["object_type", "xmin", "xmax", "ymin", "ymax", "area",
                     "inner_area", "outer_area", "overlap_area", "overlap_pct",
                     "intensity", "min_diameter", "max_diameter", "median_diameter",
                     "pla_present", "chat_present"]


def _pairs(objects):
    required = {"roi_mode", "pair_id", "job"}
    if not required.issubset(objects.columns):
        raise ValueError("Load roi_mode='all' to obtain paired box/outline metadata.")
    for box_job, outline_job in ROI_PAIRS.items():
        box = objects.loc[objects.job == box_job]
        outline = objects.loc[objects.job == outline_job]
        if box.empty or outline.empty:
            raise ValueError(f"Missing paired exports {box_job}/{outline_job}; load roi_mode='all'.")
        if set(box.roi_mode) != {"box"} or set(outline.roi_mode) != {"outline"}:
            raise ValueError("Job and ROI-mode mapping is inconsistent.")
        yield REGIONS[box_job].replace(" PLA", ""), box, outline


def comparison_summary(objects):
    """One row per paired region and object type; changes are outline minus box."""
    stats = summarize(objects).set_index(["job", "object_type"])
    rows = []
    for pair_id, box, outline in _pairs(objects):
        for typ in ("PLA", "ChAT"):
            bj, oj = str(box.job.iloc[0]), str(outline.job.iloc[0])
            b, o = stats.loc[(bj, typ)], stats.loc[(oj, typ)]
            row = {"pair_id": pair_id, "box_job": bj, "outline_job": oj, "object_type": typ,
                   "box_count": int(b.n_objects), "outline_count": int(o.n_objects),
                   "count_delta": int(o.n_objects - b.n_objects),
                   "count_delta_pct": 100 * (o.n_objects - b.n_objects) / b.n_objects}
            for key, field in [("association", "other_marker_pct"), ("mean_overlap", "mean_object_overlap_pct"),
                               ("weighted_overlap", "area_weighted_overlap_pct"),
                               ("median_area", "median_area_um2"), ("median_intensity", "median_intensity")]:
                row["box_" + key] = float(b[field])
                row["outline_" + key] = float(o[field])
                row[key + "_delta"] = float(o[field] - b[field])
            rows.append(row)
    return pd.DataFrame(rows)


def roi_summary_figure(objects):
    table = comparison_summary(objects)
    pla = table.loc[table.object_type == "PLA"].reset_index(drop=True)
    chat = table.loc[table.object_type == "ChAT"].reset_index(drop=True)
    y = np.arange(len(pla))
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    for offset, frame, color, label in [(-0.16, pla, "#0072B2", "PLA objects"),
                                       (0.16, chat, "#c5c8cc", "ChAT objects")]:
        axes[0].barh(y + offset, frame.count_delta_pct, height=0.30, color=color, label=label)
        for i, row in frame.iterrows():
            axes[0].text(row.count_delta_pct - 0.10, i + offset, f"{row.count_delta_pct:+.2f}%",
                         ha="right", va="center", fontsize=8)
    axes[0].set_xlim(min(table.count_delta_pct.min() - 2.1, -2), max(0.8, table.count_delta_pct.max() + 2.1))
    axes[0].axvline(0, color="#888888", lw=1)
    axes[0].set_title("A  Change in exported counts")
    axes[0].set_xlabel("Count change relative to box selection (%)")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(pla.pair_id.str.replace("_", " · ", regex=False))
    axes[0].legend(fontsize=8, loc="lower left")
    for ax, metric, title in [(axes[1], "association", "B  PLA objects with ChAT present"),
                               (axes[2], "weighted_overlap", "C  Area-weighted PLA overlap")]:
        for i, row in pla.iterrows():
            b, o = row["box_" + metric], row["outline_" + metric]
            ax.plot([b, o], [i - 0.12, i + 0.12], color="#b0b0b0", lw=1.5)
            ax.scatter(b, i - 0.12, s=36, marker="o", color=ROI_COLORS["box"], label="Box" if i == 0 else None)
            ax.scatter(o, i + 0.12, s=42, marker="s", color=ROI_COLORS["outline"], label="Tissue outline" if i == 0 else None)
            ax.text(103, i, f"{o-b:+.2f} pp", va="center", fontsize=9)
        ax.set_xlim(0, 128)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_yticks(y)
        ax.set_yticklabels([])
        ax.set_title(title)
        ax.set_xlabel("Percentage (%)     Δ = outline − box")
        ax.legend(fontsize=8, loc="lower left")
    for ax in axes:
        ax.set_ylim(len(pla) - 0.1, -0.6)
        ax.grid(axis="x", alpha=0.15)
    fig.suptitle("How tissue outlining changes the exported measurements", x=0.01, ha="left", fontsize=17)
    fig.subplots_adjust(left=0.15, right=0.99, bottom=0.23, top=0.83, wspace=0.24)
    fig.text(0.02, 0.035, "Paired analyses of the same regions, not independent specimens. pp = percentage points.\n"
             "Counts are not tissue-normalized density; outline polygons and analyzed tissue area are not included.", fontsize=9)
    return fig, table


def paired_grid_data(objects, bin_width=500):
    """Conserved counts on identical bin edges for each pair, both object types."""
    if not np.isfinite(bin_width) or bin_width <= 0:
        raise ValueError("bin_width must be finite and positive.")
    frames = []
    for pair_id, box, outline in _pairs(objects):
        low_x = min(box.xmin.min(), outline.xmin.min())
        high_x = max(box.xmax.max(), outline.xmax.max())
        low_y = min(box.ymin.min(), outline.ymin.min())
        high_y = max(box.ymax.max(), outline.ymax.max())
        nx = np.floor(high_x / bin_width) - np.floor(low_x / bin_width) + 1
        ny = np.floor(high_y / bin_width) - np.floor(low_y / bin_width) + 1
        if nx * ny > 2_000_000:
            raise ValueError("Comparison grid too large; increase bin_width.")
        xe, ye = _edges(low_x, high_x, bin_width), _edges(low_y, high_y, bin_width)
        xl, yt = np.meshgrid(xe[:-1], ye[:-1])
        for typ in ("PLA", "ChAT"):
            b, o = box.loc[box.object_type == typ], outline.loc[outline.object_type == typ]
            bc, oc = _hist(b, xe, ye), _hist(o, xe, ye)
            frame = pd.DataFrame({"pair_id": pair_id, "object_type": typ,
                                  "box_job": box.job.iloc[0], "outline_job": outline.job.iloc[0],
                                  "x_left": xl.ravel(), "y_top": yt.ravel(), "bin_width": bin_width,
                                  "box_count": bc.ravel().astype(int), "outline_count": oc.ravel().astype(int),
                                  "count_delta": (oc - bc).ravel().astype(int)})
            frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def paired_spatial_figure(grid, pair_id):
    """Same coordinate frame, count scale and signed difference for a pair."""
    data = grid.loc[grid.pair_id == pair_id]
    if data.empty:
        raise ValueError(f"Unknown pair_id {pair_id!r}.")
    fig, axes = plt.subplots(2, 3, figsize=(16, 12))
    for i, typ in enumerate(["PLA", "ChAT"]):
        rows = data.loc[data.object_type == typ]
        xs, ys = np.sort(rows.x_left.unique()), np.sort(rows.y_top.unique())
        width = float(rows.bin_width.iloc[0])
        extent = [xs[0], xs[-1] + width, ys[-1] + width, ys[0]]
        def matrix(column):
            return rows.pivot(index="y_top", columns="x_left", values=column).reindex(index=ys, columns=xs).to_numpy()
        b, o, d = [matrix(c) for c in ("box_count", "outline_count", "count_delta")]
        norm = LogNorm(vmin=1, vmax=max(2, b.max(), o.max()))
        delta_max = max(1, abs(d).max())
        for j, arr in enumerate([b, o, d]):
            ax = axes[i, j]
            if j < 2:
                arr = np.ma.masked_equal(arr, 0)
                im = ax.imshow(arr, extent=extent, origin="upper", interpolation="nearest", cmap="cividis", norm=norm)
                label = f"{typ} count / bin (shared log scale)"
            else:
                arr = np.ma.masked_where((b + o) == 0, arr)
                im = ax.imshow(arr, extent=extent, origin="upper", interpolation="nearest", cmap="RdBu_r",
                               norm=Normalize(vmin=-delta_max, vmax=delta_max))
                label = f"Δ {typ} count / bin (outline − box)"
            ax.set_facecolor("#e3e6e9")
            ax.set_aspect("equal")
            ax.tick_params(labelsize=8)
            ax.set_xlabel("X (export coordinate units)", fontsize=9)
            ax.set_ylabel("Y (export coordinate units)", fontsize=9)
            title = [f"Box · job {rows.box_job.iloc[0]} · {int(b.sum()):,} {typ}",
                     f"Outline · job {rows.outline_job.iloc[0]} · {int(o.sum()):,} {typ}",
                     f"Net change: {int(d.sum()):+,} {typ} objects"][j]
            ax.set_title(title, fontsize=10)
            fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.055, pad=0.13).set_label(label, fontsize=8)
    fig.suptitle(pair_id.replace("_", " · ") + " — box and tissue-outline selections", x=0.02, ha="left", fontsize=17)
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.13, top=0.91, hspace=0.39, wspace=0.26)
    fig.text(0.02, 0.025, f"Object-coordinate reconstructions; {width:g} × {width:g} coordinate-unit bins, aligned to zero; no smoothing.\n"
             "Blue differences = fewer outline detections; red = more. Grey = no detections in the displayed run (both runs for the difference).\n"
             "A count difference is not an object match or a tissue boundary. The actual outline geometry is not exported.", fontsize=9)
    return fig


def exact_signature_comparison(objects):
    """Compare multiset counts of identical exported measurement signatures.

    Equality uses all 16 non-ID object fields, including bounds. Repeated
    signatures are counted with multiplicity; no Cartesian joins or rounded
    matches. Unmatched signatures can reflect changed measurements or detection,
    not necessarily removed/new physical objects. No per-object IDs are matched.
    """
    records = []
    for pair_id, box, outline in _pairs(objects):
        for typ in ("PLA", "ChAT"):
            b = box.loc[box.object_type == typ, SIGNATURE_COLUMNS].value_counts(sort=False).rename("box")
            o = outline.loc[outline.object_type == typ, SIGNATURE_COLUMNS].value_counts(sort=False).rename("outline")
            both = pd.concat([b, o], axis=1).fillna(0)
            shared = int(both.min(axis=1).sum())
            bn, on = int(b.sum()), int(o.sum())
            records.append({"pair_id": pair_id, "object_type": typ, "box_count": bn, "outline_count": on,
                            "identical_measurement_rows": shared,
                            "box_rows_without_identical_signature": bn - shared,
                            "outline_rows_without_identical_signature": on - shared,
                            "outline_identical_signature_pct": 100 * shared / on if on else np.nan})
    return pd.DataFrame(records)

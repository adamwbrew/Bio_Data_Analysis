"""Spatial reconstructions from HALO object tables, with explicit support maps.

These plots are coordinate summaries, not reconstructed microscopy pixels or
segmentation masks. Every object contributes its bounding-box midpoint to one
square bin. Coordinate units and the analyzed tissue/ROI footprint are unknown.
The functions consume the normalized DataFrame returned by ``load_data``.

Public functions return ``(figure, tidy_grid_dataframe)`` except
``bounding_box_zoom``, which returns a figure. No objects are sampled for grids.
The tidy table always preserves raw counts separately from display smoothing.
"""

from copy import copy

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
from scipy.ndimage import gaussian_filter


METRICS = {
    "pla_count": "PLA objects per square bin",
    "pla_fraction": "PLA objects with ChAT present (%)",
    "mean_overlap": "Mean PLA object-area overlap (%)",
}
_FACE = "#e9edf1"
_INK = "#172b46"
_BLUE = "#2485ab"
_ORANGE = "#d35424"
_MAX_GRID_CELLS = 2_000_000


def _check(df, bin_width, min_support, smoothing):
    required = {
        "job", "region", "object_type", "x", "y", "xmin", "xmax",
        "ymin", "ymax", "chat_present", "overlap_pct",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError("Missing normalized columns: " + ", ".join(sorted(missing)))
    if df.empty:
        raise ValueError("There are no objects to plot.")
    if not np.isfinite(bin_width) or bin_width <= 0:
        raise ValueError("bin_width must be finite and greater than zero.")
    if not np.isfinite(smoothing) or smoothing < 0:
        raise ValueError("smoothing must be a finite Gaussian sigma >= 0, in bins.")
    if int(min_support) != min_support or min_support < 1:
        raise ValueError("min_support must be a positive integer.")
    if not np.isfinite(df[["x", "y", "xmin", "xmax", "ymin", "ymax"]].to_numpy()).all():
        raise ValueError("Spatial coordinates must be finite; resolve invalid rows first.")
    if not df["object_type"].isin(["PLA", "ChAT"]).all():
        raise ValueError("object_type must contain only the normalized labels PLA and ChAT.")
    if not df["chat_present"].isin([True, False, 0, 1]).all():
        raise ValueError("chat_present must be boolean or numeric 0/1.")


def _edges(low, high, width):
    start = np.floor(low / width) * width
    stop = np.floor(high / width) * width + width
    bins = int(round((stop - start) / width))
    return start + np.arange(bins + 1) * width


def _hist(frame, xedges, yedges, weights=None):
    return np.histogram2d(
        frame["x"].to_numpy(), frame["y"].to_numpy(),
        bins=[xedges, yedges], weights=weights,
    )[0].T


def _smooth_mean(values, sigma):
    if sigma == 0:
        return values.copy()
    # Normalize the finite-grid kernel at edges: displayed counts are a local
    # kernel-weighted mean of bin counts, not conserved raw object totals.
    normalizer = gaussian_filter(np.ones_like(values), sigma=sigma, mode="constant")
    return gaussian_filter(values, sigma=sigma, mode="constant") / normalizer


def _divide(numerator, denominator):
    result = np.full_like(denominator, np.nan, dtype=float)
    np.divide(numerator, denominator, out=result, where=denominator > 0)
    return result


def _grids(df, bin_width, min_support, smoothing):
    """Return per-job grids; fixed edges are aligned to coordinate zero."""
    _check(df, bin_width, min_support, smoothing)
    grids = []
    for job, frame in df.groupby("job", sort=True, observed=True):
        labels = frame["region"].unique()
        if len(labels) != 1:
            raise ValueError("Each job must map to exactly one analysis region.")
        # Edges cover the bounding boxes; only midpoints are counted.
        nx = int(np.floor(frame["xmax"].max() / bin_width)
                 - np.floor(frame["xmin"].min() / bin_width) + 1)
        ny = int(np.floor(frame["ymax"].max() / bin_width)
                 - np.floor(frame["ymin"].min() / bin_width) + 1)
        if nx * ny > _MAX_GRID_CELLS:
            raise ValueError(
                "This bin width would create {:,} cells for job {}. "
                "Increase bin_width (try 250, 500, or 1000).".format(nx * ny, job)
            )
        xe = _edges(frame["xmin"].min(), frame["xmax"].max(), bin_width)
        ye = _edges(frame["ymin"].min(), frame["ymax"].max(), bin_width)
        pla = frame.loc[frame["object_type"] == "PLA"]
        chat = frame.loc[frame["object_type"] == "ChAT"]
        pc = _hist(pla, xe, ye)
        cc = _hist(chat, xe, ye)
        associated = _hist(pla.loc[pla["chat_present"].astype(bool)], xe, ye)
        overlap_sum = _hist(pla, xe, ye, weights=pla["overlap_pct"].to_numpy())
        smooth_pc = _smooth_mean(pc, smoothing)
        fraction = 100 * _divide(associated, pc)
        mean_overlap = _divide(overlap_sum, pc)
        fraction_display = 100 * _divide(_smooth_mean(associated, smoothing), smooth_pc)
        overlap_display = _divide(_smooth_mean(overlap_sum, smoothing), smooth_pc)
        # Support is always the UNSMOOTHED PLA count in the exact square bin.
        supported = pc >= min_support
        fraction_display[~supported] = np.nan
        overlap_display[~supported] = np.nan
        count_display = smooth_pc.copy()
        # A Gaussian kernel must not manufacture displayed detections in a bin
        # that contains no original PLA midpoint.
        count_display[pc == 0] = np.nan
        grids.append({
            "job": str(job), "region": str(labels[0]), "xe": xe, "ye": ye,
            "pla_count": pc, "chat_count": cc,
            "pla_with_chat_count": associated, "pla_fraction_pct": fraction,
            "mean_overlap_pct": mean_overlap, "supported": supported,
            "display_pla_count": count_display,
            "display_pla_fraction": fraction_display,
            "display_mean_overlap": overlap_display,
            "bin_width": float(bin_width), "min_support": int(min_support),
            "smoothing": float(smoothing),
        })
    return grids


def _tidy(grids, metric):
    frames = []
    for grid in grids:
        xe, ye = grid["xe"], grid["ye"]
        xleft, ytop = np.meshgrid(xe[:-1], ye[:-1])
        xright, ybottom = np.meshgrid(xe[1:], ye[1:])
        frame = pd.DataFrame({
            "job": grid["job"], "region": grid["region"],
            "x_left": xleft.ravel(), "x_right": xright.ravel(),
            "y_top": ytop.ravel(), "y_bottom": ybottom.ravel(),
            "bin_width_coordinate_units": grid["bin_width"],
            "pla_count": grid["pla_count"].ravel().astype(int),
            "chat_count": grid["chat_count"].ravel().astype(int),
            "pla_with_chat_count": grid["pla_with_chat_count"].ravel().astype(int),
            "pla_fraction_pct": grid["pla_fraction_pct"].ravel(),
            "mean_overlap_pct": grid["mean_overlap_pct"].ravel(),
            "raw_pla_support_pass": grid["supported"].ravel(),
            "min_pla_support": grid["min_support"],
            "gaussian_sigma_bins": grid["smoothing"],
            "display_metric": metric,
            "display_value": grid["display_" + metric].ravel(),
        })
        # Zero counts are zero exported detections in the rectangular grid;
        # analyzed ROI membership is unavailable for every bin.
        frame["has_any_exported_detection"] = (frame["pla_count"] + frame["chat_count"]) > 0
        frame["roi_membership_known"] = False
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def spatial_grid_data(df, bin_width=1000, min_support=10, smoothing=0,
                      metric="pla_fraction"):
    """Return all bins and raw numerators/denominators, without making a plot.

    ``smoothing`` is Gaussian sigma in bins. It affects ``display_value`` only;
    proportion displays use smoothed numerators / smoothed denominators, then
    mask bins whose *raw* PLA count is below ``min_support``. Mean overlap is
    the arithmetic mean of per-object percentages, not an area-weighted ratio.
    """
    _validate_metric(metric)
    return _tidy(_grids(df, bin_width, min_support, smoothing), metric)


def _validate_metric(metric):
    if metric not in METRICS:
        raise ValueError("metric must be one of: " + ", ".join(METRICS))


def _positive_norm(grids, key):
    positives = [g[key][np.isfinite(g[key]) & (g[key] > 0)] for g in grids]
    nonempty = [values for values in positives if len(values)]
    if not nonempty:
        return LogNorm(vmin=1, vmax=2)
    minimum = min(float(values.min()) for values in nonempty)
    maximum = max(float(values.max()) for values in nonempty)
    return LogNorm(vmin=min(1.0, minimum), vmax=max(maximum, min(1.0, minimum) * 1.01))


def _norms(grids, metric):
    chat_norm = _positive_norm(grids, "chat_count")
    metric_norm = (_positive_norm(grids, "display_pla_count") if metric == "pla_count"
                   else Normalize(vmin=0, vmax=100))
    return chat_norm, metric_norm


def _extent(grid):
    return [grid["xe"][0], grid["xe"][-1], grid["ye"][-1], grid["ye"][0]]


def _base_axis(ax):
    ax.set_facecolor(_FACE)
    ax.set_aspect("equal")
    ax.set_xlabel("X (exported coordinate units)", fontsize=9)
    ax.set_ylabel("Y (exported coordinate units)", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.ticklabel_format(axis="both", style="plain", useOffset=False)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=3))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=3))
    for spine in ax.spines.values():
        spine.set_color("#aab5c1")


def _layer(ax, grid, metric, chat_norm, metric_norm, background=True):
    extent = _extent(grid)
    if background:
        cc = np.ma.masked_where(grid["chat_count"] <= 0, grid["chat_count"])
        chat_cmap = copy(plt.get_cmap("Greys"))
        chat_cmap.set_bad(alpha=0)
        ax.imshow(cc, extent=extent, origin="upper", interpolation="nearest",
                  cmap=chat_cmap, norm=chat_norm, rasterized=True)
    values = np.ma.masked_invalid(grid["display_" + metric])
    cmap = copy(plt.get_cmap("inferno" if metric == "pla_count" else "viridis"))
    cmap.set_bad(alpha=0)
    image = ax.imshow(values, extent=extent, origin="upper", interpolation="nearest",
                      cmap=cmap, norm=metric_norm, rasterized=True)
    _base_axis(ax)
    return image


def _label(metric, smoothing):
    if metric == "pla_count":
        return ("PLA count / square bin (log color scale)" if smoothing == 0 else
                "Gaussian mean PLA count / bin (log color scale)")
    label = METRICS[metric]
    if smoothing:
        label += "\n(Gaussian-weighted numerator / denominator)"
    return label


def _footer(fig, bin_width, min_support, smoothing, proportion=False):
    text = ("Each position is a bounding-box midpoint. Square bins: {0:g} × {0:g} exported coordinate units. "
            "All rows contribute; no sampling.\n"
            "Pale background = no plotted detections; analyzed ROI boundaries and coordinate calibration are unavailable. "
            "Counts are not tissue density.").format(bin_width)
    if proportion:
        text += "\nProportion colors require ≥{} raw PLA objects/bin; blank color can indicate insufficient support.".format(min_support)
    if smoothing:
        text += ("\nGaussian display smoothing: σ={:g} bins ({:g} coordinate units), finite-grid edge normalization; "
                 "raw-count support is retained.").format(smoothing, smoothing * bin_width)
    else:
        text += "\nNo smoothing. No fluorescence intensities, tissue pixels, or object contours are reconstructed."
    fig.text(0.065, 0.017, text, fontsize=8, color="#425369", va="bottom", linespacing=1.45)


def coordinate_overview(df, bin_width=250):
    """Show the four jobs in their common exported coordinate system.

    ChAT count bins form a grayscale backdrop; PLA counts are the heatmap.
    Both logarithmic color scales are shared across all jobs. The field shown
    is the union of exported bounding extents, not the original image canvas.
    """
    grids = _grids(df, bin_width, 1, 0)
    chat_norm, metric_norm = _norms(grids, "pla_count")
    fig, ax = plt.subplots(figsize=(13, 13), facecolor="white")
    fig.subplots_adjust(left=0.10, right=0.91, bottom=0.23, top=0.89)
    for grid in grids:
        _layer(ax, grid, "pla_count", chat_norm, metric_norm)
        x0, x1, y1, y0 = _extent(grid)
        ax.text(x0 + 0.02 * (x1 - x0), y0 + 0.035 * (y1 - y0),
                "{} · job {}\n{:,.0f} PLA / {:,.0f} ChAT objects".format(
                    grid["region"], grid["job"], grid["pla_count"].sum(), grid["chat_count"].sum()),
                fontsize=9, va="top", color=_INK,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.94, "pad": 4})
    xmin = min(g["xe"][0] for g in grids)
    xmax = max(g["xe"][-1] for g in grids)
    ymin = min(g["ye"][0] for g in grids)
    ymax = max(g["ye"][-1] for g in grids)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymax, ymin)
    fig.suptitle("Where the exported objects occur", fontsize=23, weight="bold", color=_INK, y=0.97)
    fig.text(0.5, 0.93, "Coordinate reconstruction; not a microscopy image.  {:,} objects in their shared coordinate frame.".format(len(df)),
             ha="center", fontsize=11, color="#425369")
    chat_bar = fig.add_axes([0.11, 0.145, 0.32, 0.014])
    pla_bar = fig.add_axes([0.57, 0.145, 0.32, 0.014])
    fig.colorbar(plt.cm.ScalarMappable(norm=chat_norm, cmap="Greys"), cax=chat_bar,
                 orientation="horizontal", label="ChAT count / square bin (log color scale)")
    fig.colorbar(plt.cm.ScalarMappable(norm=metric_norm, cmap="inferno"), cax=pla_bar,
                 orientation="horizontal", label=_label("pla_count", 0))
    _footer(fig, bin_width, 1, 0)
    return fig, _tidy(grids, "pla_count")


def spatial_atlas(df, metric="pla_count", bin_width=250, min_support=10, smoothing=0):
    """Return region panels with identical spatial and color scales.

    Each panel keeps original coordinates and an equal aspect ratio. All
    panels use the same X/Y span, so shapes and distances are comparable.
    ``pla_fraction`` is the fraction of PLA objects with ChAT-present = 1;
    ``mean_overlap`` is the mean of PLA objects' area-overlap percentages.
    """
    _validate_metric(metric)
    grids = _grids(df, bin_width, min_support, smoothing)
    chat_norm, metric_norm = _norms(grids, metric)
    ncols = 2 if len(grids) > 1 else 1
    nrows = int(np.ceil(len(grids) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(13, 5 * nrows + 2.2),
                             squeeze=False, facecolor="white")
    fig.subplots_adjust(left=0.085, right=0.96, bottom=0.23, top=0.86,
                        hspace=0.36, wspace=0.26)
    xspan = max(g["xe"][-1] - g["xe"][0] for g in grids)
    yspan = max(g["ye"][-1] - g["ye"][0] for g in grids)
    for ax, grid in zip(axes.ravel(), grids):
        _layer(ax, grid, metric, chat_norm, metric_norm)
        xc = (grid["xe"][0] + grid["xe"][-1]) / 2
        yc = (grid["ye"][0] + grid["ye"][-1]) / 2
        ax.set_xlim(xc - xspan / 2, xc + xspan / 2)
        ax.set_ylim(yc + yspan / 2, yc - yspan / 2)
        pc = int(grid["pla_count"].sum())
        associated = int(grid["pla_with_chat_count"].sum())
        fraction = "{:.2f}%".format(100 * associated / pc) if pc else "undefined"
        subtitle = "{:,.0f} PLA objects · {} with ChAT present".format(pc, fraction)
        if metric != "pla_count":
            subtitle += "\n{:,} bins pass ≥{} PLA/bin".format(int(grid["supported"].sum()), min_support)
        ax.set_title("{} · job {}\n{}".format(grid["region"], grid["job"], subtitle),
                     fontsize=10, loc="left", color=_INK, pad=9)
    for ax in axes.ravel()[len(grids):]:
        ax.set_visible(False)
    heading = {
        "pla_count": "PLA count patterns across the four exported regions",
        "pla_fraction": "Where PLA objects have ChAT present",
        "mean_overlap": "How much of each PLA object overlaps ChAT",
    }[metric]
    if len(grids) != 4:
        heading = heading.replace("the four exported regions", "exported regions")
    fig.suptitle(heading, fontsize=21, weight="bold", color=_INK, y=0.97)
    fig.text(0.5, 0.928, "Coordinate reconstruction; not a microscopy image.  Shared spatial and color scales; gray = ChAT counts.",
             ha="center", fontsize=10, color="#425369")
    chat_bar = fig.add_axes([0.12, 0.15, 0.31, 0.014])
    metric_bar = fig.add_axes([0.58, 0.15, 0.31, 0.014])
    fig.colorbar(plt.cm.ScalarMappable(norm=chat_norm, cmap="Greys"), cax=chat_bar,
                 orientation="horizontal", label="ChAT count / square bin (log color scale)")
    fig.colorbar(plt.cm.ScalarMappable(norm=metric_norm, cmap=("inferno" if metric == "pla_count" else "viridis")),
                 cax=metric_bar, orientation="horizontal", label=_label(metric, smoothing))
    _footer(fig, bin_width, min_support, smoothing, metric != "pla_count")
    return fig, _tidy(grids, metric)


def region_heatmap(df, job, metric="pla_fraction", bin_width=1000,
                   min_support=10, smoothing=0):
    """Show a selected region's overlay and its raw PLA denominator/support.

    Color limits are calculated from all jobs provided in ``df``. The right
    panel's hatching marks nonempty bins below ``min_support``; the support
    colorbar gives raw PLA objects per bin, even when the left is smoothed.
    """
    _validate_metric(metric)
    grids = _grids(df, bin_width, min_support, smoothing)
    matches = [g for g in grids if g["job"] == str(job)]
    if not matches:
        raise ValueError("Unknown job {!r}; choose from {}".format(job, [g["job"] for g in grids]))
    grid = matches[0]
    chat_norm, metric_norm = _norms(grids, metric)
    support_norm = _positive_norm(grids, "pla_count")
    fig, axes = plt.subplots(1, 2, figsize=(14, 8.8), facecolor="white")
    fig.subplots_adjust(left=0.075, right=0.97, bottom=0.29, top=0.78, wspace=0.28)
    left_image = _layer(axes[0], grid, metric, chat_norm, metric_norm)
    axes[0].set_title("{}\nGrayscale backdrop: ChAT counts".format(METRICS[metric]), fontsize=11, loc="left", pad=10)
    counts = np.ma.masked_where(grid["pla_count"] == 0, grid["pla_count"])
    right_image = axes[1].imshow(counts, extent=_extent(grid), origin="upper",
                                 interpolation="nearest", cmap="cividis", norm=support_norm, rasterized=True)
    _base_axis(axes[1])
    axes[1].set_title("Raw PLA count: the local denominator\nHatching: 1–{} PLA objects/bin".format(min_support - 1)
                      if min_support > 1 else "Raw PLA count: the local denominator\nAll nonempty PLA bins meet the support threshold",
                      fontsize=11, loc="left", pad=10)
    low_y, low_x = np.where((grid["pla_count"] > 0) & ~grid["supported"])
    if min_support > 1:
        patches = [Rectangle((grid["xe"][ix], grid["ye"][iy]), bin_width, bin_width)
                   for iy, ix in zip(low_y, low_x)]
        if patches:
            collection = PatchCollection(patches, facecolor="none", edgecolor="#f4f5f7",
                                         linewidth=0, hatch="///", alpha=0.85, rasterized=True)
            axes[1].add_collection(collection)
    # Both panels keep exactly the same bounding extent and Y-down orientation.
    for ax in axes:
        ax.set_xlim(grid["xe"][0], grid["xe"][-1])
        ax.set_ylim(grid["ye"][-1], grid["ye"][0])
    fig.suptitle("{} · job {}".format(grid["region"], grid["job"]), fontsize=23,
                 weight="bold", color=_INK, y=0.965)
    total = int(grid["pla_count"].sum())
    retained = int(grid["pla_count"][grid["supported"]].sum())
    pass_bins = int(grid["supported"].sum())
    retention = "{:.1f}%".format(100 * retained / total) if total else "undefined"
    subtitle = ("Coordinate reconstruction; not a microscopy image.\n"
                "{:,} PLA objects; {:,} bins meet ≥{} raw PLA/bin, covering {} of PLA objects.").format(
                    total, pass_bins, min_support, retention)
    fig.text(0.5, 0.865, subtitle, ha="center", fontsize=10, color="#425369", linespacing=1.65)
    cb_left = fig.add_axes([0.105, 0.21, 0.335, 0.016])
    cb_right = fig.add_axes([0.595, 0.21, 0.335, 0.016])
    fig.colorbar(left_image, cax=cb_left, orientation="horizontal", label=_label(metric, smoothing))
    fig.colorbar(right_image, cax=cb_right, orientation="horizontal", label="Raw PLA count / square bin (log color scale)")
    _footer(fig, bin_width, min_support, smoothing, metric != "pla_count")
    return fig, _tidy([grid], metric)


def bounding_box_zoom(df, job, center=None, window=150):
    """Show actual exported rectangles plus their midpoints in a local window.

    With ``center=None``, select the most PLA-populated 1000-unit bin and use
    its center. This deterministic view deliberately emphasizes a hotspot;
    it is not a representative or randomly selected tissue field. All objects
    whose exported bounding boxes intersect the window are plotted, including
    zero-width/height boxes (visible as midpoint marks).
    """
    if not np.isfinite(window) or window <= 0:
        raise ValueError("window must be finite and greater than zero.")
    frame = df.loc[df["job"].astype(str) == str(job)]
    if frame.empty:
        raise ValueError("No objects found for job {!r}.".format(job))
    automatic = center is None
    if automatic:
        grid = _grids(frame, 1000, 1, 0)[0]
        iy, ix = np.unravel_index(np.argmax(grid["pla_count"]), grid["pla_count"].shape)
        # Focus on the nearest actual PLA midpoint to the hotspot bin center
        # so small windows do not accidentally show an empty patch.
        target = np.array([(grid["xe"][ix] + grid["xe"][ix + 1]) / 2,
                           (grid["ye"][iy] + grid["ye"][iy + 1]) / 2])
        pla = frame.loc[frame["object_type"] == "PLA"]
        candidates = pla if len(pla) else frame
        distance = ((candidates[["x", "y"]].to_numpy() - target) ** 2).sum(axis=1)
        center = candidates[["x", "y"]].iloc[int(np.argmin(distance))].to_numpy()
    if len(center) != 2 or not np.isfinite(center).all():
        raise ValueError("center must contain finite (x, y) coordinates.")
    x0, x1 = center[0] - window / 2, center[0] + window / 2
    y0, y1 = center[1] - window / 2, center[1] + window / 2
    visible = frame.loc[(frame["xmax"] >= x0) & (frame["xmin"] <= x1)
                        & (frame["ymax"] >= y0) & (frame["ymin"] <= y1)]
    fig, ax = plt.subplots(figsize=(10, 9), facecolor="white")
    fig.subplots_adjust(left=0.11, right=0.91, bottom=0.20, top=0.80)
    for label, color in [("ChAT", _BLUE), ("PLA", _ORANGE)]:
        selected = visible.loc[visible["object_type"] == label]
        patches = [Rectangle((r.xmin, r.ymin), r.xmax - r.xmin, r.ymax - r.ymin)
                   for r in selected.itertuples(index=False)]
        if patches:
            ax.add_collection(PatchCollection(patches, facecolor="none", edgecolor=color,
                                               linewidth=0.7, alpha=0.8, rasterized=True))
        ax.scatter(selected["x"], selected["y"], s=8 if label == "PLA" else 4,
                   c=color, marker="o" if label == "PLA" else "+", linewidths=0.7,
                   label="{}: {:,} exported boxes".format(label, len(selected)), rasterized=True)
    _base_axis(ax)
    ax.set_facecolor("#fafbfc")
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.legend(loc="upper right", fontsize=9, framealpha=0.95)
    fig.suptitle("A close look at exported object positions", fontsize=21, color=_INK, weight="bold", y=0.97)
    fig.text(0.5, 0.89, "{} · job {}\nRectangles = exported bounding boxes; marks = their midpoints".format(frame["region"].iloc[0], job),
             ha="center", fontsize=11, color="#425369", linespacing=1.6)
    selection = ("Automatic view: nearest PLA object to the center of the highest-count 1000-unit bin."
                 if automatic else "Window center supplied by the notebook user.")
    fig.text(0.075, 0.055,
             selection + "\n"
             "Window: {:g} × {:g} exported coordinate units, centered at ({:g}, {:g}). {:,} intersecting boxes shown.\n".format(
                 window, window, center[0], center[1], len(visible)) +
             "Coordinate reconstruction; not a microscopy image. Boxes are not cell outlines or segmentation masks.\n"
             "Some boxes have zero width or height; midpoint marks keep these exported objects visible. No rows were sampled.",
             fontsize=8.5, color="#425369", linespacing=1.55)
    return fig

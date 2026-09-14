"""Read and audit the four HALO exports without removing any object rows.

Coordinates remain in the uncalibrated units exported by HALO. ``x`` and ``y``
are bounding-box midpoints, not measured object centroids. Object areas and
diameters retain the micrometer units explicitly named by the CSV headers.
"""

from pathlib import Path
import csv
import hashlib

import numpy as np
import pandas as pd


FILE_PREFIX = "20260515_173745_2_jqw0p1_Ellie_PLA_05152026_2_backSUB.ome.ome.tiff_job"
FILE_SUFFIX = "PLA and ChaT Object Colocalization FL v3.0.0_object_results.csv"
REGIONS = {
    "4594": "22M_Pfkp-TDP43 PLA",
    "4595": "24M_Pfkp-TDP43 PLA",
    "4596": "22M_Hk1-TDP43 PLA",
    "4597": "24M_Hk1-TDP43 PLA",
}
KNOWN_ROWS = {"4594": 186343, "4595": 297584, "4596": 238084, "4597": 356745}
OBJECT_TYPES = {"PLA - TRITC_E1": "PLA", "ChAt-50 - FITC_E1": "ChAT"}
MEASUREMENTS = {
    "Area in PLA and Chat coloc": "overlap_area",
    "% Area in PLA and Chat coloc": "overlap_pct",
    "Area (µm²)": "area",
    "Inner Area (µm²)": "inner_area",
    "Outer Area (µm²)": "outer_area",
    "Average Intensity": "intensity",
    "Minimum Diameter (µm)": "min_diameter",
    "Maximum Diameter (µm)": "max_diameter",
    "Median Diameter (µm)": "median_diameter",
    "PLA - TRITC_E1 present": "pla_present",
    "ChAt-50 - FITC_E1 present": "chat_present",
    "XMin": "xmin",
    "XMax": "xmax",
    "YMin": "ymin",
    "YMax": "ymax",
}
SOURCE_COLUMNS = [
    "Image Location", "Analysis Region", "Algorithm Name", "Object Id",
    "Object Type", *MEASUREMENTS,
]


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_data(root=".", verify_hashes=True):
    """Return ``(objects, manifest)`` after reading all four source CSVs.

    ``root`` is the directory containing the CSVs. Source filenames, the image
    reference, and algorithm name live in the manifest, keyed by ``job``. The
    object key is therefore ``job`` + ``object_id`` (equivalently source filename
    + original Object Id). Each original row appears exactly once in ``objects``.

    Unexpected exports, schema changes, missing/nonfinite numeric values,
    duplicate/noninteger IDs, unknown labels, and nonbinary flags raise clear
    errors. Measurement anomalies (for example a zero reported diameter) are
    retained and described by :func:`quality_control`, never silently filtered.
    ``verify_hashes`` adds reproducibility hashes to the input manifest; it does
    not compare files with a remote copy or establish biological provenance.
    """
    root = Path(root).expanduser().resolve()
    expected_names = {FILE_PREFIX + job + FILE_SUFFIX for job in REGIONS}
    missing = sorted(name for name in expected_names if not (root / name).is_file())
    if missing:
        raise FileNotFoundError(
            "Place all four original HALO CSVs in DATA_ROOT. Missing:\n"
            + "\n".join(missing)
        )
    unexpected = sorted(
        path.name for path in root.glob("*object_results.csv")
        if path.name not in expected_names
    )
    if unexpected:
        raise ValueError(
            "Unexpected object-results exports; review the file manifest before "
            "including or excluding them:\n" + "\n".join(unexpected)
        )

    frames, manifest = [], []
    for job, expected_region in REGIONS.items():
        path = root / (FILE_PREFIX + job + FILE_SUFFIX)
        with path.open(encoding="utf-8-sig", newline="") as handle:
            header = next(csv.reader(handle), [])
        if header != SOURCE_COLUMNS:
            raise ValueError(f"Job {job}: unexpected CSV column names or order: {header}")

        # Fixed numeric dtypes fail on malformed strings. ID integrality is
        # checked explicitly before conversion; source IDs are far below 2**53.
        dtypes = {name: "float64" for name in ["Object Id", *MEASUREMENTS]}
        frame = pd.read_csv(path, encoding="utf-8-sig", dtype=dtypes)
        if len(frame) != KNOWN_ROWS[job]:
            raise ValueError(
                f"Job {job}: expected {KNOWN_ROWS[job]:,} original rows, "
                f"found {len(frame):,}. Review the input manifest if data changed."
            )
        if frame.isna().any().any():
            counts = frame.isna().sum()
            raise ValueError(f"Job {job}: missing values: {counts[counts > 0].to_dict()}")
        numeric = frame[list(dtypes)].to_numpy()
        if not np.isfinite(numeric).all():
            raise ValueError(f"Job {job}: nonfinite numerical measurements are present")

        ids = frame["Object Id"].to_numpy()
        if ((ids < 0) | (ids >= 2**53) | (ids != np.floor(ids))).any():
            raise ValueError(f"Job {job}: Object Id must be a nonnegative exact integer")
        if frame["Object Id"].duplicated().any():
            raise ValueError(f"Job {job}: duplicate Object Id values")
        for name in ["PLA - TRITC_E1 present", "ChAt-50 - FITC_E1 present"]:
            if not frame[name].isin([0, 1]).all():
                raise ValueError(f"Job {job}: {name!r} contains a nonbinary flag")
        if not frame["Object Type"].isin(OBJECT_TYPES).all():
            raise ValueError(f"Job {job}: unknown object types {frame['Object Type'].unique()}")
        unique_metadata = {}
        for column in ["Analysis Region", "Image Location", "Algorithm Name"]:
            values = frame[column].unique()
            if len(values) != 1 or not str(values[0]).strip():
                raise ValueError(f"Job {job}: expected one nonempty {column}")
            unique_metadata[column] = str(values[0])
        if unique_metadata["Analysis Region"] != expected_region:
            raise ValueError(f"Job {job}: region differs from the documented input manifest")

        manifest.append({
            "job": job,
            "region": expected_region,
            "source_file": path.name,
            "rows": len(frame),
            "bytes": path.stat().st_size,
            "sha256": _sha256(path) if verify_hashes else None,
            "image_location": unique_metadata["Image Location"],
            "algorithm": unique_metadata["Algorithm Name"],
            "ids_start_at_zero_and_contiguous": bool(
                np.array_equal(np.sort(ids), np.arange(len(frame)))
            ),
        })
        frame = frame.drop(columns=["Image Location", "Algorithm Name"])
        frame = frame.rename(columns={
            **MEASUREMENTS,
            "Analysis Region": "region",
            "Object Id": "object_id",
            "Object Type": "object_type",
        })
        frame.insert(0, "job", job)
        frame["object_id"] = frame["object_id"].astype("int64")
        frame["object_type"] = frame["object_type"].map(OBJECT_TYPES)
        frame["pla_present"] = frame["pla_present"].astype(bool)
        frame["chat_present"] = frame["chat_present"].astype(bool)
        frame["x"] = (frame["xmin"] + frame["xmax"]) / 2
        frame["y"] = (frame["ymin"] + frame["ymax"]) / 2
        frames.append(frame)

    objects = pd.concat(frames, ignore_index=True)
    return objects, pd.DataFrame.from_records(manifest)


def summarize(objects):
    """Descriptive summaries, one row per job/region/object type.

    All percentages name their denominator. ``other_marker_pct`` is the number
    of objects with the other marker present divided by objects of this type.
    ``area_weighted_overlap_pct`` is 100 * sum(overlap_area) / sum(area), computed
    WITHIN type. Neither is tissue area, a cell fraction, or biological replication.
    No inferential intervals or tests treat detected objects as independent animals.
    """
    records = []
    for (job, region, object_type), group in objects.groupby(
        ["job", "region", "object_type"], sort=True, observed=True
    ):
        n = len(group)
        other_marker = group["chat_present" if object_type == "PLA" else "pla_present"]
        area_sum = group["area"].sum()
        overlap_sum = group["overlap_area"].sum()
        zero_bbox = (group["xmin"] == group["xmax"]) | (group["ymin"] == group["ymax"])
        zero_diameter = group[["min_diameter", "median_diameter", "max_diameter"]].eq(0).any(axis=1)
        # Upper-tail contribution complements medians for these skewed shapes.
        n_largest = max(1, int(np.ceil(0.01 * n)))
        records.append({
            "job": job,
            "region": region,
            "object_type": object_type,
            "n_objects": n,
            "pla_present_count": int(group["pla_present"].sum()),
            "chat_present_count": int(group["chat_present"].sum()),
            "other_marker_count": int(other_marker.sum()),
            "other_marker_pct": 100.0 * other_marker.mean(),
            "chat_present_pct": 100.0 * group["chat_present"].mean(),
            "pla_present_pct": 100.0 * group["pla_present"].mean(),
            "overlap_positive_count": int(group["overlap_area"].gt(0).sum()),
            "area_sum_um2": area_sum,
            "overlap_area_sum_um2": overlap_sum,
            "area_weighted_overlap_pct": 100.0 * overlap_sum / area_sum if area_sum > 0 else np.nan,
            "mean_object_overlap_pct": group["overlap_pct"].mean(),
            "median_object_overlap_pct": group["overlap_pct"].median(),
            "median_area_um2": group["area"].median(),
            "q25_area_um2": group["area"].quantile(0.25),
            "q75_area_um2": group["area"].quantile(0.75),
            "p99_area_um2": group["area"].quantile(0.99),
            "max_area_um2": group["area"].max(),
            "largest_1pct_area_share_pct": 100.0 * group["area"].nlargest(n_largest).sum() / area_sum if area_sum > 0 else np.nan,
            "median_intensity": group["intensity"].median(),
            "p99_intensity": group["intensity"].quantile(0.99),
            "max_intensity": group["intensity"].max(),
            "median_diameter_um": group["median_diameter"].median(),
            "zero_min_diameter_count": int(group["min_diameter"].eq(0).sum()),
            "zero_median_diameter_count": int(group["median_diameter"].eq(0).sum()),
            "zero_max_diameter_count": int(group["max_diameter"].eq(0).sum()),
            "any_zero_diameter_count": int(zero_diameter.sum()),
            "any_zero_diameter_pct": 100.0 * zero_diameter.mean(),
            "zero_bbox_count": int(zero_bbox.sum()),
            "zero_bbox_pct": 100.0 * zero_bbox.mean(),
        })
    return pd.DataFrame.from_records(records)


def quality_control(objects):
    """Report measurement checks per export; retain every flagged object.

    Nonzero findings are review prompts, not automatic exclusion criteria.
    Percentage consistency uses an absolute tolerance of 0.001 percentage
    points to allow the rounded decimal fields in the source exports.
    """
    records = []
    for (job, region), group in objects.groupby(["job", "region"], sort=True, observed=True):
        area = group["area"].to_numpy()
        overlap_area = group["overlap_area"].to_numpy()
        overlap_pct = group["overlap_pct"].to_numpy()
        recomputed = np.divide(
            100.0 * overlap_area, area, out=np.full(len(group), np.nan), where=area > 0
        )
        positive = group["overlap_area"].gt(0)
        both_present = group["pla_present"] & group["chat_present"]
        diameter = group[["min_diameter", "median_diameter", "max_diameter"]]
        records.append({
            "job": job,
            "region": region,
            "n_objects": len(group),
            "duplicate_ids": int(group["object_id"].duplicated().sum()),
            "missing_values": int(group.isna().sum().sum()),
            "nonpositive_area": int(np.count_nonzero(area <= 0)),
            "negative_overlap_area": int(np.count_nonzero(overlap_area < 0)),
            "overlap_area_exceeds_area": int(np.count_nonzero(overlap_area > area + 1e-6)),
            "overlap_pct_outside_0_100": int(np.count_nonzero((overlap_pct < 0) | (overlap_pct > 100))),
            "overlap_pct_formula_mismatch": int(np.count_nonzero(~np.isclose(overlap_pct, recomputed, atol=0.001, rtol=0))),
            "max_overlap_pct_formula_error": float(np.nanmax(np.abs(overlap_pct - recomputed))),
            "overlap_presence_mismatch": int(positive.ne(both_present).sum()),
            "own_marker_absent": int(((group["object_type"].eq("PLA") & ~group["pla_present"]) | (group["object_type"].eq("ChAT") & ~group["chat_present"])).sum()),
            "any_zero_diameter": int(diameter.eq(0).any(axis=1).sum()),
            "any_negative_diameter": int(diameter.lt(0).any(axis=1).sum()),
            "diameter_order_mismatch": int(((group["min_diameter"] > group["median_diameter"]) | (group["median_diameter"] > group["max_diameter"])).sum()),
            "zero_bbox_width_or_height": int(((group["xmin"] == group["xmax"]) | (group["ymin"] == group["ymax"])).sum()),
            "reversed_bbox": int(((group["xmin"] > group["xmax"]) | (group["ymin"] > group["ymax"])).sum()),
            "nonzero_inner_area": int(group["inner_area"].ne(0).sum()),
            "area_partition_mismatch": int((group["inner_area"] + group["outer_area"] - group["area"]).abs().gt(0.001).sum()),
        })
    return pd.DataFrame.from_records(records)


def top_objects(objects, n=5):
    """Return the largest ``n`` objects of each type/export for image review.

    Selection ranks source measurements only; large shapes are not declared
    artifacts. Preserve job and Object Id so the source rows can be retrieved.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("n must be a positive integer")
    frames = [
        group.nlargest(n, "area").assign(area_rank=lambda frame: np.arange(1, len(frame) + 1))
        for _, group in objects.groupby(["job", "object_type"], sort=True, observed=True)
    ]
    return pd.concat(frames, ignore_index=True)

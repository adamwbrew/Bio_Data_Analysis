# PLA–ChAT scientific analysis and spatial explorer

Explore HALO fluorescence measurements from **six animals, two PLA arms and twelve tissue sections**. The repository includes executed Jupyter notebooks, a portable interactive Plotly app, high-resolution figures and reproducible analysis code. The current dataset contains **3,010,776 detections: 310,024 PLA and 2,700,752 ChAT**. Three animals are Running and three are Home cage; the HK1 and PFKP arms use different physical sections from each animal.

## Start here

| What you want to do | Open this |
|---|---|
| Explore the saved results without Python or raw data | Download the repository, unzip it, then open **`outputs/spatial_handoff/app/index.html`** in a web browser. [Download latest review ZIP](https://github.com/adamwbrew/Bio_Data_Analysis/archive/refs/heads/codex/standalone-plotly-app.zip) · [App file](outputs/spatial_handoff/app/index.html) |
| Read the scientific findings, methods and next actions | [Main scientific notebook](PLA_ChAT_scientific_analysis.ipynb) |
| Inspect every section, input check and method detail | [Full spatial handoff notebook](PLA_ChAT_spatial_handoff_analysis.ipynb) |
| Recompute with the source data or troubleshoot setup | [Viewing and rerunning guide](docs/quickstart.md) |

**Publication status:** the current app and ML updates are on the review branch in [PR #7](https://github.com/adamwbrew/Bio_Data_Analysis/pull/7). Use the review ZIP above while that pull request is open; `main` will contain these updates after merge.

**GitHub shows HTML source and disables notebook JavaScript.** Download the app to use its controls. It contains the plotting library, saved measurements and figures, and works offline. To use the notebook's `localhost` links, double-click [Launch_Spatial_Explorer.command](Launch_Spatial_Explorer.command) on Mac, or run `python scripts/serve_spatial_app.py`. Then open **[the local app](http://127.0.0.1:8050/)**. Serving the saved app needs only Python's standard library; installing analysis packages and obtaining raw CSVs are necessary only for recomputation.

## A useful path through the app

1. **Running comparison:** start with animal dots and group averages. Choose either subtraction direction, then inspect individual-animal sensitivity, size standardization and the normalization audit.
2. **Machine learning:** explore learned PLA size-and-shape groups, their native-coordinate maps, animal-level comparisons and source-linked microscopy review targets. Check the stability and exclusions before interpreting a pattern.
3. **Spatial explorer:** rotate and zoom the 3D measurement surface; overlay channels; view animal layers together or one at a time. Compare native coordinates with proposed feature alignment and inspect the marked cavity anchors.
4. **Scientific figures:** open or save full-resolution 300-DPI PNGs. The notebook also exports SVGs and the underlying tables.

The app is a **saved result snapshot**. Filters change the view; they do not refit models or ingest new files. A completed notebook **Run All** refreshes the app and its generation time. The local server binds only to `127.0.0.1` and serves the app rather than the repository.

## Findings to review with the scientist

The main outcome is the **percentage of PLA detections flagged ChAT-positive**, calculated within each animal and then averaged equally across the three animals in each condition.

| Home cage minus Running | HK1 arm | PFKP arm |
|---|---:|---:|
| Primary animal-level difference | **+22.01 percentage points** | **+22.58 percentage points** |
| After matching the PLA size distribution | +23.47 points | +24.07 points |
| After matching the learned size-and-shape mixture | +22.89 points | +23.30 points |

The positive direction persists in every cross-animal contrast and delete-one-per-group check. The morphology result suggests that the measured size-and-shape mixture alone does not explain the gap. It is an exploratory sensitivity analysis, not a biological discovery or a reliable running-status predictor. Six animals remain the biological sample size.

**Priority checks:** condition coincides with horizontal source-image layout, zero-diameter measurements affect morphology eligibility, and some learned partitions are less stable across animals. The notebook and app provide native-coordinate review targets so these questions can be checked against microscopy images and segmentation settings. Layout offsets applied only after detection cannot themselves alter the association fraction; the unresolved question is whether layout reflects acquisition/processing differences or simple export placement.

**Normalization is explicit.** The primary outcome is normalized by PLA count and gives each animal equal weight. Size and morphology standardization are separate sensitivity checks. Tissue area, ChAT availability, acquisition and segmentation remain uncorrected. Reversing subtraction changes the sign, not the evidence; relative percentages also depend on which group is the baseline.

## What the spatial images and ML mean

Maps reconstruct **detected-object bounding-box midpoints**, not original microscope pixels. Counts, association fractions and mean size answer different questions; use the measurement label and support threshold when reading a heat map. The 3D surface shows measurement height. The exploded layer stack uses an animal display order. **Neither is a measured anatomical volume.**

Feature alignment uses proposed cavity traces, an internal orientation cue where visible, and discrete flips checked against section geometry. It does not stretch, rotate continuously, warp or fill missing tissue. **0/12 sections currently have anatomical approval.** Individual proposed layers can be explored, but anatomical group means and subtraction stay gated until scientist review. The approval record is [handoff_section_reviews.json](metadata/handoff_section_reviews.json); computational proposals remain separate in [handoff_landmark_proposals.json](metadata/handoff_landmark_proposals.json).

Exploratory ML uses PLA area and diameter ratio, fits each arm separately, and balances animals during fitting. ChAT association and condition are attached after fitting. Leave-one-animal-out refits and two/four-group checks assess stability. Learned groups are measurement patterns, not validated cell types, and their native-coordinate maps require no anatomical registration.

## Inputs, outputs and reproducibility

The scientist's local `HALO_spatial_analysis_raw_data_handoff.pdf` supplies the animal assignments: **Home cage = 21F, 23M, 24M; Running = 22F, 21M, 22M**. The `entire Chat dataset/` folder contains thirteen exports. Interrupted job **4978** is superseded by **5403** and excluded from primary results. `Initial Data/` preserves the original box/outline exports and four repeated outline files; copies are audited, not concatenated.

**Raw CSVs, the local handoff PDF and original microscopy images are not distributed in Git.** The saved app and executed notebooks can be viewed without them. To recompute, obtain the raw exports from the dataset owner and follow the [quickstart](docs/quickstart.md). Original microscopy pixels, verified coordinate calibration, ROI masks and tissue-depth coordinates have not been supplied.

| Location | Contents |
|---|---|
| [outputs/spatial_handoff/app/](outputs/spatial_handoff/app/) | Portable offline app |
| [outputs/spatial_handoff/figures/](outputs/spatial_handoff/figures/) | 300-DPI PNGs and SVG figures |
| [outputs/spatial_handoff/tables/](outputs/spatial_handoff/tables/) | Compressed CSV tables, condition comparisons, ML audits and review coordinates |
| [outputs/spatial_handoff/spatial_arrays/](outputs/spatial_handoff/spatial_arrays/) | Lossless grid count arrays |
| [outputs/spatial_handoff/protocols/](outputs/spatial_handoff/protocols/) | Content-addressed input/code hashes and analysis parameters |
| [bio_visuals/](bio_visuals/) · [tests/](tests/) | Analysis, visualizations and verification code |

The [main notebook](PLA_ChAT_scientific_analysis.ipynb) contains all core results. The [full audit notebook](PLA_ChAT_spatial_handoff_analysis.ipynb) adds the complete static atlas and research methods. The [original visual notebook](PLA_ChAT_visual_analysis.ipynb), [original data dictionary and experiment diagram](docs/original_data_reference.md), and [tissue-outline review](docs/tissue_outline_review.md) preserve the earlier box/outline work. Statements about missing animal assignments in that historical material describe the evidence available before the handoff.

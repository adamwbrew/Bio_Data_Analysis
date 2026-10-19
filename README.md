# PLA–ChAT results and interactive app

Explore how PLA association with ChAT differs between **Running and Home-cage mice**, using measurements from six animals and twelve tissue sections.

**Start with the [illustrated results guide](docs/results_guide.md) to read the findings on GitHub, or [download the app and results](https://github.com/adamwbrew/Bio_Data_Analysis/archive/refs/heads/main.zip) to rotate, zoom and filter the data.** Neither requires running a notebook or obtaining the raw data.

## Open the interactive app

1. **[Download the app and results ZIP](https://github.com/adamwbrew/Bio_Data_Analysis/archive/refs/heads/main.zip).** Unzip `Bio_Data_Analysis-main.zip` on your computer.
2. Open the unzipped folder, then follow **`outputs` → `spatial_handoff` → `app` → `index.html`**. Right-click `index.html` → **Open With** → Chrome, Edge, Firefox or Safari. Open it in a browser, rather than a text editor.
3. Use the navigation across the top to choose a view:

| App section | What you can do |
|---|---|
| **Running comparison** | Overlay animal values and group averages; switch subtraction direction; check normalization and size sensitivity. |
| **Machine learning** | Explore size-and-shape groups, compare animals, inspect spatial maps and download source-image review coordinates. |
| **Explore sections** | Rotate and zoom 3D measurement views, toggle PLA/ChAT channels, and overlay all animal layers or step through them individually. |
| **Figure library** | Filter figures by scientific question, open them at full resolution, and save 300-DPI PNGs. |

The downloaded app works **offline, without Python, Jupyter, a server or installation**. Allow a moment for the plots to initialize. It is a saved snapshot: filters explore the published results; they do not rerun the analysis.

**Only want the app file?** Right-click [this HTML file link](https://github.com/adamwbrew/Bio_Data_Analysis/raw/refs/heads/main/outputs/spatial_handoff/app/index.html) and choose **Save Link As…**. Keep the filename `index.html` (about 26 MB), then open it in your browser. A normal click may display source code because GitHub does not host the HTML as an app. The single file includes everything needed to explore the saved results.

**Already downloaded the repository?** Open `outputs/spatial_handoff/app/index.html` inside it. The Mac launcher and optional local-server instructions are in the [viewing guide](docs/quickstart.md#optional-use-the-notebooks-localhost-app-links). A `127.0.0.1:8050` link only works after starting that server on your own computer.

## Read the results directly on GitHub

**[Open the illustrated results guide →](docs/results_guide.md)** It explains the main findings with embedded figures, links to supporting checks, and concrete next steps for the scientist. No notebook viewer or download is needed.

For a specific question, open one of these images:

| Question | View the figure |
|---|---|
| How do Running and Home-cage averages differ? | [Animal values, means and subtraction](outputs/spatial_handoff/figures/13_condition_averages_difference.png) |
| Are the results normalized? | [Normalization and both subtraction directions](outputs/spatial_handoff/figures/17_normalization_both_directions.png) |
| What did the machine learning find? | [PLA size-and-shape patterns](outputs/spatial_handoff/figures/21_morphology_patterns.png) |
| How are the tissue layers aligned? | [Before/after alignment proposals](outputs/spatial_handoff/figures/19_feature_alignment_before_after.png) |
| Where should the scientist investigate next? | [Source-image layout check](outputs/spatial_handoff/figures/16_condition_image_layout.png) · [Review priorities](docs/results_guide.md#what-to-check-next) |

[Browse all high-resolution figures](outputs/spatial_handoff/figures/) · [Detailed viewing and troubleshooting instructions](docs/quickstart.md)

## Download a notebook for the complete analysis

**`PLA_ChAT_spatial_handoff_analysis.ipynb` is too large for GitHub's notebook preview.** This is a viewing limitation; its executed results are saved. The main notebook may also fail to preview because of its embedded figures and interactive data. Use the app or illustrated guide above for easy access.

| Notebook | Purpose | Download |
|---|---|---|
| `PLA_ChAT_scientific_analysis.ipynb` | Main report: findings, methods, figures, interactive controls and next actions. About **40 MB**. | [Download main notebook](https://github.com/adamwbrew/Bio_Data_Analysis/raw/refs/heads/main/PLA_ChAT_scientific_analysis.ipynb?download=1) |
| `PLA_ChAT_spatial_handoff_analysis.ipynb` | Full audit: every section, input checks, additional figures and detailed methods. About **76 MB**. | [Download full audit notebook](https://github.com/adamwbrew/Bio_Data_Analysis/raw/refs/heads/main/PLA_ChAT_spatial_handoff_analysis.ipynb?download=1) |

Open a downloaded notebook in **VS Code with the Jupyter extension** or **JupyterLab** to read its saved outputs. **You do not need Run All or the raw CSVs just to read the results.** Use the browser app if embedded interactive controls do not render. For the complete folder with linked figures and code, [download the repository ZIP](https://github.com/adamwbrew/Bio_Data_Analysis/archive/refs/heads/main.zip) and unzip it first.

If a notebook download link displays JSON text, right-click the link → **Save Link As…** and preserve its `.ipynb` filename, or use the ZIP instead.

To analyze the source data again, follow the separate [recomputation instructions](docs/quickstart.md#optional-recompute-from-the-source-measurements). Current results are published on `main`.

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

The downloadable notebooks above preserve the complete analysis. Use the [illustrated results guide](docs/results_guide.md) for a GitHub-readable overview. The [original visual notebook](PLA_ChAT_visual_analysis.ipynb), [original data dictionary and experiment diagram](docs/original_data_reference.md), and [tissue-outline review](docs/tissue_outline_review.md) preserve the earlier box/outline work. Statements about missing animal assignments in that historical material describe the evidence available before the handoff.

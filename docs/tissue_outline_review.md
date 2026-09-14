# Tissue-outline reanalysis of the PLA–ChAT data

The four new exports are usable with the existing analysis schema. They contain **1,067,281 object rows**, compared with **1,078,756** in the original box-selected exports. The dataset owner identifies the new runs as analyses of the same tissue sections using an outline instead of a rectangular selection. The notebook treats these as paired alternative analyses and defaults to the outlined data for its main figures.

## File correspondence and verification

| Canonical region | Box job | Outline job | Exact outline region name | Outline rows |
| --- | --- | --- | --- | ---: |
| `22M_Pfkp-TDP43` | 4594 | 4999 | `22M_TDP43-Pfkp_tight` | 184,684 |
| `24M_Pfkp-TDP43` | 4595 | 5000 | `24M_TDP43-Pfkp_tight` | 295,677 |
| `22M_Hk1-TDP43` | 4596 | 5001 | `22M_TDP43-Hk1_tight` | 237,144 |
| `24M_Hk1-TDP43` | 4597 | 5002 | `24M_TDP43-Hk1_tight` | 349,776 |

Every row in all eight files was read. The headers have the same 20 columns. Each file has one region label, the same algorithm name, finite numeric fields, no missing values, and unique IDs within that export. Areas, overlap fractions, marker flags, diameter ordering and bounding-box ordering pass the implemented consistency checks. Zero diameters and zero-span boxes remain in the tables as review features.

The source-image paths resolve to the same Windows path after normalizing repeated separators. The analysis name remains `PLA and ChaT Object Colocalization FL v3.0.0`. The actual threshold settings and algorithm configuration are not present, so an identical name does not independently verify every processing parameter. Exact filenames and SHA-256 hashes are recorded in the notebook's manifest and exported workbooks. The raw CSVs remain unchanged. [Local file inventory](../README.md#new-tissue-outline-exports)

## Changes in PLA measurements

All changes below are **outline minus box**. A percentage-point change is the difference between two percentages; it is distinct from a relative percentage change in a count.

| Region | PLA: box → outline | Relative count change | PLA with ChAT: box → outline | Change (percentage points) | Area-weighted overlap: box → outline |
| --- | ---: | ---: | ---: | ---: | ---: |
| `22M_Pfkp-TDP43` | 9,112 → 8,412 | −7.68% | 59.04% → 58.17% | −0.88 | 43.08% → 44.18% |
| `24M_Pfkp-TDP43` | 7,999 → 7,733 | −3.33% | 68.91% → 70.74% | +1.83 | 45.24% → 50.67% |
| `22M_Hk1-TDP43` | 15,282 → 14,984 | −1.95% | 42.44% → 42.73% | +0.29 | 38.02% → 37.41% |
| `24M_Hk1-TDP43` | 43,918 → 41,237 | −6.10% | 55.46% → 55.07% | −0.38 | 43.51% → 43.15% |

Across the four pairs, outlined selections contain **72,366 PLA objects**, 3,945 fewer than the original runs. ChAT objects decrease from 1,002,445 to **994,915**, a difference of 7,530. These are net differences between exports, not verified counts of objects removed by a mask.

The PLA-associated object fractions change by less than two percentage points in each region, while the `24M_Pfkp-TDP43` area-weighted overlap increases by **5.43 percentage points**. This is a useful distinction: counting objects gives each detection equal weight, whereas summed area gives large detections more influence. The notebook retains the size distributions, upper-tail area contribution, and minimum-area sensitivity plots to help review this behavior. None of these changes establishes a biological effect or validates a filtering threshold.

![Paired counts and overlap](../outputs/roi_comparison/figures/01_box_vs_outline_summary.png)

## Why the new tables are not assumed to be a simple subset

An exact comparison of each object's type, four coordinate bounds, nine area/intensity/diameter fields and two marker flags uses **16 non-ID fields**. IDs and image/region/algorithm labels are excluded. Repeated identical measurement signatures are counted with their multiplicities; they do not produce a many-to-many join.

| Region | Outline PLA rows identical to a box measurement signature | Outline PLA rows without an identical signature |
| --- | ---: | ---: |
| `22M_Pfkp-TDP43` | 8,225 / 8,412 (97.78%) | 187 |
| `24M_Pfkp-TDP43` | 7,499 / 7,733 (96.97%) | 234 |
| `22M_Hk1-TDP43` | 14,557 / 14,984 (97.15%) | 427 |
| `24M_Hk1-TDP43` | 40,793 / 41,237 (98.92%) | 444 |

Most PLA measurement rows have an exact counterpart, but some do not. An unmatched signature can reflect changed segmentation, changed measurements, a new detection, or an object excluded from the other run. This comparison cannot distinguish these explanations. In particular, it would be incorrect to label every unmatched row as a lost/new physical punctum. HALO object IDs restart within each job and are not correspondence keys across runs.

## Visual comparisons and their interpretation

The paired spatial figures show three panels for PLA and three for ChAT: box counts, outline counts, and the signed difference. Each pair uses identical square bin edges aligned to coordinate zero, equal X/Y aspect, and the same displayed coordinate extent. Box and outline count panels share a color scale within each object type; the difference panel uses a diverging scale centered on zero. The default bin width is 500 exported coordinate units. No smoothing is applied.

This design compares exactly the same coordinate bins while preserving the count totals. Red bins have more detections in the outlined analysis, blue bins have fewer, and grey difference bins contain no detections in either run. Equal counts do not imply identical detections. Empty coordinate bins do not identify tissue absence or the outline boundary.

The main notebook also generates an outlined-data coordinate reconstruction and supported PLA-fraction maps. Count maps use 250-unit bins for spatial detail; fraction maps default to 1,000-unit bins and require at least 10 raw PLA objects. The interactive explorer allows these display choices and the ROI mode to be varied. The minimum count is a visibility criterion, not a significance threshold.

The code rejects pooled box-plus-outline inputs for the main overview and spatial maps. This prevents repeat detections from being silently added or drawn over one another. Comparisons use the explicit paired functions instead.

## What a tissue outline does and does not provide here

The owner-provided selection method is new experimental context. However, these are still **object-result exports**: they have no polygon vertices, binary tissue mask, annotation area, or physical coordinate calibration. The notebook therefore does not fabricate an outline from a convex hull, claim objects per mm², or estimate the area of tissue removed. ChAT-object coverage is not a substitute for analyzed tissue area.

This distinction matters scientifically. Spatial association methods depend on the observation window and can require boundary corrections; the SODA study explicitly corrects neighbor counts near ROI boundaries. Its findings support requesting real ROI geometry before attempting tissue-aware spatial inference, rather than applying an arbitrary rectangle or point-cloud hull. [Lagache et al., 2018](https://www.nature.com/articles/s41467-018-03053-x)

HALO separately documents annotation-based spatial analyses, density heatmaps, summary data and object-pair outputs. These are useful potential additional exports, but their existence in the software does not mean they are included in the present CSVs. [Indica Labs, Spatial Analysis Module](https://indicalab.com/halo/halo-modules/spatial-analysis/)

For the next stage, retain the box and outline annotations, measured analyzed areas, and saved analysis settings alongside the original image and coordinate calibration. Review unchanged and changed measurement examples against the image, including boundaries and interior locations. These records would allow a defensible density denominator and help distinguish selection changes from other processing differences. Recording acquisition and analysis context is consistent with microscopy reproducibility guidance. [Montero Llopis et al., 2021](https://www.nature.com/articles/s41592-021-01156-w)

## Reproducing the review

Run [PLA_ChAT_visual_analysis.ipynb](../PLA_ChAT_visual_analysis.ipynb) from top to bottom. It loads the eight original CSVs, computes the comparisons directly, regenerates the selected-mode figures, and exports the values behind the charts. `ROI_MODE = 'outline'` is the default; use `'box'` to regenerate the earlier selection. The standalone loader defaults to box mode for compatibility with existing code.

The outlined analysis is saved in `outputs/outline/`; paired figures and their workbook are in `outputs/roi_comparison/`. The original artifacts at the root of `outputs/` remain the earlier box analysis. Both new workbooks include a source manifest and a Sources sheet. Run `python -m unittest discover -s tests -v` for denominator, pairing, grid-conservation, signature-multiplicity and loading checks.

## Sources

1. Dataset owner. Description of the four new exports as tissue-outline reruns of the same sections, September 14, 2026. This establishes the pairing context; the labels and source-image references provide supporting file evidence.
2. Local HALO object-result CSVs, jobs 4594–4597 and 4999–5002. Full filenames are linked in the [README inventory](../README.md). All numerical results above are calculated from the complete files, with input hashes preserved in the generated manifests.
3. Lagache, T. et al. [Mapping molecular assemblies with fluorescence microscopy and object-based spatial statistics](https://www.nature.com/articles/s41467-018-03053-x). *Nature Communications* 9, 698 (2018), particularly the description of the ROI boundary term.
4. Indica Labs. [Spatial Analysis Module](https://indicalab.com/halo/halo-modules/spatial-analysis/). Official documentation, accessed September 14, 2026.
5. Montero Llopis, P. et al. [Best practices and tools for reporting reproducible fluorescence microscopy methods](https://www.nature.com/articles/s41592-021-01156-w). *Nature Methods* 18, 1463–1476 (2021).

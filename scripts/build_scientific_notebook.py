"""Build the concise, interactive-first main notebook plus the full audit notebook."""
from pathlib import Path
from textwrap import dedent
import runpy
import copy
import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
audit = runpy.run_path(str(ROOT / 'scripts/build_spatial_handoff_notebook.py'))
source = audit['cells']
cells = []


def md(text):
    cells.append(nbf.v4.new_markdown_cell(dedent(text).strip()))


def code(text):
    cells.append(nbf.v4.new_code_cell(dedent(text).strip(), metadata={'scrolled': False}))


def reuse(prefix):
    cell = next(c for c in source if c.source.startswith(prefix))
    cells.append(copy.deepcopy(cell))


md('''
# PLA–ChAT · explore the sections and test the findings
**Main scientific notebook · interactive edition · 19 October 2026 · protocol 1.3.0**

### Open the Plotly app
**[Open the spatial app in your browser](http://127.0.0.1:8050/)** for a full-window view of the 3D surfaces, channel overlays, overlapping layers and scientific comparison figures.

**First launch:** double-click `Launch_Spatial_Explorer.command` in this repository on your Mac, or run `python scripts/serve_spatial_app.py` in a terminal. The browser opens automatically; leave the terminal open while using the link. No notebook renderer or additional Python packages are needed to explore the saved results.

**Without a server:** open the [portable app file](outputs/spatial_handoff/app/index.html) in a web browser. In VS Code, use your browser rather than opening the HTML as source text. The app is a saved snapshot; **Run All** refreshes it after the analysis finishes. The complete analysis and core results remain in this notebook.

Start with the **non-running versus Running comparison**: overlay animal values, subtract the group averages, and inspect size and individual-animal checks. Then use the spatial explorer to move the 3D view, inspect channels, and overlay or separate animal layers. This notebook then connects what you see to the animal-level findings and the next scientific decisions. The longer [spatial handoff notebook](PLA_ChAT_spatial_handoff_analysis.ipynb) holds the complete input audit, all section figures and detailed research methods.

**Open this file in JupyterLab or VS Code**, select the project Python environment, and **Run All**. Review and trust the notebook to enable embedded JavaScript. GitHub shows the saved static figures but suppresses interactive controls. The explorer is saved **inside this notebook**, with its data and plotting library; it does not require a separate website or a live kernel after rendering.

### Go directly to a result
[Non-running minus Running](#comparison) · [Interactive 3D / channels / layers](#explore) · [Size sensitivity](#sizes) · [Layer alignment review](#registration) · [Image-review coordinates and decisions](#actions)

**What is reconstructed:** 12 physical sections from six animals, measured by HALO; 3,010,776 detected objects after excluding the interrupted duplicate. HK1–TDP43 and PFKP–TDP43 are different PLA arms and stay separate. Coordinates reconstruct detections, not microscope pixels. The 3D measurement view uses signal height; the exploded stack uses an ordinal layer axis. Neither is an anatomical volume.
''')
reuse('from pathlib import Path')
md('''
## Load the experiment once

The local handoff supplies the animal assignments: **Home cage = 21F, 23M, 24M; Running = 22F, 21M, 22M**. One section per animal per PLA arm is present in this ChAT panel. This resolves the earlier notebook's missing-group information. Raw CSVs in `entire Chat dataset/` and the original outline files in `Initial Data/` are supplied separately from Git.

The next two cells verify source files, remove the documented interrupted repeat, and build unsmoothed count grids. No review decisions, rotations or tissue-depth coordinates are invented. The full row-level QC tables are in the linked audit notebook and `outputs/spatial_handoff/tables/`.
''')
load = next(c.source for c in source if c.source.startswith('objects, manifest ='))
load = load.replace("table(repeat_audit, 'interrupted_repeat_audit')", "repeat_audit.to_csv(TABLES / 'interrupted_repeat_audit.csv.gz', index=False, compression='gzip')")
load = load.replace("display(metadata[['job','animal_id','running_group','pla_arm','rows','image_label','marker_panel']])", '')
code(load)
reuse('fine = h.make_grids')
# Reuse the complete comparison block from the audit to keep both notebooks identical.
start = next(i for i,c in enumerate(source) if c.source.startswith('<a id="comparison">'))
end = next(i for i,c in enumerate(source) if c.source.startswith('<a id="explore">'))
cells.extend(copy.deepcopy(source[start:end]))
reuse('<a id="explore"></a>')
reuse('display(HTML(explorer_html')
md('''
### Saved preview for static viewers

If the control center above shows a loading/instructions message, you are probably viewing an untrusted notebook or GitHub's static preview. Open the notebook locally and Run All. The figures below are **static companions**, not the controls. In the live explorer, the top tabs switch between the actual 3D, channel and layer views.
''')
code("show(hp.native_atlas(fine, metadata, 'TDP43–HK1'), '02_native_overlay_hk1')\nshow(hp.surface_preview(surfaces, metadata, 'TDP43–HK1'), '08_surface_companion_hk1')")
md('''
<a id="sizes"></a>
## Additional sensitivity · size filters and marker availability

Larger PLA objects are more likely to be ChAT-positive. The left panels compare association within fixed size categories; the right panels show how the group summary changes when small PLA objects are removed. Thin lines retain individual animals. These are exploratory cutoffs shared across animals, not validated classes of aggregates.

Use **Mean PLA area** and **Large PLA ≥1 µm² (%)** in the 3D explorer to locate size-rich areas. Then switch to the channel overlay and ChAT-positive/negative PLA subsets. Inspect the corresponding source-image coordinates before changing segmentation settings or interpreting these areas biologically.
''')
reuse('size_bins, sensitivity =')
md('''
### Proximity is a second description, not proof of enrichment

The next curves ask how far each PLA midpoint is from the nearest ChAT midpoint. They complement the exported overlap flag. Dashed curves restrict the analysis to a conservative interior subset to reveal edge sensitivity; axes remain in uncalibrated export units. ChAT object size, density and sampled geometry can explain short distances without a specific molecular interaction. A calibrated, geometry-aware null model is needed before calling enrichment.
''')
reuse('proximity_curves, proximity_summary =')
md('''
<a id="registration"></a>
## What the overlapping layers can tell us now

**Use the third explorer tab → Feature-aligned.** All 12 sections now have cavity-anchor proposals with internal orientation cues and geometry-ranked discrete flips. Inspect flat overlays, the 3D comparison stack and one-layer controls; show cavity outlines to see the common anchor. The before/after figure and full trace atlas above make each transform inspectable. **Native image coordinates** preserves the original offsets; **Candidate-centered** retains the old 9/12 automatic proposals for comparison.

This is a visual comparison of different animals and physical sections. Display order is 21F → 23M → 24M → 22F → 21M → 22M, not a serial anatomical sequence. Vertical spacing is a display control. The proposed geometry uses no scaling, continuous rotation or warping; cropped tissue cannot move the cavity anchor.

**0/12 sections are anatomically approved.** The separate proposal file records computational visual traces and flip scores without modifying investigator reviews. Review uncertain boundaries and low-margin flips in the original images before approving anatomical means or subtraction. Individual measurement layers can be inspected in the proposed frame, but are never pooled there.
''')
code("display(h.review_status(reviews)[['animal_id','pla_arm','candidate_found','ready','needs_review']])")
md('''
### Optional landmark review

Set `SHOW_REVIEW_EDITOR = True` in setup and rerun this cell when a scientist is ready to review the original images and native plots. Saving a reviewed section updates the JSON review file only on an explicit button click. Rerun the grid/review and explorer cells to incorporate saved decisions. Detailed candidate and orientation panels remain in the full audit notebook.
''')
reuse('if SHOW_REVIEW_EDITOR:')
reuse('reviews = h.load_reviews')
md('''
<a id="actions"></a>
## Actionable image windows and next decisions

The table gives native-coordinate windows selected by two declared rules: **largest PLA count** and **highest supported ChAT-positive fraction** within each section. They are review targets, not independently tested hotspots. Use animal, PLA arm, job and source filename together to request the matching original microscopy crop.
''')
code("queue = h.spatial_review_queue(metrics, metadata)\ntable(queue.drop(columns='source_file'), 'image_review_windows')\nqueue.to_csv(TABLES / 'image_review_windows_with_sources.csv.gz', index=False, compression='gzip')")
md('''
| Decision | Evidence to inspect | Action that would change the interpretation |
|---|---|---|
| Is the positive non-running minus Running gap biological or detection-dependent? | Animal dots, all nine deletion checks, size-standardized comparison and independent ChAT counts | Review matched fluorescence crops, negative controls and HALO thresholds; keep the two PLA arms separate. |
| Does image layout contribute to the condition difference? | Running and Home sections occupy different X sides in every source image | Review original channels, composite provenance, illumination/background and consistent segmentation thresholds in blinded equivalent fields. |
| Do sections share an anatomical reference? | Flat and exploded geometry layers, candidate polygon QC and orientation panels | Validate ventricular traces and discrete flips; supply calibration and ROI masks before group spatial differences. |
| Do large objects represent relevant structures or merged detections? | Mean size / large-PLA maps and review windows | Inspect segmentation against raw images; document any common size filter before recomputation. |
| Does spatial proximity exceed geometric opportunity? | Nearest-ChAT curves and independent ChAT abundance | Add physical scale, valid tissue support and a pre-specified geometry-aware null model. |
| Does the PLA pattern indicate a functional metabolic change? | Separate HK1/PFKP findings plus independent endpoints | Measure metabolism or enzyme activity; do not infer flux or binding from overlap alone. |

The analysis follows the supplied handoff's geometry constraints. Key research informing these decisions: [Söderberg et al., PLA](https://www.nature.com/articles/nmeth947); [Dunn et al., colocalization](https://pubmed.ncbi.nlm.nih.gov/21209361/); [Lazic, biological replication](https://pmc.ncbi.nlm.nih.gov/articles/PMC2817684/); [Barone et al., TDP-43/HK1 mechanism](https://pmc.ncbi.nlm.nih.gov/articles/PMC12992470/). These sources motivate controls and follow-up; they do not establish the mechanism in the current animals. The complete references and method rationale are in the [audit notebook](PLA_ChAT_spatial_handoff_analysis.ipynb).
''')
md('''
## Reproduce or continue the work

Run `python scripts/execute_spatial_handoff_notebook.py PLA_ChAT_scientific_analysis.ipynb` with the project environment. The notebook executes from raw files and includes all explorer data and JavaScript in its saved output. Figure exports, tables, count arrays and immutable protocol snapshots live under `outputs/spatial_handoff/`. The handoff PDF and raw CSVs remain local inputs supplied separately from Git.

The [full handoff notebook](PLA_ChAT_spatial_handoff_analysis.ipynb) supplies the detailed audit and figure atlas. The [original exploratory notebook](PLA_ChAT_visual_analysis.ipynb) retains the historical workflow. The earlier version of this scientific notebook is preserved in [Git history](https://github.com/adamwbrew/Bio_Data_Analysis/blob/b1e5a841446b938788ac373f40d39a14ca2df6d8/PLA_ChAT_scientific_analysis.ipynb).
''')
code("protocol_path = h.write_protocol(ROOT, OUT, manifest, reviews, PARAMETERS)\nprint('Protocol:', protocol_path.relative_to(ROOT))\nprint('Exploration is embedded above; current anatomical review approvals:', int(h.review_status(reviews).ready.sum()))")
reuse('app_path = write_standalone_app')
notebook = nbf.v4.new_notebook(cells=cells, metadata=copy.deepcopy(audit['nb'].metadata))
nbf.write(notebook, ROOT / 'PLA_ChAT_scientific_analysis.ipynb')
print(f'Built main scientific notebook: {len(cells)} cells')

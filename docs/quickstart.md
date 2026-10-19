# Open the results or rerun the analysis

## View the saved app — no analysis setup required

1. Download the repository using GitHub **Code → Download ZIP**, then unzip it. Use the branch or release containing the analysis you want; changes in an open pull request do not appear in `main` until merged.
2. In the downloaded folder, open **`outputs/spatial_handoff/app/index.html`** with Chrome, Edge, Firefox or Safari. In VS Code, use **Reveal in Finder/File Explorer**, then open the HTML with your browser; the editor otherwise shows its source.
3. Start with **Running comparison**, then **Machine learning** and the **Spatial explorer**. The app has all required plotting code, saved measurements and figure images embedded. It needs no raw CSVs, internet connection, notebook kernel or Python installation.

The HTML file is approximately 25 MB; allow a moment for its plots to initialize. GitHub's file preview is not a hosted application. You can share the downloaded `index.html` itself to share this saved result snapshot.

### Use the local links in the notebooks

The notebook links open `http://127.0.0.1:8050/`. Start the local app first:

| Platform | From the repository folder |
|---|---|
| Mac | Double-click `Launch_Spatial_Explorer.command`. If Python is unavailable, the launcher opens the portable file instead. If Finder does not permit launch, use the portable HTML directly. |
| Windows | Run `py -3 scripts/serve_spatial_app.py` in a terminal. |
| Linux / Mac terminal | Run `python3 scripts/serve_spatial_app.py`. |

These commands use **Python 3's standard library only**. They do not install or execute the scientific analysis. Leave the terminal open while using the local URL; Ctrl+C stops the server. If port 8050 is occupied, add `--port 8051` and use the URL printed in the terminal. To print the URL without opening a browser, add `--no-browser`.

The server only serves the app and its health check, and listens only on your own computer. The app generation time identifies the saved snapshot. Refreshing the page reloads the snapshot on disk; filters do not recompute scientific results.

## Read the notebooks

The [main scientific notebook](../PLA_ChAT_scientific_analysis.ipynb) is the primary report. The [full spatial handoff notebook](../PLA_ChAT_spatial_handoff_analysis.ipynb) adds all section figures and detailed methods. GitHub can show saved text, tables and static figures; it suppresses interactive JavaScript. If GitHub cannot preview a large notebook, download it and open it locally.

For the easiest full interactive view, use the portable app. For local notebook use, open the notebook in JupyterLab or VS Code with its Jupyter extension. Review the code before trusting the file to enable embedded HTML/JavaScript. A live Python kernel is needed for **Run All** and the optional scientist review editor; it is not needed to view the standalone app.

## Recompute from the source measurements

Recomputation is optional when reviewing the published results. The repository excludes the raw files; downloading GitHub alone is not enough for Run All.

### 1. Obtain and place the inputs

Ask the dataset owner for the original exports, the data documentation and related scientific information. Preserve filenames and this structure:

```text
Bio_Data_Analysis/
├── entire Chat dataset/              # 13 HALO CSV exports, including interrupted 4978
├── Initial Data/                     # Original exports; needed for the historical copy audit
├── HALO_spatial_analysis_raw_data_handoff.pdf  # Scientific handoff; local reference
├── metadata/                         # Tracked assignments, proposals and review records
├── requirements.txt
└── PLA_ChAT_scientific_analysis.ipynb
```

The current analysis validates all thirteen expanded exports and retains twelve completed sections. Do not manually delete the interrupted file: the pipeline audits its replacement. `Initial Data/` supports the historical copy check; the primary analysis reads `entire Chat dataset/`. The handoff PDF documents the assignments and protocol; it is not parsed during every run. Keep the tracked metadata and source filenames intact. Additional microscopy images, ROI masks and calibration should be supplied separately for anatomical review.

### 2. Create the analysis environment

Use a Python version compatible with [requirements.txt](../requirements.txt); the recorded execution uses Python 3.9.6. From the repository directory:

**Mac / Linux**

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m jupyter lab PLA_ChAT_scientific_analysis.ipynb
```

**Windows PowerShell**

```powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m jupyter lab PLA_ChAT_scientific_analysis.ipynb
```

In JupyterLab or VS Code, select this `.venv` environment as the notebook kernel. The setup cell prints its Python executable so you can confirm the selection. Loading approximately three million measurements requires substantial RAM and takes longer than viewing the saved app; runtime depends on your hardware.

### 3. Execute the main notebook

Choose **Run All**, or use the activated environment:

```bash
python scripts/execute_spatial_handoff_notebook.py PLA_ChAT_scientific_analysis.ipynb
```

On Windows without activating the environment, replace `python` with `.\.venv\Scripts\python.exe`. For the longer audit notebook:

```bash
python scripts/execute_spatial_handoff_notebook.py PLA_ChAT_spatial_handoff_analysis.ipynb
```

A completed run saves figures, tables, grid arrays, a protocol snapshot and the refreshed app under `outputs/spatial_handoff/`. The execution script saves the notebook only after all cells complete. The optional review editor remains off by default and Run All does not request manual approvals.

**Developers:** `scripts/build_scientific_notebook.py` regenerates both notebook sources and clears their saved outputs. Run it only when editing notebook structure, then execute both notebooks before publishing. For an ordinary data rerun, execute the existing notebook directly.

## Common fixes

| What you see | What to do |
|---|---|
| GitHub or VS Code shows HTML source | Open the downloaded app file with a web browser. |
| `127.0.0.1:8050` cannot be reached | Start the launcher/server, or open the portable HTML directly. |
| Port 8050 is already used | Run the server with `--port 8051`; follow the printed URL. |
| Notebook controls do not render | Use the standalone app, or review/trust the notebook in a local compatible viewer. |
| Run All reports missing CSVs | Place the thirteen expanded raw files in `entire Chat dataset/`; raw data are not in Git. |
| An import is missing | Install `requirements.txt` into the same environment selected as the notebook kernel. |
| Recent results are absent in GitHub `main` | Check whether the update is still on its pull-request branch; download that branch until merged. |
| Anatomical comparison remains unavailable | This is expected while 0/12 reviews are approved. Inspect proposed layers now; scientist-reviewed correspondence is required for group anatomical maps. |

See the [repository overview](../README.md) for the key findings and interpretation limits. Model groups are exploratory measurement patterns; 3D measurement heights and stacked animals are not measured anatomical depth.

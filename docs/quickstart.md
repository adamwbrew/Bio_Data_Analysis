# View the results in your browser

**Start with the [results guide](results_guide.md) to read the findings on GitHub, or download the interactive app below.** Neither route requires Run All, raw data, a notebook kernel or Python setup.

## Download and open the app

1. **[Download the app and results ZIP](https://github.com/adamwbrew/Bio_Data_Analysis/archive/refs/heads/main.zip)** and unzip `Bio_Data_Analysis-main.zip`.
2. In the unzipped folder, open **`outputs` → `spatial_handoff` → `app` → `index.html`**. Right-click the file and choose **Open With → Chrome, Edge, Firefox or Safari**. Allow a moment for its plots to initialize.
3. Start with **Running comparison**, then **Machine learning** and **Explore sections**. Rotate and zoom the 3D views, filter the animal layers, inspect normalization and download figures or tables.

The app includes its plotting library, saved measurements and figure images. It works offline after download, without Python, Jupyter or a server. You can share just `index.html` to share the same saved snapshot. The app and analysis updates are available on `main`.

**Smaller, single-file option:** right-click [the HTML file link](https://github.com/adamwbrew/Bio_Data_Analysis/raw/refs/heads/main/outputs/spatial_handoff/app/index.html), select **Save Link As…**, and keep the name `index.html` (approximately 26 MB). Then open it with your browser. Normal-clicking the raw HTML may show source code; GitHub does not host it as an app. If your browser adds `.txt`, rename the saved file to end in `.html`.

## Optional: read the saved notebooks locally

The **[full spatial handoff notebook](https://github.com/adamwbrew/Bio_Data_Analysis/raw/refs/heads/main/PLA_ChAT_spatial_handoff_analysis.ipynb?download=1) is approximately 76 MB and cannot be previewed on GitHub**. The **[main scientific notebook](https://github.com/adamwbrew/Bio_Data_Analysis/raw/refs/heads/main/PLA_ChAT_scientific_analysis.ipynb?download=1) is approximately 40 MB and may also fail to preview**. This is a viewer limitation, not missing results. Use the lightweight [results guide](results_guide.md) for a readable overview on GitHub and the downloaded app for interactive controls.

For the complete saved notebook report, download the repository ZIP and open the `.ipynb` file in **VS Code with its Jupyter extension**, or in an existing **JupyterLab** installation. Saved text, tables and static figures can be read **without running any cells, obtaining raw CSVs or starting an analysis kernel**. Review the code before trusting a notebook to enable embedded HTML/JavaScript. If a viewer does not render those controls, use the portable app.

Run All and the optional scientist review editor require a live Python environment. They are separate from viewing the saved results; setup for recomputation is below.

### Optional: use the notebooks' localhost app links

The notebook links to `http://127.0.0.1:8050/` require a local server. If you want to use those links instead of opening the portable HTML directly, download the full repository and start the app:

| Platform | From the repository folder |
|---|---|
| Mac | Double-click `Launch_Spatial_Explorer.command`. If Python is unavailable, the launcher opens the portable file instead. If Finder does not permit launch, use the portable HTML directly. |
| Windows | Run `py -3 scripts/serve_spatial_app.py` in a terminal. |
| Linux / Mac terminal | Run `python3 scripts/serve_spatial_app.py`. |

These commands use **Python 3's standard library only**. They do not install or execute the scientific analysis. Leave the terminal open while using the local URL; Ctrl+C stops the server. If port 8050 is occupied, add `--port 8051` and use the URL printed in the terminal. To print the URL without opening a browser, add `--no-browser`.

The server only serves the app and its health check, and listens only on your own computer. The app generation time identifies the saved snapshot. Refreshing the page reloads the snapshot on disk; filters do not recompute scientific results.

## Optional: recompute from the source measurements

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
| GitHub cannot render a large notebook | Read the [results guide](results_guide.md), download the app, or open the saved notebook locally. You do not need to rerun it. |
| GitHub or VS Code shows HTML source | Open the downloaded app file with a web browser. |
| `127.0.0.1:8050` cannot be reached | Start the launcher/server, or open the portable HTML directly. |
| Port 8050 is already used | Run the server with `--port 8051`; follow the printed URL. |
| Notebook controls do not render | Use the standalone app, or review/trust the notebook in a local compatible viewer. |
| Run All reports missing CSVs | Place the thirteen expanded raw files in `entire Chat dataset/`; raw data are not in Git. |
| An import is missing | Install `requirements.txt` into the same environment selected as the notebook kernel. |
| A downloaded app shows older results | Download a fresh copy from the link above and check its generation time. Existing local downloads do not update automatically. |
| Anatomical comparison remains unavailable | This is expected while 0/12 reviews are approved. Inspect proposed layers now; scientist-reviewed correspondence is required for group anatomical maps. |

See the [repository overview](../README.md) for the key findings and interpretation limits. Model groups are exploratory measurement patterns; 3D measurement heights and stacked animals are not measured anatomical depth.

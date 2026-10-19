"""Package the notebook's Plotly explorer as a portable, local browser app.

The exported page contains its data, Plotly library and figure pixels. Serving it
never exposes the repository or starts a kernel. Re-run a notebook to refresh it.
"""
import base64
from datetime import datetime, timezone
import hashlib
from html import escape
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
from urllib.parse import unquote, urlsplit


APP_ID = "pla-chat-spatial-explorer"
HEALTH_PATH = "/__halo_health__"
DEFAULT_PORT = 8050
GALLERY = (
    ("17_normalization_both_directions.png", "Normalization and both subtraction directions",
     "Primary fractions divide by PLA count within each animal. Tissue area and ChAT availability remain uncorrected."),
    ("18_within_animal_arm_comparison.png", "Is the condition gap different between PLA arms?",
     "The same animals contribute different physical sections. Compare arm gaps without assuming spatial correspondence."),
    ("19_feature_alignment_before_after.png", "Before and after cavity-and-shape alignment",
     "Proposed alignment keeps the internal anchor fixed and ranks discrete flips. Cropped tissue remains missing; these are not anatomical approvals."),
    ("20_proposed_cavity_traces_hk1.png", "HK1 arm: inspect each proposed cavity trace",
     "Full-resolution all-object reconstructions with visible proposal polygons; uncertain edges require original-image review."),
    ("20_proposed_cavity_traces_pfkp.png", "PFKP arm: inspect each proposed cavity trace",
     "These annotations supply common cavity anchors for the preview, without claiming validated anatomical registration."),
    ("13_condition_averages_difference.png", "Non-running minus Running: averages and subtraction",
     "Each dot is one animal. Positive differences mean Home cage is higher. "
     "Black intervals are delete-one-per-group sensitivity ranges, not confidence intervals."),
    ("14_condition_robustness.png", "Does one animal explain the difference?",
     "All nine cross-animal differences and all nine delete-one-per-group checks. "
     "These cells share animals; they are not independent replicates or matched pairs."),
    ("15_condition_size_standardization.png", "Does a different PLA size mixture explain the gap?",
     "The same size distribution is applied to every animal within each arm. "
     "This descriptive sensitivity check does not control all technical or biological differences."),
    ("16_condition_image_layout.png", "A source-image check to do before interpreting a running effect",
     "Condition and horizontal image layout coincide in these exports. Review the original "
     "channels, acquisition/composite construction and HALO thresholds. This pattern does not prove an artifact."),
    ("10_size_and_filter_sensitivity.png", "Punctum size and filter sensitivity",
     "Compare individual animals and group means across fixed PLA size bins and "
     "minimum-area filters. This checks whether the descriptive comparison depends "
     "on object size; the cutoffs are not validated biological classes."),
    ("11_nearest_chat_distance.png", "Nearest ChAT midpoint",
     "Each curve summarizes distances from PLA to the closest ChAT midpoint in "
     "the same section. Solid curves use all PLA; dashed curves apply a conservative "
     "edge check. Export units are uncalibrated; these are not molecular distances."),
)


def _gallery(root):
    """Embed original-resolution PNGs, with an explicit note for absent outputs."""
    cards = []
    for index, (name, title, caption) in enumerate(GALLERY):
        path = root / "outputs" / "spatial_handoff" / "figures" / name
        if path.is_file():
            encoded = base64.b64encode(path.read_bytes()).decode("ascii")
            content = (
                f'<img id="result-{index}" src="data:image/png;base64,{encoded}" '
                f'alt="{escape(title, quote=True)}" loading="lazy">'
                '<div class="result-actions">'
                f'<button data-full="result-{index}">View full resolution</button>'
                f'<button data-save="result-{index}" data-filename="{escape(name, quote=True)}">'
                'Save 300 dpi PNG</button></div>'
            )
        else:
            content = ('<p class="missing">This figure is not in this snapshot yet. '
                       'Run the scientific notebook to completion to generate it, '
                       'then refresh this app.</p>')
        cards.append('<article class="result-card"><h3>' + escape(title) + '</h3><p>'
                     + escape(caption) + '</p>' + content + '</article>')
    return "\n".join(cards)


def write_standalone_app(explorer_fragment, root, destination=None, comparison_fragment=None):
    """Write a complete offline page and return its absolute ``Path``.

    ``explorer_fragment`` is the trusted HTML returned by ``explorer_html``. It is
    intentionally executable, whereas all added metadata and captions are escaped.
    No source CSV or PDF is included beyond the explorer's measured grid payload.
    """
    if not isinstance(explorer_fragment, str) or not explorer_fragment.strip():
        raise ValueError("A nonempty explorer HTML fragment is required.")
    root = Path(root).resolve()
    destination = (root / "outputs" / "spatial_handoff" / "app" / "index.html"
                   if destination is None else Path(destination).resolve())
    destination.parent.mkdir(parents=True, exist_ok=True)
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    page = _PAGE.replace("__GENERATED__", escape(generated)).replace(
        "__GALLERY__", _gallery(root)).replace("__EXPLORER__", explorer_fragment).replace(
        "__COMPARISON__", comparison_fragment or '<p>Run the scientific notebook to generate the condition comparison.</p>')
    # Replace atomically so a running local server never serves a partial rebuild.
    temporary = destination.with_suffix(".html.tmp")
    temporary.write_text(page, encoding="utf-8")
    temporary.replace(destination)
    return destination


def document_identity(app_file):
    """Identify the exact snapshot served, without returning local filesystem paths."""
    return {"app": APP_ID, "document_sha256": hashlib.sha256(Path(app_file).read_bytes()).hexdigest()}


def make_handler(app_file):
    """Make a route-allowlisted HTTP handler; repository files are never routable."""
    app_file = Path(app_file).resolve()

    class AppHandler(BaseHTTPRequestHandler):
        server_version = "HALOSpatialExplorer/1"

        def log_message(self, format_string, *args):
            # Browser requests are routine; CLI startup and failures remain visible.
            pass

        def do_GET(self):
            self._respond(send_body=True)

        def do_HEAD(self):
            self._respond(send_body=False)

        def _respond(self, send_body):
            path = unquote(urlsplit(self.path).path)
            try:
                if path in ("/", "/index.html"):
                    content = app_file.read_bytes()
                    content_type = "text/html; charset=utf-8"
                elif path == HEALTH_PATH:
                    content = json.dumps(document_identity(app_file)).encode("utf-8")
                    content_type = "application/json; charset=utf-8"
                else:
                    self.send_error(404, "Only the spatial explorer is available on this server.")
                    return
            except (FileNotFoundError, PermissionError):
                self.send_error(503, "App snapshot unavailable; rerun the scientific notebook.")
                return
            self.send_response(200)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(content)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.end_headers()
            if send_body:
                self.wfile.write(content)

    return AppHandler


def create_server(app_file, port=DEFAULT_PORT):
    """Bind loopback only. Port zero is supported for isolated automated tests."""
    app_file = Path(app_file).resolve()
    if not app_file.is_file():
        raise FileNotFoundError(
            "The browser app has not been generated. Run PLA_ChAT_scientific_analysis.ipynb "
            "to completion, then launch the app again. Expected: " + str(app_file))
    return ThreadingHTTPServer(("127.0.0.1", port), make_handler(app_file))


_PAGE = '''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>PLA–ChAT | Spatial Explorer</title>
<style>
html{scroll-behavior:smooth}body{margin:0;background:#edf3f5;color:#173447;font:16px/1.55 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
.app-shell{max-width:1448px;margin:auto;padding:30px 24px 64px}.app-heading{display:flex;gap:24px;align-items:flex-end;justify-content:space-between;flex-wrap:wrap;margin-bottom:24px}
.app-heading h1{font-size:clamp(27px,4vw,42px);letter-spacing:-1px;line-height:1.13;margin:5px 0 12px}.app-heading p{max-width:740px;margin:0;color:#486777}.app-eyebrow{font-size:12px;letter-spacing:1.7px;font-weight:750;color:#087e8b}
.app-nav{display:flex;gap:8px;flex-wrap:wrap}.app-nav a,.result-actions button{border:1px solid #bdd0da;border-radius:8px;padding:9px 13px;color:#17485b;background:white;text-decoration:none;font:inherit;font-size:14px;cursor:pointer}.app-nav a:first-child{background:#087e8b;color:white;border-color:#087e8b}
.app-guide{background:white;border:1px solid #d0dfe5;border-radius:12px;padding:16px 20px;margin-bottom:18px}.app-guide summary{font-weight:750;cursor:pointer}.app-guide p{margin:10px 0}.app-guide ol{padding-left:22px;margin:10px 0}.app-guide li{margin:7px 0}
.app-shell h2{font-size:27px;letter-spacing:-.3px;margin:0 0 8px}.app-results{margin-top:38px}.app-results>p{color:#486777;margin-top:0}.result-card{padding:22px;background:white;border:1px solid #d0dfe5;border-radius:14px;margin-top:18px}.result-card h3{font-size:21px;margin:0 0 6px}.result-card p{max-width:1080px;color:#486777;margin:0 0 16px}.result-card img{display:block;width:100%;height:auto;background:white}.result-actions{display:flex;gap:10px;flex-wrap:wrap;margin-top:14px}.result-actions button:hover,.app-nav a:hover{filter:brightness(.94)}.missing{padding:18px;background:#fff3dc;border-radius:8px}
.app-note{color:#486777;font-size:13px;margin:18px 0 0}.app-footer{margin-top:30px;font-size:13px;color:#486777;border-top:1px solid #cbdce4;padding-top:18px}.halo-explorer{margin-left:auto!important;margin-right:auto!important}
@media(max-width:650px){.app-shell{padding:20px 12px}.result-card{padding:14px}.app-heading{gap:16px}}
</style></head><body><main class="app-shell">
<header class="app-heading"><div><div class="app-eyebrow">SCIENTIFIC DATA · INTERACTIVE PLOTLY APP</div>
<h1>PLA–ChAT Spatial Explorer</h1><p>Compare non-running and running animals, inspect the average difference, then explore the measured sections and overlapping layers. The notebook records the methods and reproducible results.</p></div>
<nav class="app-nav" aria-label="App navigation"><a href="#comparison">Running comparison</a><a href="#explore">Explore sections</a><a href="#results">Figures &amp; next steps</a><a href="#reading-guide">Reading guide</a></nav></header>
<details class="app-guide" id="reading-guide"><summary>Start here: three ways to inspect the data</summary><ol>
<li><strong>3D measurements:</strong> select a section and measurement, drag to rotate and scroll to zoom. Height represents the selected measurement, not tissue depth.</li>
<li><strong>Channel overlay:</strong> toggle ChAT and PLA detections, select positive or negative PLA and adjust opacity.</li>
<li><strong>Overlapping layers:</strong> choose a flat overlay or separated stack, then show all layers or step through one at a time. Section checkboxes control which animals appear.</li></ol>
<p><strong>Coordinate frames:</strong> native preserves exported positions. Candidate-centered is an unreviewed geometry check; reviewed alignment appears only after the required scientific review. Stack spacing and order are display choices, not measured physical depth.</p>
<p>These views reconstruct detections, not original microscopy pixels. Keep the two PLA arms separate. Blank bins are unknown or below the selected support threshold; exported x/y units have not been physically calibrated.</p></details>
<section id="comparison" aria-label="Running versus non-running comparison">__COMPARISON__</section>
<section class="app-guide" style="margin-top:18px"><h2>What to investigate next</h2><p><strong>Completed checks:</strong> equal-animal averages, all nine delete-one-per-group contrasts and common-size-distribution standardization. A positive association gap remains in both PLA arms.</p><ol>
<li><strong>Check image layout and detection first:</strong> Running sections occupy the lower-X side and Home-cage sections the higher-X side of the exports. Obtain original channels and composite/acquisition provenance; inspect illumination, background, ChAT segmentation and consistent thresholds using blinded, equivalent fields. If offsets were added only after detection, placement alone cannot change the association fractions.</li>
<li><strong>Locate the difference anatomically:</strong> review the ventricular landmarks, orientation, coordinate scale and sampling support for all 12 sections. The notebook can then show Home mean, Running mean, their difference and contributing-animal counts. Unsupported bins remain blank.</li>
<li><strong>Test a biological explanation:</strong> use negative PLA controls, matched section levels and additional animals. ChAT-positive PLA describes detected association; it does not establish binding, neuron loss, metabolic flux or a causal running effect.</li></ol><p>See the image-layout figure and native-coordinate review windows in the notebook for concrete source-image requests. Three animals per condition remain the biological sample size.</p></section>
<section id="explore" aria-label="Interactive spatial controls">__EXPLORER__</section>
<p class="app-note">This is a saved snapshot generated __GENERATED__. Filters work locally without a notebook kernel or internet connection. Changes to source data or review metadata do not recalculate automatically: run <strong>PLA_ChAT_scientific_analysis.ipynb</strong> to completion and refresh this page.</p>
<section class="app-results" id="results"><h2>Scientific results at a glance</h2><p>Use these companion figures to connect the spatial exploration to animal-level comparisons and sensitivity checks. Figures retain their original 300 dpi PNG pixels.</p>__GALLERY__</section>
<footer class="app-footer">Methods, data exclusions, review status and analysis code: <strong>PLA_ChAT_scientific_analysis.ipynb</strong> and the accompanying Git repository. The browser app is a visual companion to that analysis. Closing this tab leaves the local launcher running; stop it with Control+C in its Terminal window.</footer>
</main><script>
(function(){
  function imageBlob(image){const parts=image.src.split(',');const raw=atob(parts[1]);const bytes=new Uint8Array(raw.length);for(let i=0;i<raw.length;i++)bytes[i]=raw.charCodeAt(i);return new Blob([bytes],{type:'image/png'});}
  document.querySelectorAll('[data-save]').forEach(button=>button.addEventListener('click',()=>{
    const image=document.getElementById(button.dataset.save),url=URL.createObjectURL(imageBlob(image)),a=document.createElement('a');a.href=url;a.download=button.dataset.filename;a.click();setTimeout(()=>URL.revokeObjectURL(url),60000);
  }));
  document.querySelectorAll('[data-full]').forEach(button=>button.addEventListener('click',()=>{
    const image=document.getElementById(button.dataset.full),url=URL.createObjectURL(imageBlob(image));window.open(url,'_blank','noopener,noreferrer');setTimeout(()=>URL.revokeObjectURL(url),60000);
  }));
})();
</script></body></html>'''

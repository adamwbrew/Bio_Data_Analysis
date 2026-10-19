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
    ("21_morphology_patterns.png", "Machine learning: size-and-shape patterns",
     "Three descriptive partitions learned without condition or ChAT overlap. Animal-level summaries show what each pattern contributes."),
    ("22_morphology_stability.png", "Do the exploratory patterns survive a change in model?",
     "Leave-one-animal-out refits and two/four-cluster sensitivity check dependence on individual animals and the chosen partition. Agreement is stability, not prediction accuracy."),
    ("23_morphology_native_hk1.png", "HK1: where the morphology patterns occur",
     "Native-coordinate maps of cluster composition, normalized by eligible PLA count in each bin. No anatomical alignment or cross-animal pooling is assumed."),
    ("23_morphology_native_pfkp.png", "PFKP: where the morphology patterns occur",
     "Use these maps to request original-image crops and inspect segmentation. Colors are learned measurement patterns, not validated cell types."),
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
GALLERY_TOPICS = {
    "21_morphology_patterns.png": "morphology",
    "22_morphology_stability.png": "morphology",
    "23_morphology_native_hk1.png": "morphology",
    "23_morphology_native_pfkp.png": "morphology",
    "17_normalization_both_directions.png": "comparison",
    "18_within_animal_arm_comparison.png": "comparison",
    "19_feature_alignment_before_after.png": "alignment",
    "20_proposed_cavity_traces_hk1.png": "alignment",
    "20_proposed_cavity_traces_pfkp.png": "alignment",
    "13_condition_averages_difference.png": "comparison",
    "14_condition_robustness.png": "comparison",
    "15_condition_size_standardization.png": "comparison",
    "16_condition_image_layout.png": "review",
    "10_size_and_filter_sensitivity.png": "review",
    "11_nearest_chat_distance.png": "review",
}


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
        topic = escape(GALLERY_TOPICS.get(name, "review"), quote=True)
        cards.append('<article class="result-card" data-topic="' + topic + '"><h3>' + escape(title) + '</h3><p>'
                     + escape(caption) + '</p>' + content + '</article>')
    return "\n".join(cards)


def write_standalone_app(explorer_fragment, root, destination=None, comparison_fragment=None, ml_fragment=None):
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
        "__COMPARISON__", comparison_fragment or '<p>Run the scientific notebook to generate the condition comparison.</p>').replace(
        "__MACHINE_LEARNING__", ml_fragment or '<p>Run the scientific notebook to generate exploratory morphology results.</p>')
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
html{scroll-behavior:smooth}body{margin:0;background:#edf3f5;color:#173447;font:16px/1.55 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}[hidden]{display:none!important}.app-shell>[id]{scroll-margin-top:86px}
.app-shell{max-width:1448px;margin:auto;padding:30px 24px 64px}.app-heading{display:flex;gap:24px;align-items:flex-end;justify-content:space-between;flex-wrap:wrap;margin-bottom:24px}
.app-heading h1{font-size:clamp(27px,4vw,42px);letter-spacing:-1px;line-height:1.13;margin:5px 0 12px}.app-heading p{max-width:740px;margin:0;color:#486777}.app-eyebrow{font-size:12px;letter-spacing:1.7px;font-weight:750;color:#087e8b}
.app-nav{display:flex;gap:8px;flex-wrap:wrap;position:sticky;top:0;z-index:100;background:#edf3f5f5;padding:12px 0;margin:0 0 18px;border-bottom:1px solid #cbdce4}.app-nav a,.result-actions button{border:1px solid #bdd0da;border-radius:8px;padding:9px 13px;color:#17485b;background:white;text-decoration:none;font:inherit;font-size:14px;cursor:pointer}.app-nav a:hover,.app-nav a[aria-current=location]{background:#087e8b;color:white;border-color:#087e8b}.app-nav a:focus-visible,.result-actions button:focus-visible,.gallery-filter select:focus-visible{outline:3px solid #63bdc5;outline-offset:3px}.gallery-filter{display:flex;gap:14px;align-items:center;flex-wrap:wrap;background:white;padding:13px 16px;border:1px solid #d0dfe5;border-radius:9px}.gallery-filter label{font-weight:650}.gallery-filter select{font:inherit;padding:8px;border:1px solid #bdd0da;border-radius:6px;color:#173447;background:white}.gallery-filter span{font-size:13px;color:#486777}
.app-guide{background:white;border:1px solid #d0dfe5;border-radius:12px;padding:16px 20px;margin-bottom:18px}.app-guide summary{font-weight:750;cursor:pointer}.app-guide p{margin:10px 0}.app-guide ol{padding-left:22px;margin:10px 0}.app-guide li{margin:7px 0}
.app-route{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px;margin:16px 0}.app-route a{display:block;text-decoration:none;color:#173447;padding:14px;background:#eff6f7;border:1px solid #cbdfe3;border-radius:9px}.app-route a:hover{background:#dceff0}.app-route strong{display:block;color:#087e8b;margin-bottom:5px}.app-route span{font-size:13px;display:block}.app-guide details{border-top:1px solid #dce7ec;padding-top:12px;margin-top:14px}.app-terms{display:grid;grid-template-columns:160px 1fr;gap:9px 18px;font-size:14px}.app-terms dt{font-weight:700}.app-terms dd{margin:0;color:#486777}.app-next-step{padding:12px 15px;background:#eff7f6;border-left:4px solid #087e8b;border-radius:5px}.app-action-table{border-collapse:collapse;width:100%;font-size:14px}.app-action-table th,.app-action-table td{padding:10px 12px;text-align:left;vertical-align:top;border-bottom:1px solid #dce7ec}.app-action-table th{background:#edf3f5}.app-table-scroll{overflow-x:auto}@media(max-width:800px){.app-route{grid-template-columns:1fr 1fr}.app-terms{grid-template-columns:1fr;gap:3px}.app-terms dd{margin-bottom:9px}}@media(max-width:480px){.app-route{grid-template-columns:1fr}}
.app-shell h2{font-size:27px;letter-spacing:-.3px;margin:0 0 8px}.app-results{margin-top:38px}.app-results>p{color:#486777;margin-top:0}.result-card{padding:22px;background:white;border:1px solid #d0dfe5;border-radius:14px;margin-top:18px}.result-card h3{font-size:21px;margin:0 0 6px}.result-card p{max-width:1080px;color:#486777;margin:0 0 16px}.result-card img{display:block;width:100%;height:auto;background:white}.result-actions{display:flex;gap:10px;flex-wrap:wrap;margin-top:14px}.result-actions button:hover,.app-nav a:hover{filter:brightness(.94)}.missing{padding:18px;background:#fff3dc;border-radius:8px}
.app-note{color:#486777;font-size:13px;margin:18px 0 0}.app-footer{margin-top:30px;font-size:13px;color:#486777;border-top:1px solid #cbdce4;padding-top:18px}.halo-explorer{margin-left:auto!important;margin-right:auto!important}
@media(max-width:650px){.app-shell{padding:20px 12px}.result-card{padding:14px}.app-heading{gap:16px}.app-nav{gap:5px}.app-nav a{padding:7px 9px;font-size:12px}.app-shell>[id]{scroll-margin-top:130px}}@media(prefers-reduced-motion:reduce){html{scroll-behavior:auto}}
</style></head><body><main class="app-shell">
<header class="app-heading"><div><div class="app-eyebrow">SCIENTIFIC DATA · INTERACTIVE PLOTLY APP</div>
<h1>PLA–ChAT Spatial Explorer</h1><p>A guided view of the measured differences between home-cage and running mice. Start with the animal-level result, then use size, shape and spatial views to decide what to check in the original microscope images.</p></div></header>
<nav class="app-nav" aria-label="App navigation"><a href="#reading-guide">Start here</a><a href="#comparison">Running comparison</a><a href="#machine-learning">Machine learning</a><a href="#explore">Explore sections</a><a href="#results">Figure library</a></nav>
<section class="app-guide" id="reading-guide"><h2>One question, four ways to investigate it</h2>
<p><strong>Does the share of PLA detections associated with the ChAT marker differ between Home cage and Running?</strong> The main result compares animals. The other views help test possible measurement explanations and identify useful image-review locations.</p>
<div class="app-route">
<a href="#comparison"><strong>1. Compare animals</strong><span>Read the group averages and their difference. Check whether one animal or a different mix of object sizes changes the answer.</span></a>
<a href="#machine-learning"><strong>2. Examine size and shape</strong><span>See which measured object patterns contribute, then choose specific detections to review in the microscope images.</span></a>
<a href="#explore"><strong>3. Inspect the sections</strong><span>Locate signal, compare coverage and inspect proposed alignment. Use the 3D viewpoint buttons to look from each side.</span></a>
<a href="#results"><strong>4. Share a figure</strong><span>Choose a scientific question, save its high-resolution figure and pair it with the result wording from the comparison panel.</span></a>
</div>
<p><strong>Keep the unit of evidence clear:</strong> thousands of detected objects describe each section, but the biological comparison is between mice. The two PLA arms use different physical sections from the same animals.</p>
<details><summary>Study setup and key terms, in plain language</summary><dl class="app-terms">
<dt>PLA detection</dt><dd>A software-detected fluorescent spot from a proximity ligation assay, intended to report close proximity of the named targets. It is not a direct molecule count or proof of binding.</dd>
<dt>ChAT</dt><dd>Choline acetyltransferase, used here as a cholinergic marker. A segmented ChAT object is not automatically a neuron.</dd>
<dt>ChAT-positive PLA</dt><dd>A PLA detection flagged by HALO as associated with ChAT signal under the image-analysis settings. This is the main association measurement.</dd>
<dt>PLA arm</dt><dd>The target pair assayed: TDP43–HK1 or TDP43–PFKP. Read each arm separately; agreement across them is not replication in an independent group of animals.</dd>
<dt>Percentage points (pp)</dt><dd>The difference between two percentages. For example, 60% minus 40% is 20 percentage points. This is different from a relative percentage change.</dd>
<dt>Normalization</dt><dd>Making a stated comparison more comparable, such as dividing positive PLA by all PLA within each animal. No single adjustment removes every imaging or biological difference.</dd>
<dt>Spatial bin</dt><dd>A square map cell containing detections. Requiring a minimum PLA count hides estimates based on very few objects; it does not establish statistical confidence.</dd>
<dt>Morphology / cluster</dt><dd>Morphology means measured size and shape. A cluster groups similar measurements; it is a review aid, not a confirmed cell type.</dd>
</dl></details>
<details><summary>What the spatial images and 3D views represent</summary>
<p>These are reconstructions from detected-object coordinates, not original microscopy pixels. In <strong>3D measurements</strong>, height is the selected measurement. In the <strong>layer stack</strong>, height separates animals for viewing. These are display representations, not measured physical depth.</p>
<p><strong>Native</strong> keeps the exported positions. <strong>Feature-aligned</strong> centers proposed cavity outlines and applies proposed flips; it needs anatomical review. <strong>Candidate-centered</strong> is an earlier geometry-only check. <strong>Reviewed alignment</strong> admits only sections with the required scientific approvals.</p>
<p>Changing the camera does not change alignment. Blank map cells can mean no detected support or too few PLA objects for the selected metric; they do not demonstrate absent tissue. Exported x/y coordinates have not been physically calibrated.</p></details>
</section>
<section id="comparison" aria-label="Running versus non-running comparison">__COMPARISON__</section>
<section class="app-guide" style="margin-top:18px"><h2>Turn the results into a lab discussion</h2>
<p class="app-next-step"><strong>Most useful next step:</strong> review original fluorescence images and HALO detection settings in comparable regions from both groups. The app supplies specific objects and coordinates to make that review concrete.</p>
<div class="app-table-scroll"><table class="app-action-table"><thead><tr><th>Question for the scientist</th><th>Where to look</th><th>What this would resolve</th></tr></thead><tbody>
<tr><td>Were the two groups imaged and detected comparably?</td><td>Figure library → Source-image &amp; measurement checks; original channels and acquisition records.</td><td>Whether background, thresholds or image assembly could contribute to the observed group difference.</td></tr>
<tr><td>Are unusual detections real image features?</td><td>Machine learning → image-review target; download eligible-object and excluded-morphology QC targets.</td><td>Whether large, unusual-shape or invalid-diameter objects reflect real signal, merged detections or measurement limits.</td></tr>
<tr><td>Can the difference be located in matching anatomy?</td><td>Explore sections → Overlapping layers; review cavity landmarks, flips, scale and tissue coverage.</td><td>Whether comparable locations can be subtracted. The proposed visual alignment alone is not sufficient.</td></tr>
<tr><td>Does the pattern hold beyond these mice?</td><td>A follow-up experiment with additional animals, matched section levels and appropriate negative PLA controls.</td><td>Reproducibility and plausible biological explanations for the measured association.</td></tr>
</tbody></table></div>
<details><summary>Why the image-layout check comes first</summary><p>In these exports, Running sections occupy lower-X positions and Home-cage sections higher-X positions. Ask how the source images and export layout were assembled. The pattern flags something to check; it does not prove an imaging artifact. Offsets added only after detection cannot by themselves change the association fractions.</p></details>
<p><strong>When presenting:</strong> state the PLA arm, animals per group, measurement and denominator, group means and signed difference. Include what was checked and what remains unresolved. The comparison panel provides wording for the current selection; the notebook preserves the full methods.</p></section>
<section id="machine-learning" aria-label="Exploratory machine learning" style="margin:30px 0">__MACHINE_LEARNING__</section>
<section id="explore" aria-label="Interactive spatial controls">__EXPLORER__</section>
<p class="app-note">This is a saved snapshot generated __GENERATED__. Filters work locally without a notebook kernel or internet connection. Changes to source data or review metadata do not recalculate automatically: run <strong>PLA_ChAT_scientific_analysis.ipynb</strong> to completion and refresh this page.</p>
<section class="app-results" id="results"><h2>Figure library</h2><p>Choose the question you want to discuss, then open or save its original 300 dpi figure. These saved figures summarize the full analysis; filters in the interactive panels above do not change them. Include the figure caption and animal-level context when sharing.</p><div class="gallery-filter"><label for="figure-topic">Show figures about</label><select id="figure-topic"><option value="all">All results</option><option value="comparison">Running comparison &amp; normalization</option><option value="morphology">Machine learning &amp; morphology</option><option value="alignment">Layer alignment &amp; cavity review</option><option value="review">Source-image &amp; measurement checks</option></select><span id="figure-count" role="status" aria-live="polite"></span></div>__GALLERY__</section>
<footer class="app-footer">Methods, data exclusions, review status and analysis code: <strong>PLA_ChAT_scientific_analysis.ipynb</strong> and the accompanying Git repository. The browser app is a visual companion to that analysis. Closing this tab leaves the local launcher running; stop it with Control+C in its Terminal window.</footer>
</main><script>
(function(){
  const gallery=document.getElementById('results'),topic=document.getElementById('figure-topic'),cards=Array.from(gallery.querySelectorAll('.result-card'));
  function filterFigures(){let shown=0;cards.forEach(card=>{card.hidden=topic.value!=='all'&&card.dataset.topic!==topic.value;if(!card.hidden)shown++;});document.getElementById('figure-count').textContent=shown+' of '+cards.length+' figures';}
  topic.addEventListener('change',filterFigures);filterFigures();
  const navigation=Array.from(document.querySelectorAll('.app-nav a'));
  function markNavigation(){navigation.forEach(link=>{if(link.hash===location.hash)link.setAttribute('aria-current','location');else link.removeAttribute('aria-current');});}
  window.addEventListener('hashchange',markNavigation);markNavigation();
  function imageBlob(image){const parts=image.src.split(',');const raw=atob(parts[1]);const bytes=new Uint8Array(raw.length);for(let i=0;i<raw.length;i++)bytes[i]=raw.charCodeAt(i);return new Blob([bytes],{type:'image/png'});}
  document.querySelectorAll('[data-save]').forEach(button=>button.addEventListener('click',()=>{
    const image=document.getElementById(button.dataset.save),url=URL.createObjectURL(imageBlob(image)),a=document.createElement('a');a.href=url;a.download=button.dataset.filename;a.click();setTimeout(()=>URL.revokeObjectURL(url),60000);
  }));
  document.querySelectorAll('[data-full]').forEach(button=>button.addEventListener('click',()=>{
    const image=document.getElementById(button.dataset.full),url=URL.createObjectURL(imageBlob(image));window.open(url,'_blank','noopener,noreferrer');setTimeout(()=>URL.revokeObjectURL(url),60000);
  }));
})();
</script></body></html>'''

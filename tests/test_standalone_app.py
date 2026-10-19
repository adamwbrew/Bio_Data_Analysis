"""Portable packaging and strict routing for the standalone Plotly app."""
import base64
import errno
import importlib.util
import json
from pathlib import Path
import tempfile
import threading
import unittest
from urllib.error import HTTPError
from urllib.request import urlopen

from bio_visuals import standalone_app as app


class StandaloneAppTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        self.addCleanup(self.temporary.cleanup)

    def test_portable_page_embeds_original_pixels_and_trusted_explorer(self):
        figures = self.root / "outputs" / "spatial_handoff" / "figures"
        figures.mkdir(parents=True)
        png = b"\x89PNG\r\n\x1a\nfixture"
        for name, _, _ in app.GALLERY:
            (figures / name).write_bytes(png)
        fragment = '<div id="test-explorer"></div><script>window.testExplorer=true;</script>'
        comparison = '<div id="condition-comparison">Home cage minus Running</div>'
        ml = '<div id="morphology-test">Exploratory morphology</div>'
        path = app.write_standalone_app(fragment, self.root, comparison_fragment=comparison, ml_fragment=ml)
        html = path.read_text()
        self.assertEqual(path, self.root / "outputs/spatial_handoff/app/index.html")
        self.assertIn(fragment, html)
        self.assertEqual(html.count("data:image/png;base64," + base64.b64encode(png).decode()), len(app.GALLERY))
        self.assertIn(comparison, html)
        self.assertIn(ml, html)
        self.assertIn('href="#machine-learning"', html)
        self.assertNotIn("__MACHINE_LEARNING__", html)
        self.assertLess(html.index(comparison), html.index(fragment))
        self.assertNotIn(str(self.root), html)
        self.assertNotIn('src="http', html)
        self.assertNotIn("__GALLERY__", html)
        self.assertIn("do not recalculate automatically", html)
        self.assertIn("not measured physical depth", html)

    def test_missing_gallery_is_explicit_and_empty_explorer_is_rejected(self):
        path = app.write_standalone_app("<div>Ready</div>", self.root)
        self.assertEqual(path.read_text().count("not in this snapshot yet"), len(app.GALLERY))
        with self.assertRaises(ValueError):
            app.write_standalone_app("", self.root)

    def test_gallery_metadata_is_html_escaped(self):
        from unittest.mock import patch
        with patch.object(app, "GALLERY", (("missing.png", "A < B", "<script>untrusted</script>"),)):
            html = app.write_standalone_app("<div>Ready</div>", self.root).read_text()
        self.assertIn("A &lt; B", html)
        self.assertIn("&lt;script&gt;untrusted&lt;/script&gt;", html)
        self.assertNotIn("<script>untrusted</script>", html)

    def test_server_exposes_only_document_and_health_and_reloads_snapshot(self):
        path = app.write_standalone_app("<div>First snapshot</div>", self.root)
        server = app.create_server(path, 0)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        self.addCleanup(thread.join, 2)
        self.addCleanup(server.server_close)
        self.addCleanup(server.shutdown)
        url = "http://127.0.0.1:{}/".format(server.server_port)
        self.assertEqual(server.server_address[0], "127.0.0.1")
        for route in ("", "index.html", "index.html?refresh=1"):
            with urlopen(url + route) as response:
                self.assertIn(b"First snapshot", response.read())
                self.assertEqual(response.headers["Cache-Control"], "no-store")
        with urlopen(url.rstrip("/") + app.HEALTH_PATH) as response:
            self.assertEqual(json.load(response), app.document_identity(path))
        for route in ("README.md", "../README.md", "%2e%2e/README.md", ".git/config", "figures/"):
            with self.assertRaises(HTTPError) as caught:
                urlopen(url + route)
            self.assertEqual(caught.exception.code, 404)
        app.write_standalone_app("<div>Updated snapshot</div>", self.root)
        with urlopen(url) as response:
            self.assertIn(b"Updated snapshot", response.read())

    def test_missing_snapshot_and_port_collision_are_detected(self):
        with self.assertRaisesRegex(FileNotFoundError, "Run PLA_ChAT_scientific_analysis"):
            app.create_server(self.root / "missing.html", 0)
        path = app.write_standalone_app("<div>Ready</div>", self.root)
        server = app.create_server(path, 0)
        self.addCleanup(server.server_close)
        with self.assertRaises(OSError) as caught:
            app.create_server(path, server.server_port)
        self.assertEqual(caught.exception.errno, errno.EADDRINUSE)

    def test_launcher_checks_snapshot_identity_and_validates_port(self):
        import argparse
        spec = importlib.util.spec_from_file_location(
            "spatial_launcher", Path(__file__).resolve().parents[1] / "scripts/serve_spatial_app.py")
        launcher = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(launcher)
        path = app.write_standalone_app("<div>Ready</div>", self.root)
        server = app.create_server(path, 0)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        self.addCleanup(thread.join, 2)
        self.addCleanup(server.server_close)
        self.addCleanup(server.shutdown)
        url = "http://127.0.0.1:{}/".format(server.server_port)
        self.assertTrue(launcher.same_app_running(url, path))
        other = self.root / "other.html"
        other.write_text("Different snapshot")
        self.assertFalse(launcher.same_app_running(url, other))
        self.assertEqual(launcher.port_number("8050"), 8050)
        for bad in ("0", "65536", "-1"):
            with self.assertRaises(argparse.ArgumentTypeError):
                launcher.port_number(bad)


if __name__ == "__main__":
    unittest.main()

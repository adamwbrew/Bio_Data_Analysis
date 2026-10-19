"""Viewing saved results works without scientific dependencies or raw files."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class AppLauncherTests(unittest.TestCase):
    def test_fresh_checkout_launcher_imports_without_site_packages(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ('scripts/serve_spatial_app.py', 'bio_visuals/__init__.py',
                         'bio_visuals/standalone_app.py'):
                target = root / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(ROOT / name, target)
            result = subprocess.run(
                [sys.executable, '-S', str(root / 'scripts/serve_spatial_app.py'), '--help'],
                cwd=root, capture_output=True, text=True, timeout=15)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('--no-browser', result.stdout)

    @unittest.skipUnless(Path('/bin/zsh').is_file(), 'Mac launcher uses zsh')
    def test_unusable_python_opens_portable_app_in_path_with_spaces(self):
        with tempfile.TemporaryDirectory(prefix='PLA ChAT launch ') as directory:
            root = Path(directory)
            shutil.copyfile(ROOT / 'Launch_Spatial_Explorer.command', root / 'launch.command')
            python = root / '.venv/bin/python'
            python.parent.mkdir(parents=True)
            python.write_text('#!/bin/sh\nexit 1\n')
            python.chmod(0o755)
            app = root / 'outputs/spatial_handoff/app/index.html'
            app.parent.mkdir(parents=True)
            app.write_text('<html>Saved snapshot</html>')
            mock_bin = root / 'mock bin'
            mock_bin.mkdir()
            opener = mock_bin / 'open'
            opener.write_text('#!/bin/sh\nprintf "%s" "$1" > "$LAUNCH_TEST_OPEN_LOG"\n')
            opener.chmod(0o755)
            opened = root / 'opened.txt'
            env = dict(os.environ, PATH=str(mock_bin) + ':/usr/bin:/bin',
                       LAUNCH_TEST_OPEN_LOG=str(opened))
            result = subprocess.run(['/bin/zsh', str(root / 'launch.command')],
                                    env=env, capture_output=True, text=True, timeout=15)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(Path(opened.read_text()).resolve(), app.resolve())
            self.assertIn('portable spatial app', result.stdout)

    @unittest.skipUnless(Path('/bin/zsh').is_file(), 'Mac launcher uses zsh')
    def test_unusable_python_and_missing_snapshot_report_action(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            shutil.copyfile(ROOT / 'Launch_Spatial_Explorer.command', root / 'launch.command')
            python = root / '.venv/bin/python'
            python.parent.mkdir(parents=True)
            python.write_text('#!/bin/sh\nexit 1\n')
            python.chmod(0o755)
            result = subprocess.run(['/bin/zsh', str(root / 'launch.command')],
                                    capture_output=True, text=True, timeout=15)
            self.assertEqual(result.returncode, 1)
            self.assertIn('Download the complete repository', result.stdout)


if __name__ == '__main__':
    unittest.main()

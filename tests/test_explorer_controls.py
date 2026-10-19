"""Run the actual inline viewer's control/masking regression when Node is present."""
from pathlib import Path
import shutil
import subprocess
import unittest


class ExplorerControlTests(unittest.TestCase):
    def test_camera_controls_and_minimum_pla_filter(self):
        # The browser template uses modern JavaScript (including Array.at).
        # Prefer PATH, then the optional bundled runtime on desktop installs.
        candidates = [shutil.which('node'), str(Path.home() / '.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node')]
        node = None
        for candidate in candidates:
            if not candidate or not Path(candidate).is_file():
                continue
            try:
                version = subprocess.check_output([candidate, '--version'], text=True, timeout=5)
                if int(version.strip().lstrip('v').split('.')[0]) >= 18:
                    node = candidate
                    break
            except (OSError, ValueError, subprocess.SubprocessError):
                continue
        if not node:
            self.skipTest('Node.js 18+ is needed for the offline JavaScript control regression')
        script = Path(__file__).with_suffix('.js')
        result = subprocess.run([node, str(script)], capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('layered/channel modes passed', result.stdout)


if __name__ == '__main__':
    unittest.main()

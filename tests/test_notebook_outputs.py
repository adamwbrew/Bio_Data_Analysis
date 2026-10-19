"""Portable notebook outputs must survive a headless fresh-kernel run."""
import base64
import io
from pathlib import Path
import subprocess
import sys
import unittest

import nbformat as nbf
from PIL import Image

from bio_visuals.notebook_validation import validate_saved_outputs


class SavedOutputTests(unittest.TestCase):
    def notebook(self, source='show(fig, "example")', data=None):
        cell = nbf.v4.new_code_cell(source, execution_count=1)
        if data is not None:
            cell.outputs = [nbf.v4.new_output('display_data', data=data)]
        return nbf.v4.new_notebook(cells=[cell])

    def test_headless_figure_repr_and_dropped_images_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'unrendered Matplotlib'):
            validate_saved_outputs(self.notebook(data={'text/plain': '<Figure size 1000x500 with 2 Axes>'}))
        with self.assertRaisesRegex(ValueError, 'no embedded PNG'):
            validate_saved_outputs(self.notebook())
        with self.assertRaisesRegex(ValueError, 'invalid PNG'):
            validate_saved_outputs(self.notebook(data={'image/png': 'not-an-image'}))

    def test_portable_png_and_interactive_formats_are_accepted(self):
        output = io.BytesIO()
        Image.new('RGB', (4, 4), 'white').save(output, format='PNG')
        nb = self.notebook(data={'image/png': base64.b64encode(output.getvalue()).decode()})
        nb.cells.append(nbf.v4.new_code_cell('display(HTML(fragment))', execution_count=2,
                         outputs=[nbf.v4.new_output('display_data', data={'text/html': '<div>app</div>'})]))
        nb.cells.append(nbf.v4.new_code_cell('def show(fig):\n    display(Image(data=fig))', execution_count=3))
        self.assertEqual(validate_saved_outputs(nb), {'code_cells': 3, 'embedded_pngs': 1})

    def test_withheld_registration_map_is_optional_but_broken_figure_is_not(self):
        nb = self.notebook('if approved:\n    show(fig, "registered")\nelse:\n    display(status)',
                           data={'text/html': '<p>Spatial maps withheld pending review</p>'})
        nb.cells[0].metadata['static_figure_optional'] = True
        self.assertEqual(validate_saved_outputs(nb)['embedded_pngs'], 0)
        nb.cells[0].outputs[0].data = {'text/plain': '<Figure size 1000x500 with 5 Axes>'}
        with self.assertRaisesRegex(ValueError, 'unrendered Matplotlib'):
            validate_saved_outputs(nb)

    def test_real_notebook_helper_embeds_png_with_agg(self):
        root = Path(__file__).resolve().parents[1]
        script = r'''
import base64, tempfile
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import nbformat
from IPython.terminal.interactiveshell import TerminalInteractiveShell
from IPython.utils.capture import capture_output
shell = TerminalInteractiveShell.instance()
nb = nbformat.read('PLA_ChAT_scientific_analysis.ipynb', as_version=4)
setup = next(c.source for c in nb.cells if c.cell_type == 'code')
shell.run_cell(setup).raise_error()
with tempfile.TemporaryDirectory() as directory:
    shell.user_ns['OUT'] = Path(directory)
    with capture_output(display=True) as captured:
        shell.run_cell("fig, ax = plt.subplots(); ax.plot([0, 1], [1, 2]); show(fig, 'test')").raise_error()
    pngs = [o.data['image/png'] for o in captured.outputs if 'image/png' in o.data]
    assert len(pngs) == 1, [o.data.keys() for o in captured.outputs]
    assert base64.b64decode(pngs[0]) == (Path(directory)/'figures/test.png').read_bytes()
    assert (Path(directory)/'figures/test.svg').exists()
'''
        result = subprocess.run([sys.executable, '-c', script], cwd=root,
                                capture_output=True, text=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == '__main__':
    unittest.main()

"""Run All must refresh local code already cached by a notebook kernel."""

from pathlib import Path
import subprocess
import sys
import unittest


class NotebookSetupTests(unittest.TestCase):
    def test_all_notebook_cells_have_valid_python(self):
        import nbformat
        from IPython.core.inputtransformer2 import TransformerManager
        root = Path(__file__).resolve().parents[1]
        for name in ('PLA_ChAT_visual_analysis.ipynb', 'PLA_ChAT_scientific_analysis.ipynb'):
            notebook = nbformat.read(root / name, as_version=4)
            for i, cell in enumerate(notebook.cells):
                if cell.cell_type == 'code':
                    with self.subTest(notebook=name, cell=i):
                        compile(TransformerManager().transform_cell(cell.source), f'{name}:cell {i}', 'exec')

    def test_setup_recovers_stale_imports_and_rebinds_functions(self):
        root = Path(__file__).resolve().parents[1]
        # Isolate IPython and module reloading from the analysis test process.
        script = r'''
import nbformat, sys
from IPython.terminal.interactiveshell import TerminalInteractiveShell
import bio_visuals.data as data

notebook = nbformat.read(sys.argv[1], as_version=4)
setup = next(cell.source for cell in notebook.cells if cell.cell_type == 'code')
shell = TerminalInteractiveShell.instance()

# Reproduce the user's old cached module: the new comparison module cannot import.
del data.ROI_PAIRS
try:
    import bio_visuals.comparison
except ImportError as error:
    assert 'ROI_PAIRS' in str(error), error
else:
    raise AssertionError('The stale-import reproduction must fail before setup runs')

for _ in range(2):
    data.load_data = lambda *args, **kwargs: 'stale implementation'
    shell.run_cell(setup).raise_error()
    import bio_visuals.comparison as comparison
    assert len(data.ROI_PAIRS) == 4
    assert comparison.ROI_PAIRS is data.ROI_PAIRS
    assert comparison.summarize is data.summarize
    assert shell.user_ns['load_data'] is data.load_data
    assert data.load_data.__name__ == 'load_data'
    assert shell.user_ns['roi_summary_figure'] is comparison.roi_summary_figure
    del data.ROI_PAIRS
'''
        for notebook in ('PLA_ChAT_visual_analysis.ipynb', 'PLA_ChAT_scientific_analysis.ipynb'):
            with self.subTest(notebook=notebook):
                result = subprocess.run(
                    [sys.executable, '-c', script, notebook], cwd=root, capture_output=True,
                    text=True, timeout=120,
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == '__main__':
    unittest.main()

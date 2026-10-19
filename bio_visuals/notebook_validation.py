"""Checks for portable saved notebook results, including headless figure output."""
import ast
import base64
import re


def _shows_static_figure(source):
    """Find executed show()/Image() calls, excluding helper definitions."""
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return False  # IPython magics are validated by the execution engine.
    def visits(node):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            return False
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            if node.func.id in ('show', 'Image'):
                return True
        return any(visits(child) for child in ast.iter_child_nodes(node))
    return visits(tree)


def validate_saved_outputs(notebook):
    """Reject incomplete runs, errors and missing or malformed static PNG output.

    A rendered figure is required in every cell that uses our static display
    helpers. Interactive Plotly/HTML outputs may use their own MIME formats.
    Returns portable-output counts for the validation record.
    """
    n_code = n_png = 0
    for index, cell in enumerate(notebook.cells):
        if cell.cell_type != 'code':
            continue
        n_code += 1
        if cell.execution_count is None:
            raise ValueError(f'Cell {index} was not executed')
        png_in_cell = 0
        for output in cell.outputs:
            if output.output_type == 'error':
                raise ValueError(f'Cell {index} contains an execution error')
            data = output.get('data', {})
            plain = data.get('text/plain', '')
            if re.search(r'<Figure size .* with \d+ Axes>', plain) and not any(
                    mime in data for mime in ('image/png', 'image/svg+xml')):
                raise ValueError(f'Cell {index} contains an unrendered Matplotlib figure')
            if 'image/png' in data:
                try:
                    png = base64.b64decode(''.join(data['image/png'].split()), validate=True)
                except (ValueError, TypeError) as error:
                    raise ValueError(f'Cell {index} has invalid PNG encoding') from error
                if not png.startswith(b'\x89PNG\r\n\x1a\n') or len(png) < 32:
                    raise ValueError(f'Cell {index} has invalid PNG bytes')
                png_in_cell += 1
        if (_shows_static_figure(cell.source) and not png_in_cell
                and not cell.metadata.get('static_figure_optional', False)):
            raise ValueError(f'Cell {index} should show a static figure but has no embedded PNG')
        n_png += png_in_cell
    return {'code_cells': n_code, 'embedded_pngs': n_png}

"""Execute the handoff notebook with this Python environment and save outputs."""
from pathlib import Path
import sys
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
ROOT=Path(__file__).resolve().parents[1]
path=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'PLA_ChAT_spatial_handoff_analysis.ipynb')
if path.resolve().parent != ROOT or path.suffix != '.ipynb':
    raise ValueError('Choose a notebook in the repository root')
nb=nbformat.read(path,as_version=4)
km=KernelManager(kernel_name='python3')
km.kernel_spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}']
def progress(cell,cell_index,**kwargs):
    if cell.cell_type=='code':print(f'Executing cell {cell_index}: {cell.source.splitlines()[0][:100]}',flush=True)
client=NotebookClient(nb,timeout=900,km=km,resources={'metadata':{'path':str(ROOT)}},on_cell_start=progress)
try:
    client.execute()
    nbformat.write(nb,path)
    errors=[o for c in nb.cells if c.cell_type=='code' for o in c.outputs if o.output_type=='error']
    assert not errors,errors
    print(f'Saved executed notebook: {path} ({path.stat().st_size/1e6:.1f} MB)',flush=True)
finally:
    if km.has_kernel:km.shutdown_kernel(now=True)

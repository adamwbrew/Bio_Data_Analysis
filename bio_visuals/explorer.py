"""Self-contained, notebook-embedded exploration of measured 2D sections.

The comparison stack uses an ordinal display axis, never an anatomical z axis.
All filtering is per section; no animal or PLA-arm averaging occurs in the UI.
"""
from pathlib import Path
import json
import uuid

import numpy as np

from .handoff import ANIMAL_ORDER, ARMS, polygon_centroid, review_status
from .handoff_plots import ANIMAL_COLORS


def _grid_payload(grid, fine=False):
    keys = ['n_pla', 'n_chat', 'n_positive']
    if not fine:
        keys += ['n_all', 'pla_area', 'n_large', 'overlap_area']
    result = {'x': ((grid['x_edges'][1:] + grid['x_edges'][:-1]) / 2).tolist(),
              'y': ((grid['y_edges'][1:] + grid['y_edges'][:-1]) / 2).tolist(),
              'width': float(grid['width'])}
    for key in keys:
        a = np.asarray(grid[key])
        if not np.isfinite(a).all() or (a < 0).any():
            raise ValueError(f'Explorer requires finite nonnegative {key}')
        result[key] = (a.astype(np.int64) if key.startswith('n_') else a).tolist()
    return result


def explorer_payload(fine, geometry, surfaces, metrics, metadata, reviews):
    """Keep measured grids separate and disclose each coordinate transform."""
    from plotly.colors import get_colorscale
    by_job = {r['job']: r for r in reviews}
    ready = review_status(reviews).set_index('job').ready.to_dict()
    sections = []
    for arm in ARMS.values():
        rows = metadata[metadata.pla_arm == arm].set_index('animal_id').loc[ANIMAL_ORDER]
        for order, (animal, row) in enumerate(rows.iterrows()):
            job = str(row.job)
            r = by_job[job]
            candidate = r.get('candidate')
            approved = bool(ready[job])
            sections.append(dict(
                job=job, animal=animal, group=row.running_group, arm=arm, order=order,
                color=ANIMAL_COLORS[animal], source=row.source_file,
                candidate_center=candidate['centroid_native'] if candidate else None,
                candidate_polygon=candidate['polygon_native'] if candidate else None,
                reviewed_center=polygon_centroid(r['reviewed_polygon_native']).tolist() if approved else None,
                reviewed_polygon=r['reviewed_polygon_native'] if approved else None,
                reviewed_flip=r['orientation'] if approved else None,
                ready=approved,
                fine=_grid_payload(fine[job], fine=True),
                grids={str(int(g[job]['width'])): _grid_payload(g[job])
                       for g in (geometry, surfaces, metrics)},
            ))
    return dict(sections=sections, arms=list(ARMS.values()), animal_order=ANIMAL_ORDER,
                colorscales={name:get_colorscale(name) for name in ('Magma','Viridis','Cividis')},
                order_meaning='Home cage then Running; display order only, not serial-section order',
                units='HALO export units; physical calibration unverified',
                version='1.1.0')


def explorer_html(fine, geometry, surfaces, metrics, metadata, reviews):
    """Return an offline HTML output embedded in a trusted Jupyter notebook.

    No server, CDN, widget extension, or live kernel is needed after rendering.
    GitHub suppresses JavaScript; the notebook supplies static companion figures.
    """
    from plotly.offline import get_plotlyjs
    payload = explorer_payload(fine, geometry, surfaces, metrics, metadata, reviews)
    data = json.dumps(payload, separators=(',', ':'), ensure_ascii=False, allow_nan=False)
    data = data.replace('<', '\\u003c')
    template = Path(__file__).with_name('explorer.html').read_text()
    return template.replace('__ROOT_ID__', 'halo-' + uuid.uuid4().hex).replace('__DATA__', data).replace(
        '__PLOTLY__', get_plotlyjs())

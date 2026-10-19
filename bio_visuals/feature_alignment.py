"""Inspectable cavity-and-shape alignment proposals, never anatomical approvals.

Translation is fixed by an annotated internal-cavity polygon. Discrete flips are
ranked against one reference using an adjacent internal structure and trimmed boundary distance. Missing
outer tissue can increase residuals but cannot move the anchor or stretch data.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from .handoff import polygon_centroid, transform_xy, ANIMAL_ORDER, ARMS

FLIPS = ('none', 'horizontal', 'vertical', 'both')


def boundary_points(grid):
    occupied = grid['n_all'] > 0
    # One-bin closing connects small sampling gaps; never fill large cavities.
    mask = ndi.binary_closing(np.pad(occupied, 2), iterations=1)[2:-2, 2:-2]
    labels, _ = ndi.label(mask)
    sizes = np.bincount(labels.ravel()); sizes[0] = 0
    mask = labels == sizes.argmax() if sizes.max() else occupied
    boundary = mask & ~ndi.binary_erosion(mask)
    y, x = np.where(boundary)
    return np.column_stack(((grid['x_edges'][x] + grid['x_edges'][x + 1]) / 2,
                            (grid['y_edges'][y] + grid['y_edges'][y + 1]) / 2))


def boundary_score(a, b, retained=.8):
    """Average the closest 80% in each direction; also retain untrimmed p90.

    This tolerates a minority of unmatched/cropped edges. It is a descriptive
    shape residual, NOT independent target-registration error or validation.
    """
    if not len(a) or not len(b):
        raise ValueError('Nonempty boundaries required')
    ab = cKDTree(b).query(a)[0]; ba = cKDTree(a).query(b)[0]
    trimmed = [np.sort(d)[:max(1, int(np.ceil(retained * len(d))))].mean() for d in (ab, ba)]
    return float(np.mean(trimmed)), float(max(np.quantile(ab, .9), np.quantile(ba, .9)))


def allowed_flips(center, reference_center, cue=None, reference_cue=None):
    """Match the visible cue's above/below side; never alter translation/scale.

    Horizontal handedness is still ranked by geometry. If the cue is absent or
    coincides with the cavity y, retain all four candidates and disclose that.
    """
    if cue is None or reference_cue is None:
        return FLIPS
    dy=float(cue[1])-float(center[1]); refdy=float(reference_cue[1])-float(reference_center[1])
    if dy == 0 or refdy == 0:
        return FLIPS
    same=dy*refdy>0
    return ('none','horizontal') if same else ('vertical','both')


def alignment_proposals(geometry, metadata, path=None):
    path = Path(path or Path(__file__).resolve().parents[1] / 'metadata/handoff_landmark_proposals.json')
    annotations = {r['job']:r for r in json.loads(path.read_text())}
    result, scores = [], []
    # 21F is a fixed, identifiable cavity reference, never selected by condition contrast.
    for arm in ARMS.values():
        rows = metadata[metadata.pla_arm.eq(arm)].set_index('animal_id').loc[ANIMAL_ORDER]
        ref_job = str(rows.loc['21F', 'job'])
        ref_center = polygon_centroid(annotations[ref_job]['polygon_native'])
        ref_cue=annotations[ref_job].get('orientation_cue_native')
        ref = transform_xy(boundary_points(geometry[ref_job]), ref_center, 'none')[1]
        ref_native = boundary_points(geometry[ref_job])
        for animal, row in rows.iterrows():
            job = str(row.job); annotation = annotations[job]
            if annotation['source_sha256'] != row.sha256 or annotation['animal_id'] != animal or annotation['pla_arm'] != arm:
                raise ValueError(f'Landmark annotation identity/hash mismatch: {job}')
            center = polygon_centroid(annotation['polygon_native'])
            points = boundary_points(geometry[job])
            allowed=allowed_flips(center,ref_center,annotation.get('orientation_cue_native'),ref_cue)
            options = []
            for flip in FLIPS:
                score, p90 = boundary_score(transform_xy(points, center, flip)[1], ref)
                options.append((score, flip, p90))
                scores.append(dict(job=job, animal_id=animal, pla_arm=arm, flip=flip,
                                   trimmed_boundary_distance=score, boundary_p90=p90, cue_compatible=flip in allowed))
            options.sort(key=lambda v:(v[0], FLIPS.index(v[1])))
            boundary_best=options[0][1]
            options=[v for v in options if v[1] in allowed]
            best, flip, p90 = options[0]
            margin = (options[1][0] - best) / max(options[1][0], 1.)
            native_score, _ = boundary_score(points, ref_native)
            center_score, _ = boundary_score(transform_xy(points, center, 'none')[1], ref)
            result.append(dict(**annotation, center=center.tolist(), flip=flip,
                               reference_job=ref_job, reference_animal='21F',
                               boundary_best_flip=boundary_best, cue_boundary_conflict=flip!=boundary_best,
                               orientation_cue_used=len(allowed)==2,
                               native_boundary_distance=native_score, centered_boundary_distance=center_score,
                               aligned_boundary_distance=best, boundary_p90=p90,
                               flip_margin_fraction=margin,
                               flip_ambiguous=bool(margin < .1 and job != ref_job),
                               method='Cavity polygon centroid; adjacent-structure above/below cue where visible; boundary ranks remaining flips; no rotation, scale or warp',
                               approved=False))
    return result, pd.DataFrame(scores)


def alignment_figure(geometry, proposals):
    import matplotlib.pyplot as plt
    from .handoff_plots import ANIMAL_COLORS
    fig, axes = plt.subplots(2, 3, figsize=(17, 12))
    for row, arm in enumerate(ARMS.values()):
        part = [p for p in proposals if p['pla_arm'] == arm]
        for col, mode in enumerate(['Native image coordinates', 'Cavity-centered only', 'Cavity + internal cue + discrete flip']):
            ax = axes[row, col]
            for p in part:
                points = boundary_points(geometry[p['job']])
                if col:
                    points = transform_xy(points, p['center'], p['flip'] if col == 2 else 'none')[1]
                ax.scatter(points[:,0], points[:,1], s=2, alpha=.5, color=ANIMAL_COLORS[p['animal_id']], label=p['animal_id'])
            ax.set_title(f'{arm}\n{mode}', fontsize=12)
            ax.set_aspect('equal'); ax.invert_yaxis(); ax.set_xlabel('x · export units'); ax.set_ylabel('y · export units')
            if col: ax.plot(0,0,'k+',ms=13,mew=2)
            if row == 0 and col == 0: ax.legend(ncol=3, fontsize=9, markerscale=3, loc='upper center')
    fig.suptitle('Align internal features first; keep missing tissue missing', fontsize=21, fontweight='bold', y=.99)
    fig.text(.05,.015,'All 12 sections have inspectable proposals. Cross = cavity anchor. Flips are geometric suggestions, not anatomical approvals.\nCropped edges do not move the anchor; remaining shape differences may reflect section level, tissue loss or imperfect landmarks.', fontsize=11)
    fig.tight_layout(rect=(0,.065,1,.95), h_pad=3.0)
    return fig


def landmark_figure(fine, proposals, arm):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator, StrMethodFormatter
    fig, axes = plt.subplots(2,3,figsize=(16,10))
    for ax,p in zip(axes.flat,[p for p in proposals if p['pla_arm']==arm]):
        g=fine[p['job']];poly=np.asarray(p['polygon_native']);center=np.asarray(p['center'])
        ax.pcolormesh(g['x_edges'],g['y_edges'],np.log1p(g['n_all']),cmap='Greys',vmin=0,vmax=5,rasterized=True)
        loop=np.vstack([poly,poly[0]])
        ax.plot(loop[:,0],loop[:,1],color='#d06738',lw=1.7);ax.plot(*center,'+',color='#d06738',ms=12,mew=2)
        cue=p.get('orientation_cue_native')
        if cue is not None:
            ax.plot(*cue,'o',ms=10,mfc='none',mec='#009dc0',mew=2)
            ax.plot([center[0],cue[0]],[center[1],cue[1]],':',color='#009dc0',lw=1.2)
        ax.set(xlim=(g['x_edges'][0],g['x_edges'][-1]),ylim=(g['y_edges'][-1],g['y_edges'][0]),aspect='equal',
               title=f'{p["animal_id"]} · proposed {p["flip"]} flip'+ ('\nBoundary uncertain' if p['ambiguous_boundary'] else '\nCavity trace proposal'))
        ax.set_xlabel('Native x · export units');ax.set_ylabel('Native y · export units')
        for axis in (ax.xaxis,ax.yaxis):
            axis.set_major_locator(MaxNLocator(nbins=4));axis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    fig.suptitle(f'{arm} | inspect the actual proposed cavity traces',fontsize=20,fontweight='bold')
    fig.text(.05,.015,'100-unit all-object reconstruction. Orange = cavity trace; blue circle = adjacent orientation cue, where identifiable. No original microscopy pixels were available.\nCavity and cue are computational visual proposals, not confirmed homologous anatomy. They do not change the scientist review file.',fontsize=10)
    fig.tight_layout(rect=(0,.065,1,.95));return fig

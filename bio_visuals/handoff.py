"""Raw-data, animal-aware analysis for the October 2026 HALO handoff.

Native descriptive maps are explicitly distinct from reviewed registration.
No atlas, continuous rotation, inferred calibration, or anatomical z is used.
"""
from pathlib import Path
import hashlib
import itertools
import json
import platform
import subprocess
from datetime import datetime, timezone
import importlib.metadata

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from matplotlib.path import Path as MplPath

from .entire_dataset import run_summary, _sha256

VERSION = '1.2.0'
GROUPS = {'22F': 'Running', '21M': 'Running', '22M': 'Running',
          '21F': 'Home cage', '23M': 'Home cage', '24M': 'Home cage'}
ANIMAL_ORDER = ['21F', '23M', '24M', '22F', '21M', '22M']
ARMS = {'Hk1': 'TDP43–HK1', 'Pfkp': 'TDP43–PFKP'}
SIZE_EDGES = [0, .03, .1, .3, 1, np.inf]
SIZE_LABELS = ['<0.03', '0.03–<0.1', '0.1–<0.3', '0.3–<1', '≥1']
PARAMETERS = dict(protocol_version=VERSION, reconstruction_bin=100., geometry_bin=250.,
                  surface_bin=500., surface_height='log(1 + PLA count per full bin); not anatomical depth',
                  metric_bin=1000., min_pla_per_bin=10, large_area_um2=1., small_area_um2=.1,
                  size_edges_um2=[0, .03, .1, .3, 1, 'infinity'],
                  edge_buffer_export_units=500., smoothing='none',
                  support='occupied all-object geometry bins only; provisional sampling support, not a tissue mask',
                  candidate_void='largest enclosed hole after 2 iterations of 3x3-cross binary closing on occupied 250-unit bins; >=6 bins',
                  orientation='none for preview; all choices require review',
                  coordinate_units='HALO export units; physical scale unverified',
                  display_exception='Native descriptive maps and whole-section summaries may be inspected before registration QC. No anatomical condition maps are inferred.',
                  min_animals_per_group_bin=3,
                  exclusions='4978 superseded by 5403; Initial Data copies not concatenated; no primary object filtering')
PARAMETERS.update(interactive_views='3D measurements, 100-unit channel overlay, flat layer overlay and ordinal exploded stack',
                  layer_order='21F,23M,24M,22F,21M,22M; display order, not anatomical section order',
                  candidate_layer_scope='Geometry only; unreviewed centroid translation, no flips',
                  explorer_bins=[250,500,1000], explorer_default_bin=500,
                  explorer_threshold='Display control 1–50 PLA per metric bin, default 10; no raw rows removed',
                  condition_difference='Home cage minus Running; equal animal means per arm',
                  size_standardization='Equal-animal reference size proportions; require 10 PLA and observed association in every animal per included stratum',
                  exact_enumeration='20 allocations of six animals; design diagnostic only, exchangeability unverified',
                  changelog='1.2.0: condition averages and signed contrasts, all nine deletion checks, common-size standardization, image-layout audit and interactive comparison. Native data and anatomical approvals unchanged.')


def section_metadata(manifest):
    """Animal assignments come ONLY from handoff p1/p8, not M/F heuristics."""
    m = manifest.loc[manifest.is_primary].copy()
    if set(m.condition.astype(str)) != set(GROUPS):
        raise ValueError('Unexpected animal labels: review the handoff assignments')
    m['animal_id'] = m.condition.astype(str)
    m['running_group'] = m.animal_id.map(GROUPS)
    m['pla_arm'] = m.pla_type.astype(str).map(ARMS)
    m['marker_panel'] = 'ChAT'
    m['section_id'] = m.animal_id + '_' + m.pla_type.astype(str)
    m['assignment_source'] = 'HALO_spatial_analysis_raw_data_handoff.pdf pp1,8'
    if m.pla_arm.isna().any() or m.section_id.duplicated().any() or len(m) != 12:
        raise ValueError('Expected exactly six animals × two distinct PLA sections')
    return m.sort_values(['pla_type', 'animal_id']).reset_index(drop=True)


def animal_summary(objects, metadata):
    result = run_summary(objects)
    result = result.merge(metadata[['job', 'animal_id', 'running_group', 'pla_arm', 'section_id']],
                          on='job', validate='one_to_one')
    extra = []
    for job, frame in objects.groupby('job', observed=True):
        p = frame[frame.object_type == 'PLA']
        c = frame[frame.object_type == 'ChAT']
        extra.append(dict(job=str(job), pla_total_area_um2=p.area.sum(),
                          chat_total_object_area_um2=c.area.sum(),
                          marker_positive_pla_area_um2=p.loc[p.chat_present, 'area'].sum(),
                          actual_overlap_area_um2=p.overlap_area.sum(),
                          pla_median_area_um2=p.area.median(),
                          large_pla_pct=100*p.area.ge(1).mean(),
                          marker_positive_area_share_pct=100*p.loc[p.chat_present, 'area'].sum()/p.area.sum()))
    return result.merge(pd.DataFrame(extra), on='job', validate='one_to_one')


def group_contrasts(summary, metric='pla_association_pct'):
    """Equal-animal means; delete-one-from-each-group sensitivity is not a CI."""
    rows = []
    for arm, g in summary.groupby('pla_arm', observed=True):
        if g.animal_id.duplicated().any():
            raise ValueError('First reduce repeated sections to one value per animal and arm')
        r = g.loc[g.running_group == 'Running', metric].to_numpy(float)
        h = g.loc[g.running_group == 'Home cage', metric].to_numpy(float)
        if len(r) != 3 or len(h) != 3 or not np.isfinite(np.r_[r, h]).all():
            raise ValueError('This protocol expects three finite animal values per condition')
        loo = [np.delete(r,i).mean()-np.delete(h,j).mean() for i,j in itertools.product(range(3),repeat=2)]
        rows.append(dict(pla_arm=arm, metric=metric, n_home=3, n_running=3,
                         home_mean=h.mean(), running_mean=r.mean(), running_minus_home=r.mean()-h.mean(),
                         leave_pair_out_min=min(loo), leave_pair_out_max=max(loo)))
    return pd.DataFrame(rows)


def size_tables(objects, metadata):
    rows, sensitivity = [], []
    for job, frame in objects.groupby('job', observed=True):
        p = frame[frame.object_type == 'PLA'].copy()
        p['size_bin'] = pd.cut(p.area, SIZE_EDGES, labels=SIZE_LABELS, right=False)
        for label, g in p.groupby('size_bin', observed=False):
            rows.append(dict(job=str(job), size_bin=str(label), n_pla=len(g),
                             n_positive=int(g.chat_present.sum()),
                             fraction_pct=100*g.chat_present.mean() if len(g) else np.nan))
        for cutoff in (0., .03, .1, .3, 1.):
            q = p[p.area >= cutoff]
            sensitivity.append(dict(job=str(job), min_area_um2=cutoff, n_pla=len(q),
                                    fraction_pct=100*q.chat_present.mean() if len(q) else np.nan))
    keys = metadata[['job', 'animal_id', 'running_group', 'pla_arm']]
    return (pd.DataFrame(rows).merge(keys,on='job',validate='many_to_one'),
            pd.DataFrame(sensitivity).merge(keys,on='job',validate='many_to_one'))


def edges_for(frame, width, padding=0):
    if not np.isfinite(width) or width <= 0:
        raise ValueError('Bin width must be positive')
    def one(lo,hi):
        start=(np.floor(lo/width)-padding)*width
        end=(np.ceil(hi/width)+padding)*width
        return start+np.arange(int(round((end-start)/width))+1)*width
    return one(frame.xmin.min(),frame.xmax.max()),one(frame.ymin.min(),frame.ymax.max())


def histogram(frame, xe, ye, weights=None):
    return np.histogram2d(frame.y.to_numpy(), frame.x.to_numpy(), bins=(ye,xe), weights=weights)[0]


def build_grid(frame, width=1000., min_pla=10):
    """Keep raw counts. A bin without observations is unknown, never zero tissue."""
    xe,ye=edges_for(frame,width)
    p=frame[frame.object_type=='PLA'];c=frame[frame.object_type=='ChAT']
    n=histogram(p,xe,ye); k=histogram(p[p.chat_present],xe,ye)
    chat=histogram(c,xe,ye); all_n=n+chat
    def ratio(a,b,minimum):
        return np.divide(a,b,out=np.full_like(a,np.nan),where=b>=minimum)
    area=histogram(p,xe,ye,p.area.to_numpy())
    associated_area=histogram(p[p.chat_present],xe,ye,p.loc[p.chat_present,'area'].to_numpy())
    overlap_area=histogram(p,xe,ye,p.overlap_area.to_numpy())
    large=histogram(p[p.area>=1],xe,ye)
    small=histogram(p[p.area<.1],xe,ye)
    support=all_n>0
    fraction=100*ratio(k,n,min_pla)
    # Descriptive counts per full bin, NOT density per measured tissue area.
    return dict(x_edges=xe,y_edges=ye,width=width,n_pla=n,n_chat=chat,n_all=all_n,
                n_positive=k,support=support,fraction=fraction,
                pla_area=area,positive_pla_area=associated_area,overlap_area=overlap_area,
                mean_pla_area=ratio(area,n,min_pla),large_fraction=100*ratio(large,n,min_pla),
                large_small_ratio=ratio(large,small,min_pla),n_large=large,n_small=small,
                min_pla=min_pla)


def make_grids(objects,width=1000.,min_pla=10):
    return {str(job):build_grid(g,width,min_pla) for job,g in objects.groupby('job',observed=True)}


def polygon_centroid(vertices):
    v=np.asarray(vertices,float)
    if v.ndim!=2 or v.shape[1]!=2 or len(v)<3 or not np.isfinite(v).all():
        raise ValueError('A polygon requires at least three finite x,y vertices')
    x,y=v.T;cross=x*np.roll(y,-1)-np.roll(x,-1)*y
    signed=cross.sum()
    if abs(signed)<1e-10:
        raise ValueError('Degenerate landmark polygon')
    centroid=np.array([((x+np.roll(x,-1))*cross).sum(),((y+np.roll(y,-1))*cross).sum()])/(3*signed)
    if not MplPath(np.vstack([v,v[0]])).contains_point(centroid):
        raise ValueError('Polygon centroid falls outside cavity: review this trace')
    return centroid


def candidate_void(grid):
    """Unreviewed *geometric* proposal, not an anatomical identification."""
    import matplotlib.pyplot as plt
    occupied=grid['n_all']>0
    closed=ndi.binary_closing(occupied,iterations=2)
    holes=ndi.binary_fill_holes(closed)&~closed
    labels,num=ndi.label(holes)
    sizes=np.bincount(labels.ravel());sizes[0]=0
    if num==0 or sizes.max()<6:
        return None
    target=labels==int(sizes.argmax())
    xc=(grid['x_edges'][:-1]+grid['x_edges'][1:])/2
    yc=(grid['y_edges'][:-1]+grid['y_edges'][1:])/2
    fig,ax=plt.subplots()
    contours=ax.contour(xc,yc,target.astype(float),levels=[.5])
    paths=contours.allsegs[0]
    plt.close(fig)
    if not paths:
        return None
    poly=max(paths,key=len)
    try:
        center=polygon_centroid(poly)
    except ValueError:
        return None
    return dict(polygon_native=poly.tolist(),centroid_native=center.tolist(),
                hole_bins=int(sizes.max()),status='unreviewed geometric candidate')


def default_reviews(metadata,geometry):
    records=[]
    for row in metadata.itertuples():
        candidate=candidate_void(geometry[str(row.job)])
        records.append(dict(job=str(row.job),section_id=row.section_id,animal_id=row.animal_id,
                            pla_arm=row.pla_arm,running_group=row.running_group,
                            source_sha256=row.sha256,candidate=candidate,
                            reviewed_polygon_native=None,landmark_reviewed=False,
                            orientation='none',orientation_reviewed=False,
                            scale_common_confirmed=False,coordinate_units='export units',
                            support_reviewed=False,support_qc='uncertain',
                            clipping_note='Expected straight cut edge is not automatically tissue loss; inspect source image/ROI.',
                            reviewer='',second_reviewer='',review_notes='',
                            protocol_version=VERSION))
    return records


def load_reviews(path,metadata,geometry):
    if not Path(path).exists():
        return default_reviews(metadata,geometry)
    records=json.loads(Path(path).read_text())
    expected=metadata.set_index('job')
    if set(r['job'] for r in records)!=set(expected.index) or len(records)!=len(expected):
        raise ValueError('Review file does not match the 12 primary sections')
    for r in records:
        source=expected.loc[r['job']]
        for field in ('section_id','animal_id','pla_arm','running_group','source_sha256'):
            truth=source['sha256'] if field=='source_sha256' else source[field]
            if r[field]!=truth:
                raise ValueError(f'Review identity/hash mismatch: {r["job"]} {field}')
    return records


def review_status(records):
    rows=[]
    for r in records:
        reasons=[]
        for flag in ('landmark_reviewed','orientation_reviewed','scale_common_confirmed','support_reviewed'):
            if not isinstance(r.get(flag),bool):
                raise ValueError(f'{flag} must be a JSON boolean, not text')
        if not r.get('landmark_reviewed') or not r.get('reviewed_polygon_native'):
            reasons.append('ventricular trace')
        else:
            polygon_centroid(r['reviewed_polygon_native'])
        if r.get('orientation') not in ('none','horizontal','vertical','both'):
            raise ValueError('Only discrete flips are permitted')
        if not r.get('orientation_reviewed'):reasons.append('orientation')
        if not r.get('scale_common_confirmed'):reasons.append('common scale')
        if not r.get('support_reviewed') or r.get('support_qc') not in ('complete','minor loss'):
            reasons.append('support/clipping')
        if not r.get('reviewer','').strip():reasons.append('reviewer')
        if r.get('orientation')!='none' and not r.get('second_reviewer','').strip():
            reasons.append('second reviewer for flip')
        rows.append(dict(job=r['job'],animal_id=r['animal_id'],pla_arm=r['pla_arm'],
                         candidate_found=r.get('candidate') is not None,
                         ready=not reasons,needs_review=', '.join(reasons)))
    return pd.DataFrame(rows)


def transform_xy(xy,centroid,orientation='none'):
    """Reflect about native origin, then translate the reflected landmark."""
    if orientation not in ('none','horizontal','vertical','both'):
        raise ValueError('Only none, horizontal, vertical, both allowed')
    signs=np.array([-1 if orientation in ('horizontal','both') else 1,
                    -1 if orientation in ('vertical','both') else 1])
    oriented=np.asarray(xy,float)*signs
    center=np.asarray(centroid,float)*signs
    return oriented,oriented-center


def provisional_coverage(geometry,records,width=1000.):
    """QC-only: candidate-centered support; excludes sections with no candidate.

    Does not produce a biological metric or claim reviewed registration.
    Each occupied geometry bin contributes once to an animal's coarse support.
    """
    xe=np.arange(-15000,16001,width);ye=np.arange(-15000,16001,width)
    maps={};used=[]
    for r in records:
        if not r.get('candidate'):continue
        grid=geometry[r['job']]
        yy,xx=np.where(grid['support'])
        xy=np.c_[grid['x_edges'][xx]+grid['width']/2,grid['y_edges'][yy]+grid['width']/2]
        _,aligned=transform_xy(xy,r['candidate']['centroid_native'],r['orientation'])
        h=np.histogram2d(aligned[:,1],aligned[:,0],bins=(ye,xe))[0]
        if int(h.sum()) != len(aligned):raise ValueError('Coverage extent clips candidate support')
        maps[r['job']]=h>0;used.append(r)
    counts={}
    for arm,group in itertools.product(ARMS.values(),('Home cage','Running')):
        selected=[maps[r['job']] for r in used if r['pla_arm']==arm and r['running_group']==group]
        counts[(arm,group)]=np.sum(selected,axis=0) if selected else np.zeros((len(ye)-1,len(xe)-1))
    return dict(x_edges=xe,y_edges=ye,counts=counts,used=used,maps=maps)


def equal_animal_maps(values,minimum=3):
    """No object pooling; missing/unsupported bins are NaN."""
    stack=np.asarray(values,float)
    n=np.isfinite(stack).sum(axis=0)
    means=np.divide(np.nansum(stack,axis=0),n,out=np.full(n.shape,np.nan),where=n>=minimum)
    return means,n


def registered_fraction_maps(objects,records,width=1000.,min_pla=10):
    """Strict gate for future reviewed maps, keeping arms/animals independent.

    Registration remains in verified common export units (no rescaling). Support
    is the explicitly accepted all-object occupied-bin proxy. Users needing true
    tissue-area density must supply ROI polygons and a new protocol version.
    """
    status=review_status(records)
    if not status.ready.all():
        raise ValueError('Group maps withheld: ventricular, orientation, scale and support review must pass for all sections')
    # Freeze to a common grid only after reviewers confirm the native scale.
    xe=np.arange(-30000,30001,width);ye=xe.copy();animal_maps={}
    for r in records:
        g=objects[objects.job==r['job']]
        center=polygon_centroid(r['reviewed_polygon_native'])
        _,xy=transform_xy(g[['x','y']].to_numpy(),center,r['orientation'])
        p=g.object_type.eq('PLA').to_numpy();k=p&g.chat_present.to_numpy()
        hist=lambda z:np.histogram2d(xy[z,1],xy[z,0],bins=(ye,xe))[0]
        n=hist(p);positive=hist(k)
        if n.sum()!=p.sum():raise ValueError('Registered grid clips PLA objects')
        animal_maps[r['job']]=100*np.divide(positive,n,out=np.full_like(n,np.nan),where=n>=min_pla)
    result={}
    for arm in ARMS.values():
        result[arm]={}
        for group in ('Home cage','Running'):
            rows=[r for r in records if r['pla_arm']==arm and r['running_group']==group]
            if len(rows)!=3 or len({r['animal_id'] for r in rows})!=3:
                raise ValueError('Expected exactly three distinct animals in each arm/group')
            mean,n=equal_animal_maps([animal_maps[r['job']] for r in rows],3)
            result[arm][group]=dict(mean=mean,n=n)
        result[arm]['difference']=result[arm]['Running']['mean']-result[arm]['Home cage']['mean']
    return dict(x_edges=xe,y_edges=ye,arms=result,animal_maps=animal_maps)


def proximity_analysis(objects,metadata,geometry,buffer=500.):
    """Within-section nearest midpoint distances + conservative border sensitivity.

    Occupied-bin proxy edges include internal gaps; the exclusion is deliberately
    conservative and does not replace reviewed tissue boundaries. No null test.
    """
    curves=[];stats=[];quantiles=np.linspace(0,1,201)
    for job,g in objects.groupby('job',observed=True):
        job=str(job);p=g[g.object_type=='PLA'];c=g[g.object_type=='ChAT']
        d=cKDTree(c[['x','y']].to_numpy()).query(p[['x','y']].to_numpy())[0]
        grid=geometry[job];support=grid['support']
        padded=np.pad(support,1)
        clearance=ndi.distance_transform_edt(padded)[1:-1,1:-1]*grid['width']-grid['width']*np.sqrt(2)/2
        ix=np.searchsorted(grid['x_edges'],p.x,side='right')-1
        iy=np.searchsorted(grid['y_edges'],p.y,side='right')-1
        ix=np.clip(ix,0,support.shape[1]-1);iy=np.clip(iy,0,support.shape[0]-1)
        edge_distance=clearance[iy,ix]
        # NN must also lie closer than the conservative distance to unsampled support.
        retained=(edge_distance>=buffer)&(d<=edge_distance)
        stats.append(dict(job=job,n_pla=len(p),n_interior=int(retained.sum()),
                          retained_pct=100*retained.mean(),median_all=np.median(d),
                          median_interior=np.median(d[retained]) if retained.any() else np.nan,
                          q95_all=np.quantile(d,.95)))
        for label,values in [('All PLA',d),('Interior sensitivity',d[retained])]:
            if len(values):
                curves.extend(dict(job=job,selection=label,probability=float(q),distance=float(x))
                              for q,x in zip(quantiles,np.quantile(values,quantiles)))
    keys=metadata[['job','animal_id','pla_arm','running_group']]
    return pd.DataFrame(curves).merge(keys,on='job'),pd.DataFrame(stats).merge(keys,on='job')


def spatial_review_queue(grids,metadata):
    rows=[]
    for r in metadata.itertuples():
        g=grids[str(r.job)];valid=g['n_pla']>=g['min_pla']
        for reason,key in [('largest PLA count','n_pla'),('highest supported ChAT fraction','fraction')]:
            z=np.where(valid,g[key],np.nan)
            if not np.isfinite(z).any():continue
            iy,ix=np.unravel_index(np.nanargmax(z),z.shape)
            rows.append(dict(job=str(r.job),animal_id=r.animal_id,pla_arm=r.pla_arm,
                             reason=reason,xmin=g['x_edges'][ix],xmax=g['x_edges'][ix+1],
                             ymin=g['y_edges'][iy],ymax=g['y_edges'][iy+1],
                             n_pla=int(g['n_pla'][iy,ix]),n_chat=int(g['n_chat'][iy,ix]),
                             positive_pct=g['fraction'][iy,ix],source_file=r.source_file))
    return pd.DataFrame(rows)


def export_objects(objects,metadata,reviews,directory):
    """Optional lossless provenance-rich per-section CSV.gz, no fake alignment."""
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    lookup={r['job']:r for r in reviews};status=review_status(reviews).set_index('job')
    for m in metadata.itertuples():
        g=objects[objects.job==m.job].copy();r=lookup[str(m.job)]
        for field in ('animal_id','running_group','pla_arm','source_file','section_id','algorithm'):
            g[field]=getattr(m,field)
        g=g.rename(columns={'x':'x_native','y':'y_native'})
        g['orientation_transform']=r['orientation'];g['protocol_version']=VERSION
        for field in ('x_oriented','y_oriented','x_ventricle','y_ventricle','x_aligned','y_aligned'):
            g[field]=np.nan
        if bool(status.loc[str(m.job),'ready']):
            center=polygon_centroid(r['reviewed_polygon_native'])
            oriented,aligned=transform_xy(g[['x_native','y_native']].to_numpy(),center,r['orientation'])
            center_oriented,_=transform_xy(center,center,r['orientation'])
            g[['x_oriented','y_oriented']]=oriented;g[['x_aligned','y_aligned']]=aligned
            g['x_ventricle'],g['y_ventricle']=center_oriented
        g.to_csv(directory/f'job_{m.job}_objects.csv.gz',index=False,compression='gzip')


def write_protocol(root,out,manifest,reviews,parameters=None):
    """Content-addressed snapshots: reruns never overwrite a previous protocol."""
    root=Path(root);out=Path(out);p=dict(parameters or PARAMETERS)
    versions={name:importlib.metadata.version(name) for name in ('numpy','pandas','scipy','matplotlib','plotly','nbformat','nbclient')}
    codepaths=['bio_visuals/handoff.py','bio_visuals/handoff_plots.py','bio_visuals/explorer.py','bio_visuals/explorer.html','bio_visuals/standalone_app.py','bio_visuals/condition_comparison.py','bio_visuals/condition_views.py','bio_visuals/condition_comparison.html','bio_visuals/condition_layout.py','bio_visuals/entire_dataset.py','bio_visuals/data.py','scripts/build_spatial_handoff_notebook.py','scripts/build_scientific_notebook.py','scripts/serve_spatial_app.py']
    code={name:_sha256(root/name) for name in codepaths}
    handoff=root/'HALO_spatial_analysis_raw_data_handoff.pdf'
    document=dict(parameters=p,section_reviews=reviews,
                  input_files=manifest[['job','source_file','sha256','is_primary','rows']].to_dict('records'),
                  code_sha256=code,software=versions,python=platform.python_version(),
                  handoff_sha256=_sha256(handoff) if handoff.exists() else 'not supplied on this rerun')
    digest=hashlib.sha256(json.dumps(document,sort_keys=True,ensure_ascii=False).encode()).hexdigest()[:12]
    directory=out/'protocols';directory.mkdir(parents=True,exist_ok=True)
    path=directory/f'v{VERSION}_{digest}.json'
    if not path.exists():
        document['execution_utc']=datetime.now(timezone.utc).isoformat()
        document['git_commit']=subprocess.run(['git','rev-parse','HEAD'],cwd=root,capture_output=True,text=True).stdout.strip()
        document['git_status']=subprocess.run(['git','status','--short'],cwd=root,capture_output=True,text=True).stdout.splitlines()
        document['investigator']='Codex computational first pass; anatomical review pending'
        path.write_text(json.dumps(document,indent=2,ensure_ascii=False)+'\n')
    (out/'latest_protocol.txt').write_text(str(path.relative_to(out))+'\n')
    return path

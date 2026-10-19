"""Notebook-first scientific figures for the HALO raw-data handoff."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, ListedColormap
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D
from .handoff import ANIMAL_ORDER, ARMS, SIZE_LABELS, transform_xy

INK='#173447';TEAL='#087E8B';ORANGE='#D06738';PALE='#edf3f5'
GROUP_COLORS={'Home cage':TEAL,'Running':ORANGE}
ANIMAL_COLORS=dict(zip(ANIMAL_ORDER,['#0072B2','#009E73','#7C52A1','#D55E00','#E69F00','#CC79A7']))


def style():
    plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','axes.edgecolor':'#bccbd1',
                         'axes.labelcolor':INK,'text.color':INK,'xtick.color':INK,'ytick.color':INK,
                         'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,
                         'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False,
                         'savefig.facecolor':'white','svg.fonttype':'none','figure.dpi':110})


def save(fig,out,name):
    target=Path(out)/'figures';target.mkdir(parents=True,exist_ok=True)
    fig.savefig(target/f'{name}.png',dpi=300,bbox_inches='tight')
    fig.savefig(target/f'{name}.svg',bbox_inches='tight')
    return target/f'{name}.png'


def extent(g):
    return [g['x_edges'][0],g['x_edges'][-1],g['y_edges'][-1],g['y_edges'][0]]


def image(ax,g,z,cmap='viridis',norm=None):
    return ax.imshow(np.ma.masked_invalid(z),extent=extent(g),origin='upper',
                     cmap=cmap,norm=norm,interpolation='nearest',aspect='equal')


def coordinates(ax):
    ax.set_xlabel('x · export units');ax.set_ylabel('y · export units')
    ax.ticklabel_format(style='sci',scilimits=(4,4),axis='both');ax.tick_params(labelsize=8)


def title(fig,heading,sub):
    height=fig.get_size_inches()[1]
    fig.text(.03,1+.45/height,heading,fontsize=21,fontweight='bold',ha='left',va='bottom',color=INK)
    fig.text(.03,1+.12/height,sub,fontsize=10,ha='left',va='bottom')


def workflow_figure():
    fig,ax=plt.subplots(figsize=(16,3.2));ax.axis('off')
    ax.set_position([.03,.02,.94,.95])
    boxes=[('01  Experiment','6 animals · 3 per condition\nTwo physical PLA sections\nChAT is the available marker'),
           ('02  Raw measurements','13 CSV exports audited\nInterrupted run excluded\n12 sections · 3.01M objects'),
           ('03  This notebook','Native overlays + heat maps\nAnimal-level summaries\nInteractive signal-height views'),
           ('04  Scientist review','Trace cavity · confirm flips\nVerify scale + tissue support\nReview candidate QC package'),
           ('05  Registered maps','After review only:\nEqual-animal condition maps\nGeometry-aware null tests')]
    for i,(h,b) in enumerate(boxes):
        x=.012+i*.199
        ax.add_patch(FancyBboxPatch((x,.12),.174,.77,boxstyle='round,pad=.007,rounding_size=.02',
                                   transform=ax.transAxes,facecolor=PALE if i!=3 else '#fff0da',edgecolor='#cbd9df'))
        ax.text(x+.009,.76,h,transform=ax.transAxes,weight='bold',fontsize=11)
        ax.text(x+.009,.58,b,transform=ax.transAxes,va='top',fontsize=9.5,linespacing=1.7)
        if i<4:ax.annotate('',xy=(x+.193,.5),xytext=(x+.182,.5),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',color=INK))
    title(fig,'From the experiment to a spatial result','All values remain traceable to raw files; the PLA arms stay separate at every stage.')
    return fig


def native_atlas(fine,metadata,arm):
    rows=metadata[metadata.pla_arm==arm].set_index('animal_id').loc[ANIMAL_ORDER]
    fig,axes=plt.subplots(2,3,figsize=(15,11),layout='constrained')
    maxp=max(float(fine[str(r.job)]['n_pla'].max()) for r in rows.itertuples())
    maxc=max(float(fine[str(r.job)]['n_chat'].max()) for r in rows.itertuples())
    pn=LogNorm(1,max(2,maxp));cn=LogNorm(1,max(2,maxc))
    for ax,(animal,r) in zip(axes.flat,rows.iterrows()):
        g=fine[str(r.job)]
        image(ax,g,np.where(g['n_chat']>0,g['n_chat'],np.nan),'Greys',cn)
        im=image(ax,g,np.where(g['n_pla']>0,g['n_pla'],np.nan),'magma',pn)
        # Color is PLA count; gray is independently detected ChAT, not raw fluorescence.
        ax.set_title(f'{animal} · {r.running_group} · job {r.job}',color=GROUP_COLORS[r.running_group],loc='left')
        coordinates(ax)
    fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.7,pad=.015,label='PLA objects / 100 × 100-unit bin · log color')
    title(fig,f'{arm} | tissue reconstructed from detections','Gray: ChAT midpoint counts. Color: PLA midpoint counts. Shared scales within this arm; no microscopy pixels or anatomical alignment.')
    return fig


def fraction_atlas(grids,metadata,arm):
    rows=metadata[metadata.pla_arm==arm].set_index('animal_id').loc[ANIMAL_ORDER]
    fig,axes=plt.subplots(2,3,figsize=(15,10.5),layout='constrained')
    for ax,(animal,r) in zip(axes.flat,rows.iterrows()):
        g=grids[str(r.job)]
        image(ax,g,np.where(g['support'],1.,np.nan),ListedColormap(['#e7ecef']))
        im=image(ax,g,g['fraction'],'viridis',Normalize(0,100))
        supported=int(np.isfinite(g['fraction']).sum())
        ax.set_title(f'{animal} · {r.running_group} | {supported} supported bins',color=GROUP_COLORS[r.running_group],loc='left')
        coordinates(ax)
    fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.7,pad=.015,label='PLA objects with ChAT present (%)')
    title(fig,f'{arm} | where detected PLA is ChAT-positive','1,000 × 1,000 export-unit bins; ≥10 PLA objects. Pale gray: objects observed but insufficient PLA. White: no object support. No smoothing.')
    return fig


def detailed_section(fine,grid,meta):
    fig,axes=plt.subplots(2,3,figsize=(16,10.5),layout='constrained')
    panels=[('PLA counts','n_pla',fine,'magma','Objects / 100-unit bin',True),
            ('ChAT counts · independent context','n_chat',fine,'Greys','Objects / 100-unit bin',True),
            ('PLA–ChAT association','fraction',grid,'viridis','ChAT-positive PLA (%)',False),
            ('PLA area with ChAT present','positive_pla_area',grid,'magma','Sum of whole positive-PLA areas (µm² / bin)',True),
            ('Mean PLA size','mean_pla_area',grid,'viridis','Mean object area (µm²)',True),
            ('Large PLA · fixed ≥1 µm²','large_fraction',grid,'cividis','PLA ≥1 µm² (%)',False)]
    for ax,(name,key,g,cmap,label,log) in zip(axes.flat,panels):
        z=g[key].copy().astype(float)
        if key in ('positive_pla_area',):z[g['n_pla']<g['min_pla']]=np.nan
        if log:
            positive=z[np.isfinite(z)&(z>0)]
            norm=LogNorm(max(positive.min(),1e-6),max(positive.max(),positive.min()*1.01))
            z[z<=0]=np.nan
        else:norm=Normalize(0,100)
        im=image(ax,g,z,cmap,norm);ax.set_title(name,loc='left');coordinates(ax)
        fig.colorbar(im,ax=ax,shrink=.75,pad=.02,label=label)
    title(fig,f'{meta.animal_id} · {meta.pla_arm} | complementary spatial readouts',
          f'{meta.running_group} · job {meta.job}. Native coordinates. Area sums are not union tissue area; metric bins are 1,000 units, unsmoothed.')
    return fig


def group_summary_figure(summary):
    fig,axes=plt.subplots(2,3,figsize=(15,8.5),layout='constrained')
    metrics=[('pla_association_pct','ChAT-positive PLA (%)'),
             ('marker_positive_area_share_pct','PLA area in ChAT-positive objects (%)'),
             ('chat_objects','ChAT detections (objects)')]
    for i,arm in enumerate(ARMS.values()):
        d=summary[summary.pla_arm==arm]
        for ax,(key,label) in zip(axes[i],metrics):
            for x,group in enumerate(('Home cage','Running')):
                g=d[d.running_group==group].sort_values('animal_id')
                offsets=np.linspace(-.13,.13,len(g))
                ax.scatter(x+offsets,g[key],s=85,c=GROUP_COLORS[group],edgecolor='white',linewidth=1,zorder=3)
                for dx,r in zip(offsets,g.itertuples()):ax.annotate(r.animal_id,(x+dx,getattr(r,key)),xytext=(5,4),textcoords='offset points',fontsize=9)
                ax.plot([x-.22,x+.22],[g[key].mean()]*2,color=INK,lw=2)
            ax.set_xticks([0,1],['Home cage','Running']);ax.set_xlim(-.5,1.55)
            ax.set_title(f'{arm}\n{label}',loc='left');ax.grid(axis='y',color='#e8eef1',zorder=0)
            if key!='chat_objects':ax.set_ylim(0,100)
            else:ax.ticklabel_format(style='sci',scilimits=(3,3),axis='y')
    title(fig,'The animal is the replicate','Each dot is one animal/section; bars are equal-animal means. n=3 per condition in each arm. Counts reflect sampling and detection, not neuron number.')
    return fig


def size_figure(size_table,sensitivity):
    fig,axes=plt.subplots(2,2,figsize=(14,9),layout='constrained')
    for i,arm in enumerate(ARMS.values()):
        for ax,data,xcol,xlabels in [(axes[i,0],size_table,'size_bin',SIZE_LABELS),
                                    (axes[i,1],sensitivity,'min_area_um2',[0.,.03,.1,.3,1.])]:
            d=data[data.pla_arm==arm]
            for animal,g in d.groupby('animal_id',observed=True):
                g=g.set_index(xcol).reindex(xlabels)
                group=g.running_group.iloc[0]
                ax.plot(range(len(xlabels)),g.fraction_pct,lw=1,alpha=.5,color=GROUP_COLORS[group])
            for group in ('Home cage','Running'):
                y=d[d.running_group==group].groupby(xcol,observed=True).fraction_pct.mean().reindex(xlabels)
                ax.plot(range(len(xlabels)),y,'o-',lw=2.8,ms=5,color=GROUP_COLORS[group],label=f'{group} · equal-animal mean')
            ax.set_xticks(range(len(xlabels)),[str(x) for x in xlabels]);ax.set_ylim(0,100)
            ax.set_ylabel('ChAT-positive PLA (%)');ax.grid(axis='y',color='#e8eef1')
            ax.set_xlabel('Fixed PLA area bin (µm²)' if xcol=='size_bin' else 'Include PLA with area ≥ threshold (µm²)')
            ax.set_title(arm+' | '+('size-specific association' if xcol=='size_bin' else 'filter sensitivity'),loc='left')
    axes[0,0].legend(frameon=False,loc='lower right',fontsize=9)
    title(fig,'Does punctum size change the comparison?','Thin lines: individual animals. Thick lines: group means. Fixed cutoffs are descriptive choices, not validated biological classes; all denominators are exported.')
    return fig


def qc_section(grid,review):
    fig,axes=plt.subplots(1,3,figsize=(15,5.4),layout='constrained')
    n=grid['n_all'];norm=LogNorm(1,max(2,n.max()))
    for ax in axes[:2]:image(ax,grid,np.where(n>0,n,np.nan),'magma',norm);coordinates(ax)
    axes[0].set_title('All detections · geometry only',loc='left')
    axes[1].set_title('Unreviewed enclosed-void proposal',loc='left')
    if review.get('candidate'):
        v=np.array(review['candidate']['polygon_native']);c=review['candidate']['centroid_native']
        axes[1].plot(v[:,0],v[:,1],color='#21ebd2',lw=1.8)
        axes[1].plot(c[0],c[1],'+',ms=12,mew=2,color='#21ebd2')
    else:
        axes[1].text(.5,.08,'No qualifying enclosed void\nManual tracing needed',transform=axes[1].transAxes,ha='center',color=INK,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
    image(axes[2],grid,np.where(grid['support'],1.,np.nan),ListedColormap(['#467888']))
    axes[2].set_title('Occupied-bin support · provisional',loc='left');coordinates(axes[2])
    title(fig,f'{review["animal_id"]} · {review["pla_arm"]} | registration review',
          f'Job {review["job"]}. 250-unit geometry grid. A detected void is a proposal, not a validated ventricle. Missing bins stay unknown.')
    return fig


def flip_figure(grid,review):
    fig,axes=plt.subplots(1,4,figsize=(16,5),layout='constrained')
    n=grid['n_all'];xc=(grid['x_edges'][1:]+grid['x_edges'][:-1])/2;yc=(grid['y_edges'][1:]+grid['y_edges'][:-1])/2
    xx,yy=np.meshgrid(xc,yc);mask=n>0
    center=review['candidate']['centroid_native'] if review.get('candidate') else [0,0]
    for ax,flip in zip(axes,('none','horizontal','vertical','both')):
        _,xy=transform_xy(np.c_[xx[mask],yy[mask]],center,flip)
        ax.scatter(xy[:,0],xy[:,1],c=np.log1p(n[mask]),cmap='Greys',s=2,rasterized=True)
        if review.get('candidate'):
            _,v=transform_xy(np.array(review['candidate']['polygon_native']),center,flip)
            ax.plot(v[:,0],v[:,1],color=TEAL);ax.plot(0,0,'+',color=ORANGE,ms=10)
        ax.set_title(flip+' · preview');ax.invert_yaxis();ax.set_aspect('equal');coordinates(ax)
    title(fig,f'{review["animal_id"]} · {review["pla_arm"]} | test the four allowed orientations',
          'Every panel is a geometry-only possibility. No choice is marked approved; no rotation, scaling or warping is applied.')
    return fig


def silhouette_overlay(geometry,reviews,arm):
    fig,ax=plt.subplots(figsize=(10,8),layout='constrained');missing=[]
    for r in [r for r in reviews if r['pla_arm']==arm]:
        g=geometry[r['job']]
        if not r.get('candidate'):missing.append(r['animal_id']);continue
        from scipy import ndimage as ndi
        # Closing only makes a readable contour; no group metric uses the filled interior.
        mask=ndi.binary_closing(g['support'],iterations=2)
        xc=(g['x_edges'][1:]+g['x_edges'][:-1])/2;yc=(g['y_edges'][1:]+g['y_edges'][:-1])/2
        xx,yy=np.meshgrid(xc,yc)
        _,xy=transform_xy(np.c_[xx.ravel(),yy.ravel()],r['candidate']['centroid_native'],r['orientation'])
        ax.contour(xy[:,0].reshape(xx.shape),xy[:,1].reshape(yy.shape),mask,levels=[.5],colors=[ANIMAL_COLORS[r['animal_id']]],linewidths=1.2,alpha=.8)
    ax.plot(0,0,'+',c=INK,ms=15,mew=2);ax.invert_yaxis();ax.set_aspect('equal');coordinates(ax)
    handles=[Line2D([0],[0],color=ANIMAL_COLORS[r['animal_id']],label=r['animal_id']+' · '+r['running_group']) for r in reviews if r['pla_arm']==arm and r.get('candidate')]
    ax.legend(handles=handles,frameon=False,loc='upper left',bbox_to_anchor=(1,1))
    title(fig,f'{arm} | do the proposed centers agree?',
          'UNREVIEWED QC overlay · candidate translation only · not an anatomical comparison')
    fig.text(.03,-.025,'No candidate, omitted: '+(', '.join(missing) if missing else 'none')+'. No scaling to improve fit. Review residual shape differences before biological group maps.',fontsize=9)
    return fig


def coverage_figure(coverage):
    fig,axes=plt.subplots(2,2,figsize=(12,10),layout='constrained')
    g={'x_edges':coverage['x_edges'],'y_edges':coverage['y_edges']}
    cmap=ListedColormap(['#edf1f3','#c4d8de','#6799aa','#154d67'])
    for i,arm in enumerate(ARMS.values()):
        for j,group in enumerate(('Home cage','Running')):
            ax=axes[i,j];z=coverage['counts'][(arm,group)]
            im=image(ax,g,z,cmap,Normalize(-.5,3.5))
            n=sum(r['pla_arm']==arm and r['running_group']==group for r in coverage['used'])
            ax.set_title(f'{arm} · {group}\n{n}/3 sections have a geometric candidate',loc='left');coordinates(ax)
    fig.colorbar(im,ax=axes.ravel().tolist(),ticks=[0,1,2,3],shrink=.65,label='Candidate sections with observed object support')
    title(fig,'Coverage is part of the result','UNREVIEWED QC only. Candidate-centered occupancy is not validated tissue coverage. No condition-effect map is authorized by this preview.')
    return fig


def proximity_figure(curves):
    fig,axes=plt.subplots(1,2,figsize=(14,5.5),layout='constrained')
    for ax,arm in zip(axes,ARMS.values()):
        data=curves[curves.pla_arm==arm]
        for (animal,sel),g in data.groupby(['animal_id','selection'],observed=True):
            ax.plot(g.distance,g.probability*100,color=ANIMAL_COLORS[animal],lw=1.7 if sel=='All PLA' else 1,
                    ls='-' if sel=='All PLA' else '--',label=animal if sel=='All PLA' else '_nolegend_')
        xmax=data[data.probability<=.99].distance.max()
        ax.set_xlim(0,xmax);ax.set_ylim(0,100);ax.grid(color='#edf1f3');ax.set_title(arm,loc='left')
        ax.set_xlabel('Distance to nearest ChAT midpoint · export units');ax.set_ylabel('PLA objects within distance (%)')
        ax.legend(ncol=3,frameon=False,fontsize=9,loc='lower right')
    title(fig,'How close are PLA detections to ChAT detections?',
          'Solid: all PLA. Dashed: ≥500-unit conservative support buffer + nearest neighbor inside observed support. Uncalibrated; no enrichment or binding claim.')
    return fig


def zoom_boxes(objects,queue,job):
    r=queue[(queue.job==job)&(queue.reason=='largest PLA count')].iloc[0]
    frame=objects[(objects.job==job)&objects.x.between(r.xmin,r.xmax)&objects.y.between(r.ymin,r.ymax)]
    fig,axes=plt.subplots(1,2,figsize=(13,6),layout='constrained')
    for typ,color,size in [('ChAT','#849da8',3),('PLA',ORANGE,10)]:
        g=frame[frame.object_type==typ]
        axes[0].scatter(g.x,g.y,s=size,color=color,alpha=.7,label=f'{typ}: {len(g):,}',rasterized=True)
        for p in g.itertuples():
            axes[1].add_patch(Rectangle((p.xmin,p.ymin),p.xmax-p.xmin,p.ymax-p.ymin,
                                       edgecolor=color,facecolor='none',lw=.7,alpha=.6))
    for ax in axes:
        ax.set_xlim(r.xmin,r.xmax);ax.set_ylim(r.ymax,r.ymin);ax.set_aspect('equal');coordinates(ax)
    axes[0].legend(frameon=False);axes[0].set_title('Object midpoints',loc='left');axes[1].set_title('Original bounding boxes',loc='left')
    title(fig,f'{r.animal_id} · {r.pla_arm} | inspect a populated window',
          f'Job {job}. Highest PLA-count 1,000-unit bin, chosen by a declared rule. Boxes describe object extent, not segmentation boundaries or microscopy.')
    return fig


def interactive_native(fine,grids,metadata,arm):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    rows=metadata[metadata.pla_arm==arm].set_index('animal_id').loc[ANIMAL_ORDER]
    fig=make_subplots(rows=1,cols=3,subplot_titles=('Reconstructed overlay','ChAT-positive PLA (%)','PLA denominator / 1,000-unit bin'),horizontal_spacing=.08)
    maxp=max(np.log1p(fine[str(r.job)]['n_pla']).max() for r in rows.itertuples())
    maxc=max(np.log1p(fine[str(r.job)]['n_chat']).max() for r in rows.itertuples())
    maxn=max(grids[str(r.job)]['n_pla'].max() for r in rows.itertuples())
    for i,(animal,r) in enumerate(rows.iterrows()):
        f=fine[str(r.job)];g=grids[str(r.job)]
        for z,source,col,colors,zmin,zmax,label,show in [
            (np.where(f['n_chat']>0,np.log1p(f['n_chat']),np.nan),f,1,'Greys',0,maxc,'log(1 + ChAT count)',False),
            (np.where(f['n_pla']>0,np.log1p(f['n_pla']),np.nan),f,1,'Magma',0,maxp,'log(1 + PLA count)',False),
            (g['fraction'],g,2,'Viridis',0,100,'ChAT-positive PLA (%)',True),
            (np.where(g['support'],g['n_pla'],np.nan),g,3,'Magma',0,maxn,'PLA count',True)]:
            x=(source['x_edges'][1:]+source['x_edges'][:-1])/2;y=(source['y_edges'][1:]+source['y_edges'][:-1])/2
            fig.add_trace(go.Heatmap(x=x,y=y,z=z,colorscale=colors,zmin=zmin,zmax=zmax,
                                    visible=i==0,showscale=show,hoverongaps=False,
                                    colorbar=dict(title=dict(text=label,side='right'),len=.7,x=.63 if col==2 else 1.01,thickness=12),
                                    hovertemplate='x=%{x:,.0f}<br>y=%{y:,.0f}<br>'+label+'=%{z:.2f}<extra>'+animal+'</extra>'),row=1,col=col)
    steps=[]
    for i,(animal,r) in enumerate(rows.iterrows()):
        visible=[i*4<=k<(i+1)*4 for k in range(len(fig.data))]
        steps.append(dict(method='update',args=[{'visible':visible},{'title.text':f'{arm} · {animal} · {r.running_group} | native coordinate explorer'}],label=animal))
    fig.update_layout(template='plotly_white',height=620,width=1200,font=dict(family='Arial',color=INK),
                      title=dict(text=f'{arm} · {ANIMAL_ORDER[0]} · Home cage | native coordinate explorer',font=dict(size=20)),
                      sliders=[dict(steps=steps,currentvalue=dict(prefix='Animal: '),pad=dict(t=50))],margin=dict(t=95,b=110,l=60,r=100))
    for col in range(1,4):
        suffix='' if col==1 else str(col)
        fig.update_xaxes(title_text='x · export units',row=1,col=col)
        fig.update_yaxes(title_text='y · export units' if col==1 else None,
                         autorange='reversed',scaleanchor='x'+suffix,scaleratio=1,row=1,col=col)
    return fig


def interactive_surface(surface_grids,metadata,arm):
    import plotly.graph_objects as go
    rows=metadata[metadata.pla_arm==arm].set_index('animal_id').loc[ANIMAL_ORDER]
    maxn=max(surface_grids[str(r.job)]['n_pla'].max() for r in rows.itertuples())
    fig=go.Figure()
    for i,(animal,r) in enumerate(rows.iterrows()):
        g=surface_grids[str(r.job)]
        x=(g['x_edges'][1:]+g['x_edges'][:-1])/2;y=(g['y_edges'][1:]+g['y_edges'][:-1])/2
        z=np.where(g['support'],g['n_pla'],np.nan)
        custom=np.stack([g['n_chat'],g['fraction'],g['n_pla']],axis=-1)
        fraction_text=np.array([f'{v:.1f}%' if np.isfinite(v) else 'not shown (PLA <10)'
                                for v in g['fraction'].ravel()]).reshape(g['fraction'].shape)
        fig.add_trace(go.Surface(x=x,y=y,z=np.log1p(z),customdata=custom,cmin=0,cmax=np.log1p(maxn),colorscale='Viridis',
                                text=fraction_text,
                                visible=i==0,connectgaps=False,name=animal,
                                colorbar=dict(title='PLA / bin · log height',thickness=18,tickvals=np.log1p([0,5,20,100,500,2500]),ticktext=['0','5','20','100','500','2500']),
                                hovertemplate='x=%{x:,.0f}<br>y=%{y:,.0f}<br>PLA=%{customdata[2]:.0f}<br>ChAT=%{customdata[0]:.0f}<br>ChAT-positive PLA=%{text}<extra>'+animal+'</extra>',
                                lighting=dict(ambient=.8,diffuse=.6,specular=.15),
                                contours=dict(z=dict(show=True,usecolormap=True,highlightcolor='#0eafc5',project_z=True))))
    steps=[]
    for i,(animal,r) in enumerate(rows.iterrows()):
        steps.append(dict(method='update',label=animal,args=[{'visible':[j==i for j in range(len(rows))]},
                     {'title.text':f'{arm} · {animal} · {r.running_group}<br><sup>Height = log(1 + PLA count / 500-unit bin) · not anatomical depth</sup>'}]))
    fig.update_layout(template='plotly_white',height=740,width=1100,font=dict(family='Arial',color=INK),
                      title=dict(text=f'{arm} · {ANIMAL_ORDER[0]} · Home cage<br><sup>Height = log(1 + PLA count / 500-unit bin) · not anatomical depth</sup>'),
                      scene=dict(xaxis_title='x · export units',yaxis=dict(title='y · export units',autorange='reversed'),
                                 zaxis=dict(title='PLA / bin · log height',range=[0,np.log1p(maxn)],tickvals=np.log1p([0,5,20,100,500,2500]),ticktext=['0','5','20','100','500','2500']),
                                 aspectmode='manual',aspectratio=dict(x=1,y=1,z=.55),
                                 camera=dict(eye=dict(x=1.3,y=1.4,z=1.1))),
                      sliders=[dict(steps=steps,currentvalue=dict(prefix='Section: '),pad=dict(t=25))],
                      updatemenus=[dict(type='buttons',direction='right',x=.02,y=1.05,buttons=[
                          dict(label='Oblique view',method='relayout',args=[{'scene.camera':dict(eye=dict(x=1.3,y=1.4,z=1.1))}]),
                          dict(label='Look down',method='relayout',args=[{'scene.camera':dict(eye=dict(x=0,y=0,z=2.5),up=dict(x=0,y=1,z=0))}])])],
                      margin=dict(l=0,r=70,t=130,b=90))
    # Preserve x:y geometry within each native section, even when extents differ.
    for step,(animal,r) in zip(steps,rows.iterrows()):
        g=surface_grids[str(r.job)];xr=np.ptp(g['x_edges']);yr=np.ptp(g['y_edges'])
        step['args'][1]['scene.aspectratio']=dict(x=1,y=yr/xr,z=.55)
    # Plotly copies dictionaries on assignment; refresh after adding each ratio.
    fig.layout.sliders[0].steps=steps
    first=surface_grids[str(rows.iloc[0].job)]
    fig.update_layout(scene_aspectratio=dict(x=1,y=np.ptp(first['y_edges'])/np.ptp(first['x_edges']),z=.55))
    return fig


def surface_preview(grids,metadata,arm):
    rows=metadata[metadata.pla_arm==arm].set_index('animal_id').loc[ANIMAL_ORDER]
    # Fixed spacing is more reliable than constrained_layout for six 3D axes
    # when the same figure is exported at both notebook and 300-DPI sizes.
    fig=plt.figure(figsize=(15,9))
    fig.subplots_adjust(left=.03,right=.97,bottom=.04,top=.92,wspace=.10,hspace=.12)
    vmax=max(grids[str(r.job)]['n_pla'].max() for r in rows.itertuples())
    for i,(animal,r) in enumerate(rows.iterrows()):
        g=grids[str(r.job)];ax=fig.add_subplot(2,3,i+1,projection='3d')
        x=(g['x_edges'][1:]+g['x_edges'][:-1])/2;y=(g['y_edges'][1:]+g['y_edges'][:-1])/2
        xx,yy=np.meshgrid(x,y);z=np.where(g['support'],g['n_pla'],np.nan)
        ax.plot_surface(xx,yy,np.log1p(z),cmap='viridis',vmin=0,vmax=np.log1p(vmax),linewidth=0,antialiased=False,rcount=z.shape[0],ccount=z.shape[1])
        ax.set_zlim(0,np.log1p(vmax));ax.set_zticks(np.log1p([0,5,20,100,500,2500]),['0','5','20','100','500','2500']);ax.invert_yaxis();ax.view_init(40,-60)
        ax.set_box_aspect((1,np.ptp(y)/np.ptp(x),.55));ax.tick_params(labelsize=6)
        ax.set_title(f'{animal} · {r.running_group}',fontsize=11,color=GROUP_COLORS[r.running_group])
        ax.set_xlabel('x',fontsize=8);ax.set_ylabel('y',fontsize=8);ax.set_zlabel('PLA / bin · log',fontsize=8)
    title(fig,f'{arm} | spatial signal in three dimensions','Static companion to the interactive view. x/y are native export units; height is log(1 + count / 500-unit bin), not tissue depth. Shared height scale.')
    return fig


def landmark_atlas(geometry,reviews,arm):
    rows={r['animal_id']:r for r in reviews if r['pla_arm']==arm}
    fig,axes=plt.subplots(2,3,figsize=(15,11),layout='constrained')
    for ax,animal in zip(axes.flat,ANIMAL_ORDER):
        r=rows[animal];g=geometry[r['job']];z=g['n_all']
        image(ax,g,np.where(z>0,z,np.nan),'Greys',LogNorm(1,max(2,z.max())))
        candidate=r.get('candidate')
        if candidate:
            v=np.array(candidate['polygon_native']);c=candidate['centroid_native']
            ax.plot(v[:,0],v[:,1],color=ORANGE,lw=2)
            ax.plot(*c,'+',color=ORANGE,ms=12,mew=2)
            note='unreviewed void proposal'
        else:note='no enclosed-void proposal'
        ax.set_title(f'{animal} · {note}',fontsize=11,loc='left');coordinates(ax)
    title(fig,f'{arm} | inspect the cavity before registration',
          'Orange polygon: automated trace of an enclosed geometric void. It may be only part of the ventricular cavity or a different gap; review every proposal.')
    return fig


def registration_widget(geometry,records,path):
    """Live-kernel review aid; saves only on an explicit button click."""
    import copy,json
    import ipywidgets as W
    from IPython.display import display,clear_output
    state=copy.deepcopy(records)
    choices=[(r['animal_id']+' · '+r['pla_arm']+' · '+r['job'],i) for i,r in enumerate(state)]
    select=W.Dropdown(options=choices,description='Section:',layout=W.Layout(width='550px'))
    polygon=W.Textarea(description='Polygon:',layout=W.Layout(width='95%',height='100px'))
    flip=W.Dropdown(options=['none','horizontal','vertical','both'],description='Flip:')
    landmark=W.Checkbox(description='Cavity trace reviewed');orientation=W.Checkbox(description='Orientation reviewed')
    scale=W.Checkbox(description='Common coordinate scale verified');support=W.Checkbox(description='Support proxy reviewed')
    support_qc=W.Dropdown(options=['uncertain','complete','minor loss','substantial loss'],description='Support:')
    reviewer=W.Text(description='Reviewer:');second=W.Text(description='2nd review:')
    notes=W.Textarea(description='Notes:',layout=W.Layout(width='95%',height='70px'))
    button=W.Button(description='Save section review',button_style='info');preview=W.Button(description='Preview draft polygon')
    output=W.Output();message=W.Output()
    def fill(change=None):
        r=state[select.value]
        v=r.get('reviewed_polygon_native') or (r.get('candidate') or {}).get('polygon_native')
        polygon.value=json.dumps(v)
        flip.value=r['orientation'];landmark.value=r['landmark_reviewed'];orientation.value=r['orientation_reviewed']
        scale.value=r['scale_common_confirmed'];support.value=r['support_reviewed'];support_qc.value=r['support_qc']
        reviewer.value=r['reviewer'];second.value=r['second_reviewer'];notes.value=r['review_notes']
        show()
    def draft():
        from .handoff import polygon_centroid
        r=copy.deepcopy(state[select.value]);v=json.loads(polygon.value)
        if v is not None:polygon_centroid(v)
        r.update(reviewed_polygon_native=v,orientation=flip.value,landmark_reviewed=landmark.value,
                 orientation_reviewed=orientation.value,scale_common_confirmed=scale.value,
                 support_reviewed=support.value,support_qc=support_qc.value,reviewer=reviewer.value,
                 second_reviewer=second.value,review_notes=notes.value)
        return r
    def show(_=None):
        with output:
            clear_output(wait=True)
            try:
                r=draft()
                if r['reviewed_polygon_native']:
                    from .handoff import polygon_centroid
                    r['candidate']=dict(polygon_native=r['reviewed_polygon_native'],centroid_native=polygon_centroid(r['reviewed_polygon_native']).tolist())
                fig=qc_section(geometry[r['job']],r);display(fig);plt.close(fig)
            except Exception as error:print('Draft not displayed:',error)
    def persist(_):
        with message:
            clear_output(wait=True)
            try:
                from .handoff import review_status
                r=draft();review_status([r]);state[select.value]=r
                Path(path).write_text(json.dumps(state,indent=2)+'\n')
                print('Saved section record. Rerun the review-status and gated-map cells to use it. No approval flags are set automatically.')
            except Exception as error:print('Not saved:',error)
    select.observe(fill,names='value');button.on_click(persist);preview.on_click(show);fill()
    return W.VBox([W.HTML('<b>Scientist review aid.</b> Edit native x,y polygon vertices, inspect geometry, and record actual review. A proposal is not evidence of approval.'),
                   select,polygon,W.HBox([flip,support_qc]),W.HBox([landmark,orientation]),W.HBox([scale,support]),
                   W.HBox([reviewer,second]),notes,W.HBox([preview,button]),message,output])

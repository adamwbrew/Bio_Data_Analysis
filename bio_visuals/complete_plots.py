"""High-resolution native-coordinate views and explicit export-level summaries."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm, ListedColormap
from matplotlib.patches import Rectangle
from matplotlib.ticker import PercentFormatter, MaxNLocator

TYPE_COLORS = {'Pfkp':'#0072B2', 'Hk1':'#D55E00'}


def _finish(fig, title, caption, bottom=.15):
    fig.suptitle(title, x=.025, ha='left', fontsize=17, fontweight='bold')
    fig.text(.025,.025,caption,fontsize=9,ha='left',va='bottom',color='#465365')
    fig.subplots_adjust(left=.07,right=.96,top=.86,bottom=bottom,wspace=.30)
    return fig


def overview_figure(summary):
    """Every export shown; no blank run-by-group matrix or inferred treatment."""
    data = summary.sort_values(['condition','pla_type']).reset_index(drop=True)
    fig, axes = plt.subplots(1,3,figsize=(16,7))
    y = np.arange(len(data))
    labels = data.condition.astype(str)+' · '+data.pla_type.astype(str)+'  ['+data.job.astype(str)+']'
    for ax, metric, title in zip(axes, ['pla_objects','pla_association_pct','area_weighted_overlap_pct'],
                                ['PLA detections','PLA objects with ChAT present','Area-weighted PLA overlap']):
        for i,row in data.iterrows():
            color=TYPE_COLORS[str(row.pla_type)]
            ax.plot([0,row[metric]],[i,i],color=color,alpha=.3)
            ax.scatter(row[metric],i,c=color,s=40,zorder=3)
            ax.annotate(f'{row[metric]:,.0f}' if metric=='pla_objects' else f'{row[metric]:.2f}%',
                        (row[metric],i),xytext=(5,0),textcoords='offset points',va='center',fontsize=8)
        ax.set_yticks(y)
        ax.set_yticklabels(labels if ax is axes[0] else ['']*len(data))
        ax.invert_yaxis();ax.grid(axis='x',alpha=.18);ax.set_title(title)
        if metric!='pla_objects':
            ax.set_xlim(0,108);ax.set_xlabel('Percentage (%)')
        else:ax.set_xlim(0,data[metric].max()*1.25);ax.set_xlabel('Exported objects')
    _finish(fig,'Twelve primary analyses: counts and overlap',
            'Blue = Pfkp–TDP43; orange = Hk1–TDP43. Labels identify exports; they do not assign running groups.\nCounts lack a tissue-area denominator; percentages describe segmented objects, not protein binding.')
    fig.subplots_adjust(left=.17)
    return fig


def paired_type_figure(summary):
    """Link literal labels across assays, without implying same-cell matching."""
    fig,axes=plt.subplots(1,3,figsize=(16,6.5))
    colors=plt.cm.tab10(np.linspace(0,.8,summary.condition.nunique()))
    for ax,metric,title in zip(axes,['pla_association_pct','area_weighted_overlap_pct','pla_objects'],
                              ['Object association','Area-weighted overlap','PLA detections']):
        for (label,g),color in zip(summary.groupby('condition',observed=True,sort=True),colors):
            if set(g.pla_type)!=set(TYPE_COLORS):raise ValueError('Each label needs both PLA types')
            values=g.set_index('pla_type').loc[['Hk1','Pfkp'],metric].to_numpy()
            ax.plot([0,1],values,'o-',label=str(label),color=color,lw=1.4,markersize=5)
        ax.set_xticks([0,1],['Hk1–TDP43','Pfkp–TDP43']);ax.set_xlim(-.15,1.15)
        ax.set_title(title);ax.grid(axis='y',alpha=.2)
        if metric!='pla_objects':ax.set_ylim(0,100);ax.set_ylabel('Percentage (%)')
        else:ax.set_yscale('log');ax.set_ylabel('Exported PLA objects (log scale)')
    axes[0].legend(title='Literal label',ncol=2,fontsize=8,loc='lower right')
    return _finish(fig,'PLA-type contrasts within the same literal label',
        'Lines connect labels across separate exports. They do not establish the same tissue section, cell, or independent animal.\nDifferent antibody pairs and section coverage can affect counts and overlap; no significance tests are shown.')


def averaging_figure(type_table):
    fig,ax=plt.subplots(figsize=(9,4.5))
    y=np.arange(len(type_table))
    for key,shift,color,label in [('mean_export_association_pct',-.13,'#0072B2','Equal weight per export'),
                                 ('pooled_object_association_pct',.13,'#D55E00','Weight by PLA detections')]:
        ax.scatter(type_table[key],y+shift,c=color,label=label,s=70)
        for i,v in enumerate(type_table[key]):ax.annotate(f'{v:.2f}%',(v,i+shift),xytext=(8,0),textcoords='offset points',va='center')
    ax.set_yticks(y,type_table.pla_type);ax.set_xlim(0,100);ax.set_ylim(-.6,len(y)-.4)
    ax.set_xlabel('PLA objects with ChAT present (%)');ax.grid(axis='x',alpha=.2)
    ax.legend(loc='lower left',fontsize=9)
    return _finish(fig,'“Average” depends on the unit receiving equal weight',
        'Six exports per PLA type. These are descriptive summaries of the available exports;\nneither weighting estimates a running effect or an animal-level mean without a sample map.',bottom=.24)


def layout_figure(summary):
    images=list(summary.image_id.drop_duplicates())
    fig,axes=plt.subplots(1,len(images),figsize=(5.5*len(images),7),squeeze=False)
    for ax,image_id in zip(axes[0],images):
        g=summary.loc[summary.image_id==image_id]
        for row in g.itertuples():
            ax.add_patch(Rectangle((row.x_min,row.y_min),row.x_max-row.x_min,row.y_max-row.y_min,
                                  fill=False,edgecolor=TYPE_COLORS[str(row.pla_type)],lw=1.6))
            ax.text((row.x_min+row.x_max)/2,(row.y_min+row.y_max)/2,
                    f'{row.condition}\n{row.pla_type}\njob {row.job}',ha='center',va='center',fontsize=10)
        ax.set_xlim(g.x_min.min()-1000,g.x_max.max()+1000);ax.set_ylim(g.y_max.max()+1000,g.y_min.min()-1000)
        ax.set_aspect('equal');ax.set_title(str(g.image_label.iloc[0]));ax.set_xlabel('Native X (export units)')
        ax.set_ylabel('Native Y (export units)')
    return _finish(fig,'The exports occupy different positions in each source image',
        'Rectangles show measurement extents, not tissue-outline polygons. Each panel keeps its own source coordinates.\nSimilar slide layout does not establish corresponding anatomy; no translation, rotation, scaling or deformation has been applied.')


def native_grids(objects, bin_width=100, min_support=10):
    """Full-data histograms, with per-export edges and explicit support."""
    if not np.isfinite(bin_width) or bin_width<=0:raise ValueError('bin_width must be positive')
    if int(min_support)!=min_support or min_support<1:raise ValueError('min_support must be a positive integer')
    result={}
    for job,g in objects.groupby('job',observed=True,sort=True):
        if g.image_id.nunique()!=1:raise ValueError('A spatial grid cannot combine source images')
        edges=[]
        for lo,hi in [('xmin','xmax'),('ymin','ymax')]:
            start=np.floor(g[lo].min()/bin_width)*bin_width
            end=max(start+bin_width,np.ceil(g[hi].max()/bin_width)*bin_width)
            edges.append(np.arange(start,end+bin_width*.5,bin_width))
        xe,ye=edges
        def hist(frame,weights=None):
            return np.histogram2d(frame.y,frame.x,bins=(ye,xe),weights=weights)[0]
        pla=g.loc[g.object_type=='PLA'];chat=g.loc[g.object_type=='ChAT']
        n=hist(pla);k=hist(pla.loc[pla.chat_present]);c=hist(chat)
        fraction=np.divide(100*k,n,out=np.full_like(n,np.nan),where=n>=min_support)
        mean=np.divide(hist(pla,pla.overlap_pct),n,out=np.full_like(n,np.nan),where=n>=min_support)
        if int(n.sum())!=len(pla) or int(c.sum())!=len(chat):raise ValueError('Histogram lost objects')
        result[str(job)]={'job':str(job),'image_id':str(g.image_id.iloc[0]),'region':str(g.region.iloc[0]),
            'label':f'{g.condition.iloc[0]} · {g.pla_type.iloc[0]}', 'x_edges':xe,'y_edges':ye,
            'pla_count':n,'chat_count':c,'associated_count':k,'fraction':fraction,'mean_overlap':mean,
            'min_support':int(min_support),'bin_width':bin_width}
    return result


def grid_table(grids):
    frames=[]
    for job,g in grids.items():
        yy,xx=np.meshgrid(g['y_edges'][:-1],g['x_edges'][:-1],indexing='ij')
        frames.append(pd.DataFrame({'job':job,'image_id':g['image_id'],'x_left':xx.ravel(),'y_top':yy.ravel(),
            'bin_width':g['bin_width'],'pla_count':g['pla_count'].ravel().astype(int),
            'chat_count':g['chat_count'].ravel().astype(int),'associated_pla_count':g['associated_count'].ravel().astype(int),
            'pla_association_pct':g['fraction'].ravel(),'mean_overlap_pct':g['mean_overlap'].ravel(),
            'supported':(g['pla_count']>=g['min_support']).ravel()}))
    return pd.concat(frames,ignore_index=True)


def shared_count_scales(grids):
    return {field:max(2,max(float(g[field].max()) for g in grids.values())) for field in ['pla_count','chat_count']}


def _extent(g):return [g['x_edges'][0],g['x_edges'][-1],g['y_edges'][-1],g['y_edges'][0]]


def _count(ax,g,field,cmap,vmax):
    if field == 'pla_count' and cmap == 'inferno':
        # Keep sparse PLA bins visibly colored against the grayscale background.
        cmap = ListedColormap(plt.get_cmap('inferno')(np.linspace(.32, 1, 256)))
    return ax.imshow(np.ma.masked_less(g[field],1),extent=_extent(g),origin='upper',interpolation='nearest',
                     cmap=cmap,norm=LogNorm(vmin=1,vmax=vmax),aspect='equal')


def section_figure(count_grid, fraction_grid, scales):
    """Four panels: background, fine overlay, supported fraction, denominator."""
    c,f=count_grid,fraction_grid
    fig,axes=plt.subplots(1,4,figsize=(20,6.2))
    m=_count(axes[0],c,'chat_count','gray_r',scales['chat_count'])
    fig.colorbar(m,ax=axes[0],orientation='horizontal',pad=.13,fraction=.05).set_label('ChAT objects / bin (log)')
    _count(axes[1],c,'chat_count','gray_r',scales['chat_count'])
    m=_count(axes[1],c,'pla_count','inferno',scales['pla_count'])
    fig.colorbar(m,ax=axes[1],orientation='horizontal',pad=.13,fraction=.05).set_label('PLA objects / bin (log)')
    low=(f['pla_count']>0)&(f['pla_count']<f['min_support'])
    axes[2].imshow(np.ma.masked_where(~low,np.ones_like(f['pla_count'])),extent=_extent(f),origin='upper',
                   interpolation='nearest',cmap='Greys',vmin=0,vmax=2,aspect='equal')
    m=axes[2].imshow(np.ma.masked_invalid(f['fraction']),extent=_extent(f),origin='upper',interpolation='nearest',
                     cmap='viridis',vmin=0,vmax=100,aspect='equal')
    fig.colorbar(m,ax=axes[2],orientation='horizontal',pad=.13,fraction=.05).set_label('PLA with ChAT (%)')
    m=_count(axes[3],f,'pla_count','cividis',max(2,float(f['pla_count'].max())))
    fig.colorbar(m,ax=axes[3],orientation='horizontal',pad=.13,fraction=.05).set_label('Raw PLA support (log)')
    supported=f['pla_count']>=f['min_support'];retained=100*f['pla_count'][supported].sum()/f['pla_count'].sum()
    for ax,title in zip(axes,['ChAT coordinate reconstruction','PLA counts over ChAT background',
                              f'Local PLA association (N ≥ {f["min_support"]})','Denominator for the fraction map']):
        ax.set_title(title,fontsize=11);ax.set_xlabel('Native X');ax.set_ylabel('Native Y');ax.tick_params(labelsize=8)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.set_facecolor('white')
    return _finish(fig,f'{c["label"]} · job {c["job"]}: native-coordinate reconstruction',
        f'All objects contribute. Fine bins: {c["bin_width"]:g} units; fraction bins: {f["bin_width"]:g} units. Supported bins cover {retained:.1f}% of PLA detections.\n'
        'Gray fraction bins have too few PLA objects; white bins have none. Count scales are shared across the 12 sections; support scale is local.\n'
        'Coordinates are box midpoints in export units. No microscopy pixels, tissue contours, smoothing, physical scale or anatomical registration are inferred.',bottom=.26)


def fraction_atlas(grids, summary):
    labels=sorted(summary.condition.unique());fig,axes=plt.subplots(2,len(labels),figsize=(3.7*len(labels),9),squeeze=False)
    for row,typ in enumerate(['Hk1','Pfkp']):
        for col,label in enumerate(labels):
            job=str(summary.loc[(summary.condition==label)&(summary.pla_type==typ),'job'].iloc[0]);g=grids[job];ax=axes[row,col]
            low=(g['pla_count']>0)&(g['pla_count']<g['min_support'])
            ax.imshow(np.ma.masked_where(~low,np.ones_like(g['pla_count'])),origin='upper',extent=_extent(g),
                      cmap='Greys',vmin=0,vmax=2,interpolation='nearest',aspect='equal')
            m=ax.imshow(np.ma.masked_invalid(g['fraction']),origin='upper',extent=_extent(g),vmin=0,vmax=100,
                        cmap='viridis',interpolation='nearest',aspect='equal')
            ax.set_title(f'{label} · {typ}\njob {job}',fontsize=10);ax.tick_params(labelsize=7)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.set_xlabel('Native X',fontsize=8);ax.set_ylabel('Native Y',fontsize=8)
    fig.subplots_adjust(left=.045,right=.96,top=.85,bottom=.22,wspace=.35,hspace=.5)
    cax=fig.add_axes([.30,.12,.40,.022]);fig.colorbar(m,cax=cax,orientation='horizontal').set_label('PLA objects with ChAT present (%)')
    fig.suptitle('All twelve sections: supported local fractions on a common color scale',fontsize=18,x=.025,ha='left')
    fig.text(.025,.025,'Each panel uses its own native coordinates and extent. Panels are not anatomically aligned.\nGray = nonempty but below support; white = no PLA observations. Individual section figures below provide the full-resolution view.',fontsize=10)
    return fig


def support_sensitivity(objects, widths=(500,1000,2000), thresholds=(5,10,20)):
    rows=[]
    for width in widths:
        grids=native_grids(objects,bin_width=width,min_support=1)
        for job,g in grids.items():
            n=g['pla_count']
            for threshold in thresholds:
                rows.append(dict(job=job,label=g['label'],bin_width=width,min_support=threshold,
                    occupied_bins=int((n>0).sum()),visible_bins=int((n>=threshold).sum()),
                    visible_pla_pct=100*n[n>=threshold].sum()/n.sum()))
    return pd.DataFrame(rows)


def support_figure(table):
    labels=table[['job','label']].drop_duplicates().sort_values('label');widths=sorted(table.bin_width.unique())
    fig,axes=plt.subplots(1,len(widths),figsize=(16,7),sharey=True)
    for ax,width in zip(axes,widths):
        for threshold,color in zip(sorted(table.min_support.unique()),['#56B4E9','#0072B2','#D55E00']):
            t=table.loc[(table.bin_width==width)&(table.min_support==threshold)].set_index('job').loc[labels.job]
            ax.plot(t.visible_pla_pct,np.arange(len(t)),'o-',label=f'N ≥ {threshold}',color=color,lw=.8,markersize=4)
        ax.set_yticks(np.arange(len(labels)),labels.label);ax.invert_yaxis();ax.set_xlim(0,101)
        ax.set_title(f'{width:g}-unit bins');ax.set_xlabel('PLA objects in visible bins (%)');ax.grid(axis='x',alpha=.2)
    axes[0].legend(fontsize=8,loc='lower left')
    return _finish(fig,'Resolution–support tradeoff: how much of each section remains visible?',
        'Finer spatial bins contain fewer objects. Visibility thresholds control display support, not statistical significance.\nUse these curves alongside the atlas; a changing map can reflect binning and support rather than a biological change.')


def distribution_figure(objects):
    fig,axes=plt.subplots(2,3,figsize=(17,9));frames=[]
    labels=sorted(objects.condition.unique());colors=dict(zip(labels,plt.cm.tab10(np.linspace(0,.8,len(labels)))))
    for ri,typ in enumerate(['Hk1','Pfkp']):
        for ci,(metric,xlabel) in enumerate([('area','PLA area (µm², as exported; log axis)'),('intensity','PLA intensity (export units; log axis)'),('overlap_pct','PLA area overlapping ChAT (%)')]):
            ax=axes[ri,ci]
            for (label,job),g in objects.loc[(objects.object_type=='PLA')&(objects.pla_type==typ)].groupby(['condition','job'],observed=True):
                values,counts=np.unique(g[metric],return_counts=True);cdf=100*counts.cumsum()/counts.sum()
                # symlog preserves possible zero intensity observations.
                ax.step(values,cdf,where='post',color=colors[label],label=str(label),lw=1.2)
                frames.append(pd.DataFrame(dict(job=str(job),metric=metric,value=values,object_count=counts,cumulative_pct=cdf)))
            if metric!='overlap_pct':
                ax.set_xscale('symlog',linthresh=1e-6) if metric=='intensity' else ax.set_xscale('log')
            else:ax.set_xlim(0,100)
            ax.set_ylim(0,100);ax.set_title(f'{typ}–TDP43: {metric}');ax.set_xlabel(xlabel);ax.set_ylabel('PLA at or below value (%)');ax.grid(alpha=.15)
    axes[0,0].legend(title='Literal label',fontsize=8,ncol=2)
    _finish(fig,'Read the full PLA distributions, separately for each assay',
        'Exact cumulative distributions preserve ties and all exported objects. The intensity axis includes zero if present.\nIntensity comparability across images is unverified. No object-level confidence intervals imply animal-level precision.',bottom=.16)
    fig.subplots_adjust(hspace=.45)
    return fig,pd.concat(frames,ignore_index=True)


def quality_figure(summary):
    s=summary.sort_values(['condition','pla_type']);labels=s.condition.astype(str)+' · '+s.pla_type.astype(str)
    fig,axes=plt.subplots(1,3,figsize=(17,7),sharey=True)
    for ax,key,title in zip(axes,['any_zero_diameter_pct','zero_bbox_pct','largest_1pct_area_share_pct'],
                           ['PLA with a zero diameter','PLA with a collapsed bounding box','PLA area carried by largest 1% of objects']):
        ax.barh(np.arange(len(s)),s[key],color=[TYPE_COLORS[str(t)] for t in s.pla_type])
        for i,value in enumerate(s[key]):ax.text(value+.15,i,f'{value:.1f}%',va='center',fontsize=8)
        ax.set_yticks(np.arange(len(s)),labels);ax.set_title(title,fontsize=11);ax.set_xlabel('Percentage (%)');ax.grid(axis='x',alpha=.15)
        ax.set_xlim(0,max(5,s[key].max()*1.20))
    axes[0].invert_yaxis()
    return _finish(fig,'Measurement review: geometry flags and large-object influence',
        'Zero geometry is retained and requires image review. Large objects are not automatically artifacts.\nArea concentration helps explain why area-weighted overlap can differ from the typical-object result.')


def interrupted_figure(short_grid,full_grid):
    # Grid edges must be identical to show a genuine spatial difference.
    for key in ['x_edges','y_edges']:
        if not np.array_equal(short_grid[key],full_grid[key]):raise ValueError('Rerun grids need identical edges')
    a=short_grid['chat_count']+short_grid['pla_count'];b=full_grid['chat_count']+full_grid['pla_count'];delta=b-a
    fig,axes=plt.subplots(1,3,figsize=(16,6.7));maximum=max(2,b.max(),a.max())
    for ax,array,title in zip(axes,[a,b,delta],['Interrupted job 4978','Replacement job 5403','Replacement minus interrupted']):
        if ax is axes[2]:
            bound=max(1,np.abs(delta).max());m=ax.imshow(array,extent=_extent(full_grid),origin='upper',interpolation='nearest',
                  cmap='RdBu_r',norm=TwoSlopeNorm(vcenter=0,vmin=-bound,vmax=bound),aspect='equal')
        else:m=ax.imshow(np.ma.masked_less(array,1),extent=_extent(full_grid),origin='upper',interpolation='nearest',
                        cmap='cividis',norm=LogNorm(1,maximum),aspect='equal')
        ax.set_title(title);ax.set_xlabel('Native X');ax.set_ylabel('Native Y')
        fig.colorbar(m,ax=ax,orientation='horizontal',pad=.13,fraction=.05).set_label('All objects / bin'+(' (log)' if ax is not axes[2] else ' (signed)'))
    return _finish(fig,'Interrupted-run audit: 22F · Pfkp–TDP43',
        'Both object types included for coverage review. Identical coordinate bins and shared count scales.\nOnly the replacement contributes to the main analysis. Exact row-signature matching is reported in the accompanying table.',bottom=.24)


def aligned_rerun_grids(short,full,bin_width=250):
    result=[]
    base=native_grids(full,bin_width=bin_width,min_support=1)[str(full.job.iloc[0])]
    for objects in [short,full]:
        g=base.copy()
        for typ,field in [('PLA','pla_count'),('ChAT','chat_count')]:
            part=objects.loc[objects.object_type==typ]
            g[field]=np.histogram2d(part.y,part.x,bins=(g['y_edges'],g['x_edges']))[0]
            if int(g[field].sum())!=len(part):raise ValueError('Rerun extends beyond replacement frame')
        result.append(g)
    return tuple(result)


def animal_figure(animal,groups):
    regions=list(animal.anatomical_region.unique());fig,axes=plt.subplots(1,len(regions),figsize=(6*len(regions),5),squeeze=False)
    for ax,region in zip(axes[0],regions):
        part=animal.loc[animal.anatomical_region==region];labels=sorted(part.running_group.unique())
        for i,label in enumerate(labels):
            for typ,shift in [('Hk1',-.12),('Pfkp',.12)]:
                g=part.loc[(part.running_group==label)&(part.pla_type==typ)]
                ax.scatter(np.full(len(g),i+shift),g.association_pct,color=TYPE_COLORS[typ],alpha=.7,label=typ if i==0 else None)
                if len(g):ax.plot(i+shift,g.association_pct.mean(),marker='_',ms=20,mew=3,color='black')
        ax.set_xticks(range(len(labels)),labels);ax.set_ylim(0,100);ax.set_ylabel('Animal mean PLA association (%)');ax.set_title(region);ax.legend()
    return _finish(fig,'Running-group summaries from confirmed sample metadata',
        'Each dot is an animal after averaging its sections within anatomy and PLA type. Black marks show group means.\nNo significance test is fitted; sample design, group balance and assay comparability still require review.')


def area_sensitivity_figure(table):
    fig,axes=plt.subplots(2,2,figsize=(13,9));labels=sorted(table.condition.unique())
    colors=dict(zip(labels,plt.cm.tab10(np.linspace(0,.8,len(labels)))))
    for row,typ in enumerate(['Hk1','Pfkp']):
        for col,metric,title in [(0,'retained_pct','Objects retained'),(1,'association_pct','PLA with ChAT present')]:
            ax=axes[row,col]
            for label,g in table.loc[table.pla_type==typ].groupby('condition',observed=True):
                ax.plot(g.minimum_area_um2,g[metric],'o-',color=colors[label],label=label,ms=4,lw=1)
            ax.set_xscale('symlog',linthresh=.03);ax.set_ylim(0,100)
            ax.set_title(f'{typ}: {title}');ax.set_xlabel('Minimum PLA area (µm², as exported)');ax.set_ylabel('Percentage (%)');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8,ncol=2)
    _finish(fig,'Size-filter sensitivity: results depend on which detections remain',
        'Illustrative inclusive area thresholds; no threshold is selected or applied to the primary results.\nReview retention and association together. A shifted fraction after filtering does not establish improved assay accuracy.')
    fig.subplots_adjust(hspace=.45)
    return fig


def workflow_figure():
    fig,ax=plt.subplots(figsize=(16,6));ax.set_xlim(0,16);ax.set_ylim(0,6);ax.axis('off')
    boxes=[(.2,3.5,'Tissue + fluorescence imaging\nExperimental metadata needed'),(4.2,3.5,'HALO object exports\nInitial: 8  |  Expanded: 13'),
           (8.2,3.5,'Validate + select analyses\n4 copied files; 1 interrupted job'),(12.2,3.5,'12 primary analyses\n3 source images; 6 literal labels'),
           (12.2,.8,'Describe each section\nCounts · overlap · native maps'),(8.2,.8,'Confirm sample assignments\nAnimal · running · section · anatomy'),
           (4.2,.8,'Validate anatomical registration\nImages · landmarks · masks'),(.2,.8,'Answer the scientific question\nAnimal summaries or spatial contrasts')]
    for x,y,text in boxes:
        ax.add_patch(Rectangle((x,y),3.5,1.2,facecolor='#f4f6f8',edgecolor='#657487',lw=1.2))
        ax.text(x+1.75,y+.6,text,ha='center',va='center',fontsize=10)
    for a,b in [((3.7,4.1),(4.2,4.1)),((7.7,4.1),(8.2,4.1)),((11.7,4.1),(12.2,4.1)),((13.95,3.5),(13.95,2)),
                ((12.2,1.4),(11.7,1.4)),((8.2,1.4),(7.7,1.4)),((4.2,1.4),(3.7,1.4))]:
        ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle='->',color='#34465b',lw=1.3))
    # Section/animal averages have a direct route that bypasses registration.
    ax.plot([9.95,9.95,1.95],[2,2.75,2.75],color='#34465b',lw=1.3)
    ax.annotate('',xy=(1.95,2),xytext=(1.95,2.75),arrowprops=dict(arrowstyle='->',color='#34465b',lw=1.3))
    ax.text(5.8,2.88,'Section / animal summaries: no registration needed',ha='center',fontsize=9)
    ax.text(6,.54,'Spatial comparisons: validate corresponding anatomy',ha='center',fontsize=9)
    ax.text(8,5.4,'Scientific workflow: from exported measurements to testable comparisons',ha='center',fontsize=18,fontweight='bold')
    ax.text(8,.15,'Registration is needed for anatomical spatial comparisons; whole-section summaries require comparable sampling and metadata, but no image warping.',ha='center',fontsize=10)
    fig.tight_layout()
    return fig

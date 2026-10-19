"""Explicit denominators, two subtraction directions and animal-level checks."""
import numpy as np
import pandas as pd
from .condition_comparison import HOME, RUNNING, PRIMARY


def audit(results):
    a=results['animal_metrics'].copy(); rows=[]; directions=[]
    for arm,part in a.groupby('pla_arm',observed=True):
        h=part[part.running_group.eq(HOME)];r=part[part.running_group.eq(RUNNING)]
        primary=results['contrasts'].query('pla_arm == @arm and metric == @PRIMARY').iloc[0]
        standard=results['standardization'].query("pla_arm == @arm and estimator == 'standardized'").iloc[0]
        estimates=[('PLA-count normalized · equal animals',primary.home_mean,primary.running_mean,'Primary: 100 × positive PLA / all PLA per animal, then mean of three animals'),
                   ('Common PLA size mix · equal animals',standard.home_mean,standard.running_mean,'Same five size weights per animal; descriptive sensitivity analysis'),
                   ('Actual overlap / PLA area · equal animals',h.overlap_area_ratio_pct.mean(skipna=False),r.overlap_area_ratio_pct.mean(skipna=False),'100 × summed overlap area / summed PLA area per animal; a different estimand'),
                   ('Pooled PLA objects · diagnostic only',100*h.n_positive.sum(min_count=3)/h.pla_count.sum(),100*r.n_positive.sum(min_count=3)/r.pla_count.sum(),'Counterfactual object weighting: large exports dominate; not the biological primary estimate')]
        for name,hv,rv,definition in estimates:
            rows.append(dict(pla_arm=arm,normalization=name,home_mean=hv,running_mean=rv,home_minus_running=hv-rv,running_minus_home=rv-hv,definition=definition))
        hv=h.association_pct.to_numpy();rv=r.association_pct.to_numpy()
        # n=3 biological animals per group, never object-level binomial SE.
        se=np.sqrt(np.var(hv,ddof=1)/len(hv)+np.var(rv,ddof=1)/len(rv)) if np.isfinite(np.r_[hv,rv]).all() else np.nan
        for direction,first,second,lo,hi in [
            ('Home cage − Running',primary.home_mean,primary.running_mean,primary.leave_pair_out_min,primary.leave_pair_out_max),
            ('Running − Home cage',primary.running_mean,primary.home_mean,-primary.leave_pair_out_max,-primary.leave_pair_out_min)]:
            delta=first-second
            directions.append(dict(pla_arm=arm,direction=direction,first_mean=first,reference_mean=second,difference_pp=delta,
                                   relative_change_pct=100*delta/second if second>0 else np.nan,
                                   ratio=first/second if second>0 else np.nan,
                                   between_animal_se_pp=se,difference_over_se=delta/se if se>0 else np.nan,
                                   leave_pair_out_min=lo,leave_pair_out_max=hi))
    # Per-animal paired assay comparison, not paired anatomy or group pairing.
    arms=list(sorted(a.pla_arm.unique()))
    paired=a.pivot(index=['animal_id','running_group'],columns='pla_arm',values='association_pct').reset_index()
    paired['pfkp_minus_hk1_pp']=paired[arms[1]]-paired[arms[0]]
    denom=a[['pla_arm','animal_id','running_group','n_positive','pla_count','total_overlap_area_um2','total_pla_area_um2','association_pct','overlap_area_ratio_pct']].copy()
    denom['equal_animal_weight_pct']=100/3
    denom['pooled_object_weight_pct']=100*denom.pla_count/denom.groupby(['pla_arm','running_group']).pla_count.transform('sum')
    return dict(normalization=pd.DataFrame(rows),directions=pd.DataFrame(directions),denominators=denom,paired_arms=paired)


def normalization_figure(results):
    import matplotlib.pyplot as plt
    data=results['normalization'];fig,axes=plt.subplots(1,2,figsize=(16,6.5),sharey=True)
    labels=['Per-animal PLA fraction\nPRIMARY','Common size mix','Actual overlap / PLA area','Pooled objects\nWEIGHTING DIAGNOSTIC']
    for ax,(arm,part) in zip(axes,data.groupby('pla_arm',observed=True)):
        y=np.arange(len(part));d=part.home_minus_running.to_numpy()
        ax.barh(y-.14,d,height=.25,color='#087e8b',label='Home cage − Running')
        ax.barh(y+.14,-d,height=.25,color='#d06738',label='Running − Home cage')
        for i,v in enumerate(d):
            ax.text(v+.5,i-.14,f'+{v:.2f}',va='center',fontsize=10);ax.text(-v-.5,i+.14,f'{-v:.2f}',ha='right',va='center',fontsize=10)
        ax.set(yticks=y,yticklabels=labels,xlim=(-34,34),title=arm,xlabel='Difference · percentage points')
        ax.set_ylim(3.6,-.6);ax.axvline(0,color='#173447',lw=1);ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
    axes[0].legend(loc='lower center',bbox_to_anchor=(1.05,-.31),ncol=2,frameon=False)
    fig.suptitle('Is it normalized? Yes to PLA count; not to tissue area or ChAT availability',fontsize=19,fontweight='bold')
    fig.text(.04,.025,'Each bar is an observed descriptive contrast, not a significance test. Reversing subtraction reverses its sign only.\nThe overlap-area row asks a different question. Pooled objects are shown only to expose unequal animal weights.',fontsize=11)
    fig.tight_layout(rect=(0,.12,1,.92));return fig


def paired_arm_figure(results):
    import matplotlib.pyplot as plt
    p=results['paired_arms'];arms=[c for c in p if c.startswith('TDP43')]
    fig,axes=plt.subplots(1,2,figsize=(14,6))
    for _,r in p.iterrows():
        color='#087e8b' if r.running_group==HOME else '#d06738'
        axes[0].plot([0,1],[r[arms[0]],r[arms[1]]],'-o',color=color,alpha=.75,lw=1.6)
    last=-np.inf
    for _,r in p.sort_values(arms[1]).iterrows():
        label_y=max(float(r[arms[1]]),last+4.3);last=label_y
        color='#087e8b' if r.running_group==HOME else '#d06738'
        axes[0].annotate(r.animal_id,(1,r[arms[1]]),xytext=(1.09,label_y),fontsize=10,color=color,
                         va='center',arrowprops=dict(arrowstyle='-',color=color,lw=.6))
    axes[0].set(xticks=[0,1],xticklabels=['HK1 arm','PFKP arm'],xlim=(-.15,1.25),ylabel='PLA with ChAT signal (%)',ylim=(0,100),title='Same animals, different physical sections')
    for i,group in enumerate([HOME,RUNNING]):
        v=p[p.running_group.eq(group)];color=['#087e8b','#d06738'][i]
        axes[1].scatter(np.array([i]*3)+[-.12,0,.12],v.pfkp_minus_hk1_pp,color=color,s=55)
        axes[1].hlines(v.pfkp_minus_hk1_pp.mean(),i-.3,i+.3,color=color,lw=4)
        for j,r in enumerate(v.itertuples()):axes[1].annotate(r.animal_id,(i+[-.12,0,.12][j],r.pfkp_minus_hk1_pp),xytext=(0,8),textcoords='offset points',ha='center',fontsize=10)
    hd=p[p.running_group.eq(HOME)].pfkp_minus_hk1_pp.mean();rd=p[p.running_group.eq(RUNNING)].pfkp_minus_hk1_pp.mean()
    axes[1].set(xticks=[0,1],xticklabels=[HOME,RUNNING],xlim=(-.5,1.5),ylim=(0,35),ylabel='PFKP-arm fraction − HK1-arm fraction (pp)',title=f'Between-condition difference of arm gaps: {hd-rd:+.2f} pp')
    fig.suptitle('Both PLA arms show a similar condition gap',fontsize=20,fontweight='bold')
    fig.text(.04,.025,'Lines link animal identities, not spatially matched tissue. Similar gaps motivate checking shared ChAT/detection effects;\nthey do not prove a shared mechanism or equivalent assays. Three animals per condition limit arm-selectivity claims.',fontsize=11)
    fig.tight_layout(rect=(0,.10,1,.91));return fig

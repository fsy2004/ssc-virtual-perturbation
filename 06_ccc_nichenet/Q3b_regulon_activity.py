#!/usr/bin/env python
"""Q3b: HES1 regulon ACTIVITY per cell type (decoupler ULM + CollecTRI), lognorm input.
Positive control SMAD3, negative control MEF2C. Mainline-consistent readout."""
import os, json, time
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ.setdefault(k,'4')
import numpy as np, pandas as pd, anndata as ad, h5py
import decoupler as dc
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
try: from anndata.io import read_elem
except Exception: from anndata.experimental import read_elem
np.random.seed(0)

FULL='/data/ssc/powered/full_integrated_annotated.h5ad'
FIB='/data/ssc/powered/fibroblast_integrated.h5ad'
COND='/data/ssc/powered/sample_condition.csv'
NET='/data/ssc/data/collectri_net.csv'
AQ='/data/ssc/atlas_queries'; AF=AQ+'/figures'; os.makedirs(AF, exist_ok=True)
TFS=['HES1','SMAD3','MEF2C']; RECEIVER='Myofibroblast'
def log(*a): print('[Q3b]',*a,flush=True)

# ---- metadata: ct_work + condition ----
obs=ad.read_h5ad(FULL, backed='r').obs.copy()
fibsub=ad.read_h5ad(FIB, backed='r').obs['fib_subtype'].astype(str)
cond=pd.read_csv(COND); s2c=dict(zip(cond['sample'].astype(str), cond['condition'].astype(str)))
ctw=obs['celltype'].astype(str).copy()
common=obs.index.intersection(fibsub.index); ctw.loc[common]=fibsub.loc[common].values
ctw[(obs['celltype'].astype(str)=='Fibroblast') & (~obs.index.isin(common))]='Fibroblast_other'
obs['ct_work']=ctw.values; obs['condition']=obs['sample'].astype(str).map(s2c).fillna('unknown')

# ---- lognorm X ----
log('loading lognorm X...'); t0=time.time()
with h5py.File(FULL,'r') as f:
    X=read_elem(f['X']); var=read_elem(f['var'])
vn=np.asarray(var.index if hasattr(var,'index') else var['_index'])
adata=ad.AnnData(X=X, obs=obs[['ct_work','condition']].copy(), var=pd.DataFrame(index=vn))
adata.var_names_make_unique()
log('adata', adata.shape, '%.0fs'%(time.time()-t0))

# ---- OOM fix: stratified per-cell-type subsample (per-celltype mean is robust; dense ULM on 424k x 37556 densifies to ~88GB -> OOM-killed twice, dmesg-confirmed) ----
CAP=5000
_rng=np.random.RandomState(0)
_ct=obs['ct_work'].values
_idx=[]
for _c in pd.unique(_ct):
    _ii=np.where(_ct==_c)[0]
    if len(_ii)>CAP: _ii=_rng.choice(_ii, CAP, replace=False)
    _idx.append(_ii)
_idx=np.sort(np.concatenate(_idx))
adata=adata[_idx].copy()
obs=obs.iloc[_idx].copy()
log('subsampled to', adata.shape, 'cap/ct=%d'%CAP, 'n_ct=%d'%len(pd.unique(_ct)))


# ---- CollecTRI net, filtered to the 3 TFs ----
net=pd.read_csv(NET)[['source','target','weight']]
net=net[net['source'].isin(TFS)].copy()
log('net targets:', net.groupby('source').size().to_dict())

# ---- ULM regulon activity (per cell) ----
dc.mt.ulm(adata, net, tmin=5, verbose=True)
act=adata.obsm['score_ulm'].copy()            # cells x TFs (t-values)
act.columns=[str(c) for c in act.columns]
present=[t for t in TFS if t in act.columns]
log('ULM done; TFs scored:', present)
act['ct_work']=obs['ct_work'].values; act['condition']=obs['condition'].values

# ---- per cell type (all conditions) ----
rows=[]
for tf in present:
    g=act.groupby('ct_work')[tf].agg(mean='mean', median='median', std='std', n='count').reset_index()
    g['TF']=tf; rows.append(g)
per_ct=pd.concat(rows, ignore_index=True)
per_ct.to_csv(AQ+'/Q3b_HES1_regulon_activity_by_celltype.csv', index=False)

# ---- per cell type x condition ----
rows2=[]
for tf in present:
    g=act.groupby(['ct_work','condition'])[tf].agg(mean='mean', median='median', n='count').reset_index()
    g['TF']=tf; rows2.append(g)
per_ctcond=pd.concat(rows2, ignore_index=True)
per_ctcond.to_csv(AQ+'/Q3b_regulon_activity_by_celltype_condition.csv', index=False)

# ---- ranks / sanity ----
def rank_of(tf, ct):
    d=per_ct[per_ct.TF==tf].sort_values('mean', ascending=False).reset_index(drop=True)
    r=d.index[d.ct_work==ct]
    return (int(r[0])+1, len(d), float(d.loc[r[0],'mean'])) if len(r) else (None,len(d),None)
KEY={}
for tf in present:
    r,n,val=rank_of(tf, RECEIVER)
    top=per_ct[per_ct.TF==tf].sort_values('mean',ascending=False).head(5)[['ct_work','mean']]
    KEY[tf]={'myofib_rank':r,'n_celltypes':n,'myofib_mean_activity':round(val,3) if val is not None else None,
             'top5':[(x.ct_work, round(x.mean,3)) for x in top.itertuples()]}
# HES1 SSc vs HC in myofib
h=per_ctcond[(per_ctcond.TF=='HES1')&(per_ctcond.ct_work==RECEIVER)].set_index('condition')['mean']
KEY['HES1_myofib_SSc_vs_HC']={c: round(float(h[c]),3) for c in h.index if c in ['SSc','HC']}
json.dump(KEY, open(AQ+'/Q3b_KEY.json','w'), indent=2, default=str)
log('KEY:\n'+json.dumps(KEY, indent=2, default=str))

# ---- figure: 3-panel lollipop (HES1 | SMAD3 | MEF2C), myofib highlighted ----
order_cts=['Myofibroblast','SFRP4_proFib','Adipogenic','SFRP2_DPP4','FMO1_LSP1','LGR5_Gur','Inflammatory',
           'Fibroblast_other','Endothelial','Pericyte','SmoothMuscle','Myeloid','Lymphatic','Mast','Tcell',
           'NK','Bcell','Plasma','pDC','Keratinocyte','Melanocyte']
fig,axes=plt.subplots(1, len(present), figsize=(4.2*len(present), 6.4), sharey=True)
if len(present)==1: axes=[axes]
for ax,tf in zip(axes, present):
    d=per_ct[per_ct.TF==tf].set_index('ct_work')
    cts=[c for c in order_cts if c in d.index]
    y=np.arange(len(cts))[::-1]
    vals=[d.loc[c,'mean'] for c in cts]
    cols=['#C44E52' if c==RECEIVER else ('#4C72B0' if c in ['SFRP4_proFib','Adipogenic','SFRP2_DPP4','FMO1_LSP1','LGR5_Gur','Inflammatory','Fibroblast_other'] else '#999999') for c in cts]
    ax.hlines(y, 0, vals, color=cols, lw=1.8)
    ax.scatter(vals, y, color=cols, s=42, zorder=3)
    ax.axvline(0, color='k', lw=0.6, ls='-')
    ax.set_yticks(y); ax.set_yticklabels(cts, fontsize=7.5)
    role={'HES1':'(mainline effector)','SMAD3':'(positive control)','MEF2C':'(negative control)'}.get(tf,'')
    ax.set_title(f'{tf} regulon activity\n{role}', fontsize=9)
    ax.set_xlabel('ULM activity (mean t-value)'); ax.grid(True, axis='x', ls=':', alpha=0.4)
axes[0].scatter([],[],color='#C44E52',label='Myofibroblast'); axes[0].scatter([],[],color='#4C72B0',label='other fibroblast'); axes[0].scatter([],[],color='#999999',label='non-fibroblast')
axes[0].legend(fontsize=7, loc='lower right')
fig.suptitle('HES1 regulon activity vs SMAD3 (pos) / MEF2C (neg) controls  -  CollecTRI + decoupler ULM', fontsize=10)
plt.tight_layout(); plt.savefig(AF+'/Q3b_regulon_activity_lollipop.png', dpi=200, bbox_inches='tight'); plt.savefig(AF+'/Q3b_regulon_activity_lollipop.pdf', bbox_inches='tight'); plt.close()
log('figure saved'); log('Q3b_DONE')

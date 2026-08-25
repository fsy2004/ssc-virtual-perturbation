#!/usr/bin/env python
"""Corrected niche->myofibroblast Notch CCC: prep + LIANA (SSc vs HC)."""
import os, sys, json, time
os.environ.setdefault('OMP_NUM_THREADS','4'); os.environ.setdefault('MKL_NUM_THREADS','4')
os.environ.setdefault('OPENBLAS_NUM_THREADS','4'); os.environ.setdefault('NUMBA_NUM_THREADS','4')
os.environ.setdefault('NUMEXPR_NUM_THREADS','4')
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
import liana as li

OUT='/data/ssc/ccc_nichenet'
FULL='/data/ssc/powered/full_integrated_annotated.h5ad'
FIB='/data/ssc/powered/fibroblast_integrated.h5ad'
COND='/data/ssc/powered/sample_condition.csv'
os.makedirs(OUT+'/liana', exist_ok=True); os.makedirs(OUT+'/export', exist_ok=True)

SEED=0; CAP=8000; EXPR_PROP=0.10; N_PERMS=1000
NICHE=['Endothelial','Pericyte','SmoothMuscle','Myeloid','Lymphatic','Mast']
RECEIVER='Myofibroblast'
NOTCH_LIG=['DLL1','DLL3','DLL4','JAG1','JAG2']
NOTCH_REC=['NOTCH1','NOTCH2','NOTCH3','NOTCH4']
PROGRAM=['ACTA2','TAGLN','POSTN','COL1A1','COL1A2','COL3A1','COL11A1','COMP','CTHRC1','FN1']
NOTCH_TGT=['HES1','HEY1','HEYL','HES5','NRARP','DTX1']
def log(*a): print('[prep]',*a,flush=True)

# ---- metadata (backed) ----
log('reading obs (backed)...')
full_b=ad.read_h5ad(FULL, backed='r')
obs=full_b.obs.copy()
fib_b=ad.read_h5ad(FIB, backed='r')
fibsub=fib_b.obs['fib_subtype'].astype(str)
cond=pd.read_csv(COND); s2c=dict(zip(cond['sample'].astype(str), cond['condition'].astype(str)))

ctw=obs['celltype'].astype(str).copy()
common=obs.index.intersection(fibsub.index)
ctw.loc[common]=fibsub.loc[common].values
fibmask=(obs['celltype'].astype(str)=='Fibroblast') & (~obs.index.isin(common))
ctw[fibmask]='Fibroblast_other'
obs['ct_work']=ctw.values
obs['condition']=obs['sample'].astype(str).map(s2c).fillna('unknown')
log('ct_work counts:\n'+obs['ct_work'].value_counts().to_string())
log('condition counts:\n'+obs['condition'].value_counts().to_string())
# cross tab receiver/niche by condition
ct_tab=pd.crosstab(obs['ct_work'], obs['condition'])
ct_tab.to_csv(OUT+'/export/ctwork_by_condition.csv')
log('crosstab saved')

# ---- pick cells: SSc/HC only, cap large types, keep receiver+niche full ----
rng=np.random.default_rng(SEED)
KEEP_FULL=set([RECEIVER])|set(NICHE)
keep=[]
for cn in ['SSc','HC']:
    sub_obs=obs[obs['condition']==cn]
    for ct, g in sub_obs.groupby('ct_work'):
        idx=g.index.to_numpy()
        if len(idx)==0: continue
        if ct in KEEP_FULL or len(idx)<=CAP: keep.append(idx)
        else: keep.append(rng.choice(idx, CAP, replace=False))
keep=np.concatenate(keep); keep_set=set(keep.tolist())
log('total kept cells:', len(keep))
mask=obs.index.isin(keep_set)

log('loading subset into memory (X + obs only)...')
t0=time.time()
sub_full=full_b[mask.values if hasattr(mask,'values') else mask].to_memory()
X=sub_full.X
sub=ad.AnnData(X=X, obs=sub_full.obs.copy(), var=sub_full.var.copy())
del sub_full
sub.obs['ct_work']=obs.loc[sub.obs_names,'ct_work'].values
sub.obs['condition']=obs.loc[sub.obs_names,'condition'].values
log('subset shape', sub.shape, 'in %.1fs'%(time.time()-t0))
# sanity: lognorm?
d=sub.X[:500].data if sp.issparse(sub.X) else np.asarray(sub.X[:500]).ravel()
log('X max=%.3f allint=%s (expect lognorm, non-integer)'%(d.max(), bool(np.allclose(d,np.round(d)))))

# ---- expressed-gene exports for NicheNet ----
vn=np.array(sub.var_names)
def frac_expr(a):
    Xc=a.X;
    return np.asarray((Xc>0).mean(axis=0)).ravel()
# receiver background (both conditions)
rec=sub[sub.obs['ct_work']==RECEIVER]
log('receiver cells:', rec.n_obs)
rf=frac_expr(rec)
bg=vn[rf>=EXPR_PROP]
pd.Series(bg).to_csv(OUT+'/export/receiver_background.txt', index=False, header=False)
# geneset of interest: program+HES1+targets expressed in receiver (>=5%)
gset_all=PROGRAM+NOTCH_TGT
gset=[g for g in gset_all if (g in set(vn)) and rf[np.where(vn==g)[0][0]]>=0.05]
pd.Series(gset).to_csv(OUT+'/export/geneset_myofib.txt', index=False, header=False)
log('geneset (expressed>=5%% in receiver):', gset)
# sender expressed (union over niche senders present)
present_niche=[c for c in NICHE if (sub.obs['ct_work']==c).sum()>=10]
sender_union=set()
sender_frac={}
for c in present_niche:
    fc=frac_expr(sub[sub.obs['ct_work']==c])
    sender_frac[c]=fc
    sender_union|=set(vn[fc>=EXPR_PROP].tolist())
pd.Series(sorted(sender_union)).to_csv(OUT+'/export/sender_expressed.txt', index=False, header=False)
log('present niche senders:', present_niche, '| union expressed genes:', len(sender_union))

# diagnostics: Notch ligand/receptor fractions per ct_work per condition
diag=[]
for cn in ['SSc','HC']:
    for ct in [RECEIVER]+present_niche:
        a=sub[(sub.obs['ct_work']==ct)&(sub.obs['condition']==cn)]
        if a.n_obs<10: continue
        fc=frac_expr(a)
        row={'condition':cn,'ct_work':ct,'n_cells':a.n_obs}
        for g in NOTCH_LIG+NOTCH_REC+['HES1']:
            row[g]=float(fc[np.where(vn==g)[0][0]]) if g in set(vn) else np.nan
        diag.append(row)
pd.DataFrame(diag).to_csv(OUT+'/export/notch_expr_fraction.csv', index=False)
log('notch expr fraction diagnostic saved')

# ---- LIANA per condition ----
summary={}
for cn in ['SSc','HC']:
    ac=sub[sub.obs['condition']==cn].copy()
    vc=ac.obs['ct_work'].value_counts(); keepct=vc[vc>=10].index.tolist()
    ac=ac[ac.obs['ct_work'].isin(keepct)].copy()
    ac.obs['ct_work']=ac.obs['ct_work'].astype('category')
    log('LIANA %s: %d cells, %d celltypes'%(cn, ac.n_obs, len(keepct)))
    li.mt.rank_aggregate(ac, groupby='ct_work', resource_name='consensus',
                         expr_prop=EXPR_PROP, use_raw=False, n_perms=N_PERMS,
                         seed=SEED, verbose=True)
    res=ac.uns['liana_res'].copy()
    res.to_csv(OUT+'/liana/liana_%s_all.csv'%cn, index=False)
    nf=res[res['ligand_complex'].isin(NOTCH_LIG) & res['receptor_complex'].isin(NOTCH_REC)].copy()
    nf.to_csv(OUT+'/liana/liana_%s_notch.csv'%cn, index=False)
    nf_rec=nf[(nf['target']==RECEIVER) & (nf['source'].isin(present_niche))].copy()
    nf_rec.to_csv(OUT+'/liana/liana_%s_notch_niche2myo.csv'%cn, index=False)
    summary[cn]={'n_cells':int(ac.n_obs),'n_celltypes':len(keepct),
                 'n_notch_pairs':int(len(nf)),'n_notch_niche2myo':int(len(nf_rec))}
    log('LIANA %s done. notch pairs=%d, niche->myo notch=%d'%(cn,len(nf),len(nf_rec)))

json.dump(summary, open(OUT+'/liana/liana_summary.json','w'), indent=2)
log('ALL LIANA DONE'); log('LIANA_OK')

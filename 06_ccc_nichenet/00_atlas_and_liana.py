#!/usr/bin/env python
"""One-load driver: atlas druggability queries (1-4) + NicheNet inputs + LIANA (SSc vs HC).
Loads full lognorm X once (memory-frugal) and serves both the quick queries and LIANA."""
import os, sys, json, time
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ.setdefault(k,'4')
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad, h5py
from scipy.stats import mannwhitneyu
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
try: from anndata.io import read_elem
except Exception: from anndata.experimental import read_elem

FULL='/data/ssc/powered/full_integrated_annotated.h5ad'
FIB='/data/ssc/powered/fibroblast_integrated.h5ad'
COND='/data/ssc/powered/sample_condition.csv'
AQ='/data/ssc/atlas_queries'; AF=AQ+'/figures'; os.makedirs(AF, exist_ok=True)
CO='/data/ssc/ccc_nichenet'; os.makedirs(CO+'/liana', exist_ok=True); os.makedirs(CO+'/export', exist_ok=True)
def log(*a): print('[atlas]',*a,flush=True)

SEED=0; CAP=8000; EXPR_PROP=0.10; N_PERMS=1000
NICHE=['Endothelial','Pericyte','SmoothMuscle','Myeloid','Lymphatic','Mast']
RECEIVER='Myofibroblast'
NOTCH_LIG=['DLL1','DLL3','DLL4','JAG1','JAG2']; NOTCH_REC=['NOTCH1','NOTCH2','NOTCH3','NOTCH4']
PROGRAM=['ACTA2','TAGLN','POSTN','COL1A1','COL1A2','COL3A1','COL11A1','COMP','CTHRC1','FN1']
NOTCH_TGT=['HES1','HEY1','HEYL','HES5','NRARP','DTX1']
DRUG_GENES=['NOTCH1','NOTCH2','NOTCH3','PSEN1','PSEN2','NCSTN','APH1A','APH1B','PSENEN']

# ---------- metadata ----------
log('reading obs (backed)...')
full_b=ad.read_h5ad(FULL, backed='r'); obs=full_b.obs.copy()
fibsub=ad.read_h5ad(FIB, backed='r').obs['fib_subtype'].astype(str)
cond=pd.read_csv(COND); s2c=dict(zip(cond['sample'].astype(str), cond['condition'].astype(str)))
ctw=obs['celltype'].astype(str).copy()
common=obs.index.intersection(fibsub.index); ctw.loc[common]=fibsub.loc[common].values
ctw[(obs['celltype'].astype(str)=='Fibroblast') & (~obs.index.isin(common))]='Fibroblast_other'
obs['ct_work']=ctw.values
obs['condition']=obs['sample'].astype(str).map(s2c).fillna('unknown')
log('ct_work:\n'+obs['ct_work'].value_counts().to_string())
log('condition:\n'+obs['condition'].value_counts().to_string())

# ---------- load full lognorm X once ----------
log('loading full lognorm X (read_elem)...'); t0=time.time()
with h5py.File(FULL,'r') as f:
    X=read_elem(f['X']); umap=read_elem(f['obsm']['X_umap'])
    var=read_elem(f['var'])
vn=np.asarray(var.index if hasattr(var,'index') else var['_index'])
gidx={g:i for i,g in enumerate(vn)}
log('X', X.shape, 'loaded %.0fs; X max=%.2f'%(time.time()-t0, X[:500].data.max()))
def group_stats(genes, gcol, sub_mask=None):
    present=list(dict.fromkeys([g for g in genes if g in gidx]))  # dedupe, keep order
    cols=[gidx[g] for g in present]
    D=np.asarray(X[:, cols].todense())
    df=pd.DataFrame(D, columns=present); df['g']=gcol.values
    if sub_mask is not None: df=df[sub_mask.values]
    mean=df.groupby('g')[present].mean()
    pct=(df[present]>0).assign(g=df['g'].values).groupby('g')[present].mean()
    return present, mean, pct

# ---------- Q1: druggability targets in myofibroblast, SSc vs HC ----------
log('=== Q1 druggability targets ===')
q1_rows=[]
for cn in ['SSc','HC']:
    m=(obs['condition']==cn)
    pres,mean,pct=group_stats(DRUG_GENES, obs['ct_work'], sub_mask=m)
    for g in pres:
        for ct in mean.index:
            q1_rows.append({'condition':cn,'gene':g,'ct_work':ct,
                            'mean_lognorm':float(mean.loc[ct,g]),'pct_expr':float(pct.loc[ct,g])})
q1=pd.DataFrame(q1_rows); q1.to_csv(AQ+'/Q1_druggability_target_expression.csv', index=False)
myo1=q1[q1.ct_work==RECEIVER].pivot_table(index='gene',columns='condition',values=['mean_lognorm','pct_expr'])
log('Q1 Myofibroblast target expression (mean | pct):\n'+myo1.round(3).to_string())

# ---------- Q3: HES1 cell-type specificity ----------
log('=== Q3 HES1 specificity ===')
pres,mean,pct=group_stats(['HES1']+NOTCH_TGT, obs['ct_work'])
q3=pd.DataFrame({'ct_work':mean.index,'HES1_mean':mean['HES1'].values,'HES1_pct':pct['HES1'].values}).sort_values('HES1_mean',ascending=False)
q3.to_csv(AQ+'/Q3_HES1_by_celltype.csv', index=False)
log('Q3 HES1 by ct_work (top):\n'+q3.head(12).round(3).to_string(index=False))

# ---------- Q2: myofibroblast proportion SSc vs HC ----------
log('=== Q2 celltype proportions ===')
keep_s=obs['condition'].isin(['SSc','HC'])
ob2=obs[keep_s].copy()
tot=ob2.groupby('sample').size(); good=tot[tot>=100].index
ob2=ob2[ob2['sample'].isin(good)]
ct_counts=pd.crosstab(ob2['sample'], ob2['ct_work'])            # samples x cts, missing=0
prop_mat=ct_counts.div(ct_counts.sum(1), axis=0)
prop=prop_mat.reset_index().melt(id_vars='sample', var_name='ct_work', value_name='prop')
smeta=ob2.drop_duplicates('sample')[['sample','condition']]
prop=prop.merge(smeta,on='sample')
prop.to_csv(AQ+'/Q2_per_sample_proportion.csv', index=False)
q2_rows=[]
for ct in prop['ct_work'].unique():
    a=prop[(prop.ct_work==ct)&(prop.condition=='SSc')]['prop']; b=prop[(prop.ct_work==ct)&(prop.condition=='HC')]['prop']
    if len(a)>=3 and len(b)>=3:
        try: p=mannwhitneyu(a,b,alternative='two-sided').pvalue
        except Exception: p=np.nan
        q2_rows.append({'ct_work':ct,'median_SSc':float(a.median()),'median_HC':float(b.median()),
                        'n_SSc':int(len(a)),'n_HC':int(len(b)),'p_MWU':float(p),
                        'log2FC':float(np.log2((a.median()+1e-4)/(b.median()+1e-4)))})
q2=pd.DataFrame(q2_rows).sort_values('log2FC',ascending=False)
q2.to_csv(AQ+'/Q2_proportion_stats.csv', index=False)
log('Q2 proportion SSc vs HC (myofib & top):\n'+q2[q2.ct_work.isin([RECEIVER,'SFRP4_proFib','Adipogenic','Endothelial','Pericyte'])].round(4).to_string(index=False))

# ---------- Q1/Q3 figures ----------
def dotplot(df_long, genes, cts, title, fn):
    df=df_long[df_long.gene.isin(genes) & df_long.ct_work.isin(cts)]
    conds=['SSc','HC']; fig,axes=plt.subplots(1,2,figsize=(1.1*len(cts)+2.5,0.45*len(genes)+1.6),sharey=True)
    vmax=df.mean_lognorm.max() or 1
    for ax,cn in zip(axes,conds):
        d=df[df.condition==cn]
        for _,r in d.iterrows():
            ax.scatter(cts.index(r.ct_work), genes.index(r.gene), s=8+r.pct_expr*260,
                       c=[r.mean_lognorm], cmap='Reds', vmin=0, vmax=vmax, edgecolors='k', linewidths=0.3)
        ax.set_xticks(range(len(cts))); ax.set_xticklabels(cts, rotation=45, ha='right', fontsize=7)
        ax.set_title(cn, fontsize=9); ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=8)
        ax.grid(True, ls=':', alpha=0.3)
    fig.suptitle(title, fontsize=10); sm=plt.cm.ScalarMappable(cmap='Reds',norm=plt.Normalize(0,vmax))
    fig.colorbar(sm, ax=axes, fraction=0.03, pad=0.02, label='mean lognorm')
    plt.savefig(fn+'.png', dpi=200, bbox_inches='tight'); plt.savefig(fn+'.pdf', bbox_inches='tight'); plt.close()
cts_show=[RECEIVER,'SFRP4_proFib','Adipogenic','SFRP2_DPP4','Endothelial','Pericyte','SmoothMuscle','Myeloid','Keratinocyte','Tcell']
cts_show=[c for c in cts_show if c in set(q1.ct_work)]
dotplot(q1, DRUG_GENES, cts_show, 'Nirogacestat targets (gamma-secretase / Notch receptors)', AF+'/Q1_druggability_dotplot')
log('Q1 dotplot saved')
# Q3 lollipop
q3s=q3.set_index('ct_work').loc[[c for c in cts_show if c in q3.ct_work.values]].reset_index() if len(cts_show) else q3
q3p=q3.sort_values('HES1_mean').tail(14)
fig,ax=plt.subplots(figsize=(4.4,5)); y=np.arange(len(q3p))
cols=['#C44E52' if c==RECEIVER else '#888' for c in q3p.ct_work]
ax.hlines(y,0,q3p.HES1_mean,color=cols,lw=1.8); ax.scatter(q3p.HES1_mean,y,color=cols,s=45,zorder=3)
ax.set_yticks(y); ax.set_yticklabels(q3p.ct_work,fontsize=8); ax.set_xlabel('HES1 mean lognorm')
ax.set_title('HES1 expression by cell type',fontsize=10); ax.grid(True,axis='x',ls=':',alpha=0.4)
plt.tight_layout(); plt.savefig(AF+'/Q3_HES1_lollipop.png',dpi=200); plt.savefig(AF+'/Q3_HES1_lollipop.pdf'); plt.close()
# Q2 dumbbell
q2p=q2[q2.n_SSc>=3].sort_values('median_SSc').tail(14)
fig,ax=plt.subplots(figsize=(5,5)); y=np.arange(len(q2p))
ax.hlines(y,q2p.median_HC,q2p.median_SSc,color='#bbb',lw=2)
ax.scatter(q2p.median_HC,y,color='#4C72B0',s=40,label='HC',zorder=3); ax.scatter(q2p.median_SSc,y,color='#C44E52',s=40,label='SSc',zorder=3)
for i,(_,r) in enumerate(q2p.iterrows()):
    if r.p_MWU<0.05: ax.text(max(r.median_SSc,r.median_HC)+0.005,i,'*',fontsize=12,va='center')
ax.set_yticks(y); ax.set_yticklabels(q2p.ct_work,fontsize=8); ax.set_xlabel('per-sample proportion'); ax.legend(fontsize=8)
ax.set_title('Cell-type proportion SSc vs HC\n(* MWU p<0.05)',fontsize=10); ax.grid(True,axis='x',ls=':',alpha=0.4)
plt.tight_layout(); plt.savefig(AF+'/Q2_proportion_dumbbell.png',dpi=200); plt.savefig(AF+'/Q2_proportion_dumbbell.pdf'); plt.close()
log('Q2/Q3 figures saved')

# ---------- Q4: UMAP overview ----------
log('=== Q4 UMAP overview ===')
hes1=np.asarray(X[:, gidx['HES1']].todense()).ravel() if 'HES1' in gidx else np.zeros(X.shape[0])
pd_umap=pd.DataFrame({'umap1':umap[:,0],'umap2':umap[:,1],'celltype':obs['celltype'].values,
                      'ct_work':obs['ct_work'].values,'condition':obs['condition'].values,
                      'HES1':hes1,'myo_score':obs['myo'].values})
pd_umap.to_csv(AQ+'/Q4_umap_plotdata.csv.gz', index=False, compression='gzip')
rng=np.random.default_rng(0); n=pd_umap.shape[0]; sel=rng.choice(n, min(n,140000), replace=False)
u=pd_umap.iloc[sel]
fig,axes=plt.subplots(2,2,figsize=(13,11))
# a celltype
cts=sorted(u.celltype.unique()); cmap=plt.cm.tab20(np.linspace(0,1,len(cts)))
for c,col in zip(cts,cmap):
    d=u[u.celltype==c]; axes[0,0].scatter(d.umap1,d.umap2,s=1.2,color=col,label=c,rasterized=True)
axes[0,0].legend(markerscale=6,fontsize=6,ncol=2,loc='best'); axes[0,0].set_title('a  Cell type')
# b condition
for c,col in zip(['SSc','HC','unknown'],['#C44E52','#4C72B0','#cccccc']):
    d=u[u.condition==c]; axes[0,1].scatter(d.umap1,d.umap2,s=1.2,color=col,label=c,rasterized=True)
axes[0,1].legend(markerscale=6,fontsize=8); axes[0,1].set_title('b  Condition')
# c HES1
sc2=axes[1,0].scatter(u.umap1,u.umap2,s=1.2,c=u.HES1,cmap='magma',vmax=np.quantile(u.HES1,0.99),rasterized=True)
fig.colorbar(sc2,ax=axes[1,0],fraction=0.04); axes[1,0].set_title('c  HES1 (lognorm)')
# d myo score
sc3=axes[1,1].scatter(u.umap1,u.umap2,s=1.2,c=u.myo_score,cmap='viridis',vmax=np.quantile(u.myo_score,0.99),rasterized=True)
fig.colorbar(sc3,ax=axes[1,1],fraction=0.04); axes[1,1].set_title('d  Myofibroblast program score')
for ax in axes.ravel(): ax.set_xticks([]); ax.set_yticks([])
plt.tight_layout(); plt.savefig(AF+'/Q4_umap_overview.png',dpi=170); plt.savefig(AF+'/Q4_umap_overview.pdf'); plt.close()
log('Q4 UMAP saved')

# ---------- key numbers for coordinator ----------
myo_prop=q2[q2.ct_work==RECEIVER]
KEY={'Q1_myofib_targets':{g:{'SSc_mean':float(q1[(q1.gene==g)&(q1.ct_work==RECEIVER)&(q1.condition=='SSc')]['mean_lognorm'].values[0]) if len(q1[(q1.gene==g)&(q1.ct_work==RECEIVER)&(q1.condition=='SSc')]) else None,
                              'SSc_pct':float(q1[(q1.gene==g)&(q1.ct_work==RECEIVER)&(q1.condition=='SSc')]['pct_expr'].values[0]) if len(q1[(q1.gene==g)&(q1.ct_work==RECEIVER)&(q1.condition=='SSc')]) else None} for g in DRUG_GENES},
     'Q2_myofib_proportion':(myo_prop.round(4).to_dict('records')[0] if len(myo_prop) else None),
     'Q3_HES1_top3':q3.head(3)[['ct_work','HES1_mean','HES1_pct']].round(3).to_dict('records'),
     'Q3_HES1_myofib_rank':int(q3.reset_index(drop=True).index[q3.reset_index(drop=True).ct_work==RECEIVER][0])+1 if RECEIVER in set(q3.ct_work) else None}
json.dump(KEY, open(AQ+'/KEY_NUMBERS.json','w'), indent=2, default=str)
log('KEY_NUMBERS:\n'+json.dumps(KEY, indent=2, default=str))
log('QUERIES_DONE')

# ================= LIANA (subset from in-memory X) =================
log('=== building LIANA subset ===')
rng=np.random.default_rng(SEED); KEEP_FULL=set([RECEIVER])|set(NICHE); keep=[]
for cn in ['SSc','HC']:
    so=obs[obs['condition']==cn]
    for ct,g in so.groupby('ct_work'):
        idx=g.index.to_numpy()
        if len(idx)==0: continue
        keep.append(idx if (ct in KEEP_FULL or len(idx)<=CAP) else rng.choice(idx,CAP,replace=False))
keep=np.concatenate(keep); pos=obs.index.get_indexer(keep)
sub=ad.AnnData(X=X[pos].copy(), obs=obs.loc[keep,['ct_work','condition','sample']].copy(),
               var=pd.DataFrame(index=vn))
log('LIANA subset', sub.shape)
# NicheNet exports
def frac(a): return np.asarray((a.X>0).mean(axis=0)).ravel()
rec=sub[sub.obs.ct_work==RECEIVER]; rf=frac(rec)
bg=vn[rf>=EXPR_PROP]; pd.Series(bg).to_csv(CO+'/export/receiver_background.txt',index=False,header=False)
gset=[g for g in PROGRAM+NOTCH_TGT if g in gidx and rf[list(vn).index(g)]>=0.05]
pd.Series(gset).to_csv(CO+'/export/geneset_myofib.txt',index=False,header=False)
present_niche=[c for c in NICHE if (sub.obs.ct_work==c).sum()>=10]
su=set()
for c in present_niche: su|=set(vn[frac(sub[sub.obs.ct_work==c])>=EXPR_PROP].tolist())
pd.Series(sorted(su)).to_csv(CO+'/export/sender_expressed.txt',index=False,header=False)
log('exported nichenet inputs: bg=%d geneset=%s sender=%d'%(len(bg),gset,len(su)))
# notch expr fraction diagnostic
diag=[]
for cn in ['SSc','HC']:
    for ct in [RECEIVER]+present_niche:
        a=sub[(sub.obs.ct_work==ct)&(sub.obs.condition==cn)]
        if a.n_obs<10: continue
        fc=frac(a); row={'condition':cn,'ct_work':ct,'n_cells':int(a.n_obs)}
        for g in NOTCH_LIG+NOTCH_REC+['HES1']: row[g]=float(fc[list(vn).index(g)]) if g in gidx else np.nan
        diag.append(row)
pd.DataFrame(diag).to_csv(CO+'/export/notch_expr_fraction.csv',index=False)
# free full X
del X; import gc; gc.collect()
import liana as li
summ={}
for cn in ['SSc','HC']:
    ac=sub[sub.obs.condition==cn].copy()
    vc=ac.obs.ct_work.value_counts(); kc=vc[vc>=10].index.tolist()
    ac=ac[ac.obs.ct_work.isin(kc)].copy(); ac.obs['ct_work']=ac.obs.ct_work.astype('category')
    log('LIANA %s: %d cells %d cts'%(cn,ac.n_obs,len(kc)))
    li.mt.rank_aggregate(ac, groupby='ct_work', resource_name='consensus', expr_prop=EXPR_PROP,
                         use_raw=False, n_perms=N_PERMS, seed=SEED, verbose=True)
    res=ac.uns['liana_res'].copy(); res.to_csv(CO+'/liana/liana_%s_all.csv'%cn,index=False)
    nf=res[res.ligand_complex.isin(NOTCH_LIG)&res.receptor_complex.isin(NOTCH_REC)].copy()
    nf.to_csv(CO+'/liana/liana_%s_notch.csv'%cn,index=False)
    nfr=nf[(nf.target==RECEIVER)&(nf.source.isin(present_niche))].copy(); nfr.to_csv(CO+'/liana/liana_%s_notch_niche2myo.csv'%cn,index=False)
    summ[cn]={'n_cells':int(ac.n_obs),'n_notch_pairs':int(len(nf)),'n_notch_niche2myo':int(len(nfr))}
    log('LIANA %s notch=%d niche2myo=%d'%(cn,len(nf),len(nfr)))
json.dump(summ, open(CO+'/liana/liana_summary.json','w'), indent=2)
log('LIANA_OK')

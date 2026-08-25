#!/usr/bin/env python
"""Combine LIANA + NicheNet, SSc-vs-HC comparison, figures (dot/lollipop/heatmap)."""
import os, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT='/data/ssc/ccc_nichenet'; FIG=OUT+'/figures'; os.makedirs(FIG, exist_ok=True)
NOTCH_LIG=['DLL1','DLL3','DLL4','JAG1','JAG2']; NICHE=['Endothelial','Pericyte','SmoothMuscle','Myeloid','Lymphatic','Mast']
def rd(p): return pd.read_csv(p) if os.path.exists(p) else pd.DataFrame()

# ---- 1. LIANA niche->myo Notch, SSc vs HC ----
liana={}
for cn in ['SSc','HC']:
    df=rd(f'{OUT}/liana/liana_{cn}_notch_niche2myo.csv')
    if not df.empty: df['condition']=cn
    liana[cn]=df
lia=pd.concat([liana['SSc'],liana['HC']], ignore_index=True) if any(not v.empty for v in liana.values()) else pd.DataFrame()
if not lia.empty:
    lia['pair']=lia['source']+': '+lia['ligand_complex']+'->'+lia['receptor_complex']
    lia['mag_strength']=-np.log10(lia['magnitude_rank'].clip(lower=1e-12))
    lia['spec_strength']=-np.log10(lia['specificity_rank'].clip(lower=1e-12))
    lia.to_csv(f'{OUT}/liana_notch_niche2myo_SSc_vs_HC.csv', index=False)
    # dotplot: rows=pair, cols=condition, size=spec_strength, color=mag_strength
    pairs=sorted(lia['pair'].unique())
    conds=['SSc','HC']; ymap={p:i for i,p in enumerate(pairs)}; xmap={c:i for i,c in enumerate(conds)}
    fig,ax=plt.subplots(figsize=(4.2, max(2.5,0.42*len(pairs)+1)))
    sc=None
    for _,r in lia.iterrows():
        s=sc=ax.scatter(xmap[r['condition']], ymap[r['pair']],
            s=30+r['spec_strength']*55, c=[r['mag_strength']], cmap='viridis',
            vmin=lia['mag_strength'].min(), vmax=lia['mag_strength'].max(),
            edgecolors='k', linewidths=0.4, zorder=3)
    ax.set_xticks([0,1]); ax.set_xticklabels(conds); ax.set_yticks(range(len(pairs))); ax.set_yticklabels(pairs, fontsize=8)
    ax.set_xlim(-0.5,1.5); ax.set_ylim(-0.7,len(pairs)-0.3); ax.set_title('Niche -> Myofibroblast\nNotch ligand-receptor (LIANA)', fontsize=9)
    ax.grid(True, axis='y', ls=':', alpha=0.4)
    cb=fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.02); cb.set_label('magnitude strength\n-log10(rank)', fontsize=7)
    ax.text(1.02,-0.09,'dot size = specificity', transform=ax.transAxes, fontsize=6.5, ha='right')
    plt.tight_layout(); plt.savefig(f'{FIG}/Fig_liana_notch_niche2myo.png', dpi=200); plt.savefig(f'{FIG}/Fig_liana_notch_niche2myo.pdf'); plt.close()
    print('[fig] LIANA notch dotplot saved,', len(pairs),'pairs')
else:
    print('[fig] no LIANA niche->myo Notch pairs')

# ---- 2. Notch ligand expression fraction per niche sender, SSc vs HC (dumbbell) ----
ef=rd(f'{OUT}/export/notch_expr_fraction.csv')
if not ef.empty:
    lig_present=[g for g in NOTCH_LIG if g in ef.columns and ef[g].fillna(0).max()>0.02]
    senders=[s for s in NICHE if s in ef['ct_work'].unique()]
    rows=[]
    for s in senders:
        for g in lig_present:
            ssc=ef[(ef.ct_work==s)&(ef.condition=='SSc')][g]
            hc =ef[(ef.ct_work==s)&(ef.condition=='HC')][g]
            rows.append({'label':f'{s}: {g}','sender':s,'ligand':g,
                         'SSc':float(ssc.iloc[0]) if len(ssc) else np.nan,
                         'HC':float(hc.iloc[0]) if len(hc) else np.nan})
    dd=pd.DataFrame(rows).dropna(subset=['SSc','HC'], how='all')
    dd=dd[(dd[['SSc','HC']].max(axis=1)>0.02)].reset_index(drop=True)
    dd.to_csv(f'{OUT}/notch_ligand_fraction_SSc_vs_HC.csv', index=False)
    if not dd.empty:
        dd=dd.sort_values('SSc')
        fig,ax=plt.subplots(figsize=(5, max(2.5,0.4*len(dd)+1)))
        y=np.arange(len(dd))
        ax.hlines(y, dd['HC'], dd['SSc'], color='#bbbbbb', lw=2, zorder=1)
        ax.scatter(dd['HC'], y, color='#4C72B0', s=45, label='HC', zorder=3)
        ax.scatter(dd['SSc'], y, color='#C44E52', s=45, label='SSc', zorder=3)
        ax.set_yticks(y); ax.set_yticklabels(dd['label'], fontsize=8)
        ax.set_xlabel('fraction of cells expressing ligand'); ax.legend(fontsize=8, loc='lower right')
        ax.set_title('Notch-ligand expression in niche senders\nSSc vs HC', fontsize=9); ax.grid(True, axis='x', ls=':', alpha=0.4)
        plt.tight_layout(); plt.savefig(f'{FIG}/Fig_notch_ligand_fraction_dumbbell.png', dpi=200); plt.savefig(f'{FIG}/Fig_notch_ligand_fraction_dumbbell.pdf'); plt.close()
        print('[fig] notch ligand dumbbell saved,', len(dd),'sender-ligand')
else:
    print('[fig] no notch expr fraction table')

# ---- 3. NicheNet ligand activity lollipop + notch highlight ----
la=rd(f'{OUT}/nichenet/ligand_activities.csv')
if not la.empty:
    la=la.sort_values('aupr_corrected', ascending=False).reset_index(drop=True)
    top=la.head(25).iloc[::-1].reset_index(drop=True)
    isnotch=top['test_ligand'].isin(NOTCH_LIG)
    fig,ax=plt.subplots(figsize=(4.6, 6.2))
    y=np.arange(len(top)); cols=['#C44E52' if n else '#888888' for n in isnotch]
    ax.hlines(y, 0, top['aupr_corrected'], color=cols, lw=1.6)
    ax.scatter(top['aupr_corrected'], y, color=cols, s=40, zorder=3)
    labels=[f'{l} *' if n else l for l,n in zip(top['test_ligand'], isnotch)]
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel('ligand activity (aupr_corrected)'); ax.set_title('NicheNet ligand activity for\nmyofibroblast/HES1 program\n(* = Notch ligand)', fontsize=9)
    ax.grid(True, axis='x', ls=':', alpha=0.4)
    plt.tight_layout(); plt.savefig(f'{FIG}/Fig_nichenet_ligand_activity.png', dpi=200); plt.savefig(f'{FIG}/Fig_nichenet_ligand_activity.pdf'); plt.close()
    print('[fig] nichenet lollipop saved')

# ---- 4. Notch ligand->target heatmap ----
lt=rd(f'{OUT}/nichenet/notch_ligand_target_links.csv')
if not lt.empty:
    piv=lt.pivot_table(index='ligand', columns='target', values='weight', aggfunc='max').fillna(0)
    fig,ax=plt.subplots(figsize=(max(4,0.5*piv.shape[1]+2), max(2,0.5*piv.shape[0]+1.5)))
    im=ax.imshow(piv.values, cmap='magma', aspect='auto')
    ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index, fontsize=9)
    ax.set_title('Notch-ligand -> myofibroblast target\nregulatory potential (NicheNet)', fontsize=9)
    cb=fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02); cb.set_label('prior regulatory potential', fontsize=7)
    plt.tight_layout(); plt.savefig(f'{FIG}/Fig_notch_ligand_target_heatmap.png', dpi=200); plt.savefig(f'{FIG}/Fig_notch_ligand_target_heatmap.pdf'); plt.close()
    print('[fig] notch ligand-target heatmap saved,', piv.shape)

# ---- verdict / combined summary ----
V={}
if not lia.empty:
    ssc_m=lia[lia.condition=='SSc']['mag_strength'].mean(); hc_m=lia[lia.condition=='HC']['mag_strength'].mean()
    top_send=lia.groupby('source')['mag_strength'].mean().sort_values(ascending=False)
    V['liana_notch_niche2myo_pairs_SSc']=int((lia.condition=='SSc').sum())
    V['liana_notch_niche2myo_pairs_HC']=int((lia.condition=='HC').sum())
    V['liana_dominant_notch_sender']=top_send.index[0] if len(top_send) else None
    V['liana_mag_strength_SSc_vs_HC']=[round(float(ssc_m),3), round(float(hc_m),3)]
st=OUT+'/nichenet/nichenet_status.txt'
if os.path.exists(st): V['nichenet_status']=open(st).read().strip().splitlines()
json.dump(V, open(OUT+'/RESULT_SUMMARY.json','w'), indent=2)
print('[report] RESULT_SUMMARY:', json.dumps(V, indent=2))
print('[report] DONE')

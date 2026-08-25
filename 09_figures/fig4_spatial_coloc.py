"""Fig 5 — HES1 co-localises with the myofibroblast programme in situ.
(a,b) tissue maps of HES1 regulon activity and the myofibroblast programme on a representative SSc
Visium section (HES1 regulon activity = decoupler ULM + CollecTRI, the signal used for the spatial
statistics; raw HES1 gene is dropout-prone). (c) bivariate Moran's I by platform (Visium primary,
Xenium is a targeted panel). (d) control ordering across Visium slices (SMAD3 > FOSB > HES1 > MEF2C).
Verified. NO bar charts. Reuses figstyle. Vector PDF.
"""
import os, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, scanpy as sc, decoupler as dc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import figstyle as fs

H5 = r'C:\Users\fsy\Desktop\SSc_virtualKO_MR\server_archive\tier2_large\spatial_extract\raw\spatial\zenodo_14577696\visium_all.h5ad'
NET = r'C:\Users\fsy\Desktop\SSc_virtualKO_MR\02_analysis_modules\rigor_fixes\collectri_net.csv'
DLS = r'C:\Users\fsy\Desktop\SSc_virtualKO_MR\04_manuscript\plot_data_local\powered\spatial'
OUT = r'C:\Users\fsy\Desktop\SSc_virtualKO_MR\03_results_figures\figures\fig4_spatial_coloc'
MYO = ['ACTA2', 'TAGLN', 'COL1A1', 'COL1A2', 'POSTN', 'FN1']; SLICE = 'SSc-HL35'

# ---- tissue: HES1 regulon activity + myofib on the representative SSc section ----
a = sc.read_h5ad(H5)
if float(a.X.max()) > 30: sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
sc.tl.score_genes(a, [g for g in MYO if g in a.var_names], score_name='myofib')
sub = a[a.obs['sample'] == SLICE].copy()
dc.mt.ulm(sub, pd.read_csv(NET), tmin=5, verbose=False)
xy = sub.obsm['spatial']; hes_act = sub.obsm['score_ulm']['HES1'].values; myo = sub.obs['myofib'].values

fig = plt.figure(figsize=(13.4, 8.4))
gs = GridSpec(2, 2, height_ratios=[1.05, 0.85], wspace=0.28, hspace=0.42, figure=fig)

def tmap(ax, vals, title, cmap):
    o = np.argsort(vals)
    s = ax.scatter(xy[o, 0], -xy[o, 1], c=np.asarray(vals)[o], s=10, cmap=cmap, linewidths=0, rasterized=True)
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9.6)
    for sp in ('top', 'right', 'left', 'bottom'): ax.spines[sp].set_visible(False)
    cb = plt.colorbar(s, ax=ax, fraction=0.046, pad=0.02); cb.ax.tick_params(labelsize=6.5)

axA = fig.add_subplot(gs[0, 0]); tmap(axA, hes_act, f'{SLICE} — HES1 regulon activity', 'magma'); fs.panel_label(axA, 'a')
axB = fig.add_subplot(gs[0, 1]); tmap(axB, myo, f'{SLICE} — myofibroblast programme', 'viridis'); fs.panel_label(axB, 'b')

# ---- quant: platform + control ordering (verified pooled results) ----
d = pd.read_csv(os.path.join(DLS, 'pooled_coupling.csv'))
d['platform'] = d['slice'].apply(lambda s: 'Xenium' if 'Xenium' in s else 'Visium')
prog = d[d.substate == 'myofib_program']

# (c) HES1 moranBV by platform
axC = fig.add_subplot(gs[1, 0]); h = prog[prog.regulon == 'HES1']
rng = np.random.RandomState(0); pcol = {'Visium': fs.LEAD_BLUE, 'Xenium': '#BDBDBD'}
for i, pf in enumerate(['Visium', 'Xenium']):
    g = h[h.platform == pf]; x = i + (rng.rand(len(g)) - 0.5) * 0.28
    axC.scatter(x, g.moranBV, s=34, color=pcol[pf], alpha=0.85, edgecolor='white', linewidth=0.4, zorder=3)
    m = g.moranBV.mean(); axC.plot([i - 0.26, i + 0.26], [m, m], color='#333', lw=2.2, zorder=4)
    axC.text(i + 0.32, m, f"mean {m:.2f}\n{100*(g.perm_p_coord<0.05).mean():.0f}% sig", ha='left', va='center', fontsize=7.4, fontweight='bold', color='#333')
axC.axhline(0, color='#888', lw=0.7); axC.set_xticks([0, 1]); axC.set_xticklabels([f"Visium\n(n={h[h.platform=='Visium'].shape[0]})", f"Xenium\n(n={h[h.platform=='Xenium'].shape[0]})"], fontsize=7.8)
axC.set_ylabel("HES1 bivariate Moran's I", fontsize=8.5); axC.set_title('Spatial co-localisation by platform', fontsize=9.6); axC.set_xlim(-0.6, 1.95)
fs.despine(axC); fs.panel_label(axC, 'c')

# (d) control ordering (Visium)
axD = fig.add_subplot(gs[1, 1]); vis = prog[prog.platform == 'Visium']
order = ['SMAD3', 'FOSB', 'HES1', 'MEF2C']; role = {'SMAD3': 'positive control', 'FOSB': 'co-lead', 'HES1': 'lead', 'MEF2C': 'negative control'}
cmap = {'SMAD3': fs.CTRL_COLOR, 'FOSB': fs.LEAD_RED, 'HES1': fs.LEAD_BLUE, 'MEF2C': '#BDBDBD'}
rows = [(r, vis[vis.regulon == r].moranBV.mean(), 100*(vis[vis.regulon == r].perm_p_coord < 0.05).mean()) for r in order][::-1]
y = np.arange(len(rows))
for yy, (reg, m, fsig) in zip(y, rows):
    axD.hlines(yy, 0, m, color=cmap[reg], lw=2.8); axD.scatter(m, yy, color=cmap[reg], s=90, zorder=3, edgecolor='white', linewidth=0.8)
    axD.text(m + (0.02 if m >= 0 else -0.02), yy, f"{m:.2f} ({fsig:.0f}% sig)", ha='left' if m >= 0 else 'right', va='center', fontsize=7.6, fontweight='bold', color='#333')
axD.axvline(0, color='#888', lw=0.7); axD.set_yticks(y); axD.set_yticklabels([f"{r}\n({role[r]})" for r, _, _ in rows], fontsize=8)
axD.set_xlabel("mean bivariate Moran's I (Visium, 14 slices)", fontsize=8.5); axD.set_title('In-situ control ordering', fontsize=9.6); axD.set_xlim(-0.2, 0.85)
fs.despine(axD); fs.panel_label(axD, 'd')

fig.suptitle('HES1 co-localises with the myofibroblast programme in situ (SSc dermal Visium)', fontsize=11, fontweight='bold', y=1.0)
fs.save_fig(fig, OUT)
print('Visium HES1 moranBV', round(h[h.platform=='Visium'].moranBV.mean(),3), '| tissue slice', SLICE)

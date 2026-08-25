"""Fig 5 (advanced) — in-situ tissue maps: HES1 REGULON ACTIVITY (not raw expression, which is
dropout-prone for a low-abundance TF) co-localises with the myofibroblast programme on a
representative SSc dermal Visium section. HES1 regulon activity = decoupler ULM + CollecTRI,
exactly the signal used for the bivariate Moran's I result. Controls SMAD3 (positive) / MEF2C
(negative). Verified. Vector PDF.
"""
import os, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, scanpy as sc, decoupler as dc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import figstyle as fs

H5 = r'C:\Users\fsy\Desktop\SSc_virtualKO_MR\server_archive\tier2_large\spatial_extract\raw\spatial\zenodo_14577696\visium_all.h5ad'
NET = r'C:\Users\fsy\Desktop\SSc_virtualKO_MR\02_analysis_modules\rigor_fixes\collectri_net.csv'
OUT = r'C:\Users\fsy\Desktop\SSc_virtualKO_MR\03_results_figures\figures\Fig5b_spatial_tissue'
MYO = ['ACTA2', 'TAGLN', 'COL1A1', 'COL1A2', 'POSTN', 'FN1']
SLICE = 'SSc-HL35'

a = sc.read_h5ad(H5)
if float(a.X.max()) > 30:
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
sc.tl.score_genes(a, [g for g in MYO if g in a.var_names], score_name='myofib')
net = pd.read_csv(NET)
sub = a[a.obs['sample'] == SLICE].copy()
dc.mt.ulm(sub, net, tmin=5, verbose=False)
act = sub.obsm['score_ulm']
xy = sub.obsm['spatial']; myo = sub.obs['myofib'].values

def lag(v, k=18):
    nn = NearestNeighbors(n_neighbors=min(k, len(xy))).fit(xy); _, idx = nn.kneighbors(xy)
    return np.asarray(v)[idx].mean(1)

def tmap(ax, vals, title, cmap):
    o = np.argsort(vals)
    s = ax.scatter(xy[o, 0], -xy[o, 1], c=np.asarray(vals)[o], s=9, cmap=cmap, linewidths=0, rasterized=True)
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9.6)
    for sp in ('top', 'right', 'left', 'bottom'): ax.spines[sp].set_visible(False)
    cb = plt.colorbar(s, ax=ax, fraction=0.046, pad=0.02); cb.ax.tick_params(labelsize=6.5)

fig = plt.figure(figsize=(13.2, 4.6))
gs = GridSpec(1, 3, width_ratios=[1, 1, 0.85], wspace=0.30, figure=fig)

axA = fig.add_subplot(gs[0, 0]); tmap(axA, act['HES1'].values, f'{SLICE} — HES1 regulon activity', 'magma'); fs.panel_label(axA, 'a')
axB = fig.add_subplot(gs[0, 1]); tmap(axB, myo, f'{SLICE} — myofibroblast programme', 'viridis'); fs.panel_label(axB, 'b')

# co-localisation lollipop: spatial-lag rho vs myofib for HES1 / SMAD3(pos) / MEF2C(neg)
axC = fig.add_subplot(gs[0, 2])
order = ['SMAD3', 'HES1', 'MEF2C']; role = {'SMAD3': 'positive control', 'HES1': 'lead', 'MEF2C': 'negative control'}
cmap = {'HES1': fs.LEAD_BLUE, 'SMAD3': fs.CTRL_COLOR, 'MEF2C': '#BDBDBD'}
myo_s = lag(myo)
rows = [(tf, spearmanr(lag(act[tf].values), myo_s)[0]) for tf in order if tf in act.columns][::-1]
y = np.arange(len(rows))
for yy, (tf, r) in zip(y, rows):
    axC.hlines(yy, 0, r, color=cmap[tf], lw=2.8); axC.scatter(r, yy, color=cmap[tf], s=95, zorder=3, edgecolor='white', linewidth=0.8)
    axC.text(r + (0.02 if r >= 0 else -0.02), yy, f'{r:.2f}', va='center', ha='left' if r >= 0 else 'right', fontsize=8, fontweight='bold', color='#333')
axC.axvline(0, color='#888', lw=0.7)
axC.set_yticks(y); axC.set_yticklabels([f'{tf}\n({role[tf]})' for tf, _ in rows], fontsize=8)
axC.set_xlabel('spatial co-localisation with\nmyofibroblast programme (spatial-lag ' + r'$\rho$)', fontsize=8)
axC.set_title('In-situ control ordering', fontsize=9.6); axC.set_xlim(-0.55, 1.05)
fs.despine(axC); fs.panel_label(axC, 'c')

fig.suptitle('In-situ co-localisation of HES1 regulon activity with the myofibroblast programme (SSc Visium)',
             fontsize=11, fontweight='bold', y=1.0)
fs.save_fig(fig, OUT)
print('slice', SLICE, '| co-localisation (spatial-lag rho vs myofib):', {tf: round(r, 3) for tf, r in rows})

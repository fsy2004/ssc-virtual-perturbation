"""Figure 2: CellOracle in silico knockout scores across GRN priors and
perturbation vector fields for SMAD3 and HES1.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DL = os.path.join(PROJECT, '04_manuscript', 'plot_data_local', 'powered')
VF = os.path.join(PROJECT, '03_results_figures', 'figures', 'fig2_vectorfield.png')
OUT = os.path.join(PROJECT, '03_results_figures', 'figures', 'Fig2_engine')

skin = pd.read_csv(os.path.join(DL, 'KO_ranking_skinatac.csv'))
prom = pd.read_csv(os.path.join(DL, 'KO_ranking.csv'))
LEAD = 'HES1'

fig = plt.figure(figsize=(12.6, 9.4))
gs = GridSpec(2, 2, height_ratios=[1.0, 1.02], width_ratios=[1.0, 1.05], hspace=0.30, wspace=0.34, figure=fig)

# ===== A: knockout perturbation scores =====
axA = fig.add_subplot(gs[0, 0])
show = ['MYC', 'SNAI2', 'SRF', 'STAT3', 'SMAD3', 'TWIST1', 'HES1']
d = skin[skin.gene.isin(show)].set_index('gene').loc[show].reset_index().iloc[::-1]
cols = [fs.LEAD_BLUE if g == LEAD else fs.CTRL_COLOR for g in d.gene]
y = np.arange(len(d))
axA.hlines(y, 0, d.ps_score, color=cols, lw=2.4, alpha=0.9)
axA.scatter(d.ps_score, y, color=cols, s=64, zorder=3, edgecolor='white', linewidth=0.7)
axA.set_yticks(y); axA.set_yticklabels(d.gene, fontsize=9)
for tk, g in zip(axA.get_yticklabels(), d.gene):
    if g == LEAD: tk.set_fontweight('bold'); tk.set_color(fs.LEAD_BLUE)
    elif g == 'SMAD3': tk.set_color(fs.CTRL_COLOR)
for yy, (_, r) in zip(y, d.iterrows()):
    axA.text(r.ps_score - 1.0, yy, fs.star(r.ps_qbh), ha='right', va='center', fontsize=9, fontweight='bold', color='#333')
axA.axvline(0, color='#888', lw=0.6)
axA.set_xlabel('KO perturbation score (skin scATAC GRN)')
axA.set_title('Knockout perturbation scores across candidate TFs', fontsize=9.6)
axA.invert_xaxis()
axA.legend(handles=[Line2D([0],[0], marker='o', color='w', markerfacecolor=fs.CTRL_COLOR, ms=8, label='canonical driver (control)'),
                    Line2D([0],[0], marker='o', color='w', markerfacecolor=fs.LEAD_BLUE, ms=8, label='HES1')],
           loc='lower left', fontsize=7.2)
fs.despine(axA); fs.panel_label(axA, 'a')

# ===== B: skin vs promoter network priors =====
axB = fig.add_subplot(gs[0, 1])
m = skin[['gene', 'ps_score']].merge(prom[['gene', 'ps_score']], on='gene', suffixes=('_skin', '_prom'))
r = np.corrcoef(m.ps_score_prom, m.ps_score_skin)[0, 1]
axB.scatter(m.ps_score_prom, m.ps_score_skin, s=34, color='#BDBDBD', edgecolor='#888', linewidth=0.4, zorder=2)
lim = [min(m.ps_score_prom.min(), m.ps_score_skin.min()) - 3, 2]
axB.plot(lim, lim, ls='--', color='#999', lw=0.9, zorder=1)
for g, c in {'HES1': fs.LEAD_BLUE, 'SMAD3': fs.CTRL_COLOR, 'MEF2C': '#BDBDBD'}.items():
    row = m[m.gene == g]
    if len(row):
        axB.scatter(row.ps_score_prom, row.ps_score_skin, s=95, color=c, edgecolor='black', linewidth=0.8, zorder=4)
        axB.annotate(g, (row.ps_score_prom.iloc[0], row.ps_score_skin.iloc[0]), xytext=(6, -2),
                     textcoords='offset points', fontsize=8.5, fontweight='bold', color=c if g != 'MEF2C' else '#555')
axB.set_xlabel('KO score — generic promoter prior')
axB.set_ylabel('KO score — skin scATAC prior')
axB.set_title(f'Perturbation scores across network priors (r = {r:.2f})', fontsize=9.6)
axB.annotate('HES1 significant in both analyses\n(q=1.1e-10 skin, q=8.7e-7 promoter)', xy=(0.03, 0.06),
             xycoords='axes fraction', fontsize=7.2, style='italic', color='#333', va='bottom')
fs.despine(axB); fs.panel_label(axB, 'b')

# ===== C/D: perturbation vector field (rendered CellOracle output) =====
axV = fig.add_subplot(gs[1, :])
if os.path.exists(VF):
    axV.imshow(mpimg.imread(VF)); axV.axis('off')
    axV.set_title('In silico knockout perturbation fields — SMAD3 (positive control) and HES1',
                  fontsize=9.6, pad=6)
else:
    axV.text(0.5, 0.5, 'vector field PNG missing', ha='center'); axV.axis('off')
fs.panel_label(axV, 'c')

fig.suptitle('Engine 1 — single-cell in silico knockout identifies HES1 in the myofibroblast state',
             fontsize=11.5, fontweight='bold', y=0.995)
fs.save_fig(fig, OUT)
print('prior robustness r =', round(r, 3))

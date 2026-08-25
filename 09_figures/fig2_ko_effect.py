"""Fig 2 (data panels) — In-silico knockout nominates HES1, robustly across GRN priors.
Honest, verified (plot_data_local powered, 2026-07-10). NO bar charts. Reuses figstyle.
(Vector-field quiver panels are rendered server-side from the CellOracle object; assembled later.)

A  Engine recovers canonical myofibroblast drivers (positive control) and nominates the
   druggable Notch effector HES1: perturbation score on the skin-specific scATAC GRN.
B  Robustness of the KO effect to the GRN prior: perturbation score on a generic promoter
   prior vs the skin-specific scATAC prior (all 43 panel TFs). HES1 is stable and significant
   on both — the nomination is not an artefact of the prior.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

DL = r'C:\Users\fsy\Desktop\SSc_virtualKO_MR\04_manuscript\plot_data_local\powered'
OUT = r'C:\Users\fsy\Desktop\SSc_virtualKO_MR\03_results_figures\figures\fig2_ko_effect'

skin = pd.read_csv(os.path.join(DL, 'KO_ranking_skinatac.csv'))
prom = pd.read_csv(os.path.join(DL, 'KO_ranking.csv'))

CANON = ['MYC', 'SNAI2', 'SRF', 'SMAD3', 'STAT3', 'TWIST1', 'FOS', 'EGR1']  # canonical fibrotic/myofib drivers
LEAD = 'HES1'

fig = plt.figure(figsize=(12.6, 5.4))
gs = GridSpec(1, 2, width_ratios=[1.0, 1.05], wspace=0.34, figure=fig)

# ===== A: perturbation score of canonical drivers (control) + HES1 (skin GRN) =====
axA = fig.add_subplot(gs[0, 0])
show = ['MYC', 'SNAI2', 'SRF', 'STAT3', 'SMAD3', 'TWIST1', 'HES1']
d = skin[skin.gene.isin(show)].copy().set_index('gene').loc[show].reset_index()
d = d.iloc[::-1]  # HES1 at bottom
def kcolor(g):
    if g == LEAD:  return fs.LEAD_BLUE
    return fs.CTRL_COLOR
cols = [kcolor(g) for g in d.gene]
y = np.arange(len(d))
# ps_score is negative (KO reduces myofibroblast identity); plot magnitude to the left
axA.hlines(y, 0, d.ps_score, color=cols, lw=2.4, alpha=0.9)
axA.scatter(d.ps_score, y, color=cols, s=64, zorder=3, edgecolor='white', linewidth=0.7)
axA.set_yticks(y); axA.set_yticklabels(d.gene, fontsize=9)
for tick, g in zip(axA.get_yticklabels(), d.gene):
    if g == LEAD: tick.set_fontweight('bold'); tick.set_color(fs.LEAD_BLUE)
    elif g == 'SMAD3': tick.set_color(fs.CTRL_COLOR)
for yy, (_, r) in zip(y, d.iterrows()):
    axA.text(r.ps_score - 1.0, yy, fs.star(r.ps_qbh), ha='right', va='center',
             fontsize=9, fontweight='bold', color='#333')
axA.axvline(0, color='#888', lw=0.6)
axA.set_xlabel('KO perturbation score (skin-specific scATAC GRN)')
axA.set_title('Engine recovers canonical drivers; nominates HES1', fontsize=9.8)
axA.invert_xaxis()  # more-negative (stronger) to the right
from matplotlib.lines import Line2D
axA.legend(handles=[Line2D([0],[0], marker='o', color='w', markerfacecolor=fs.CTRL_COLOR, ms=8, label='canonical driver (control)'),
                    Line2D([0],[0], marker='o', color='w', markerfacecolor=fs.LEAD_BLUE, ms=8, label='HES1 (druggable Notch lead)')],
           loc='lower left', fontsize=7.4)
fs.despine(axA)
fs.panel_label(axA, 'a')

# ===== B: robustness — skin vs promoter prior (all 43 TFs) =====
axB = fig.add_subplot(gs[0, 1])
m = skin[['gene', 'ps_score']].merge(prom[['gene', 'ps_score']], on='gene', suffixes=('_skin', '_prom'))
r = np.corrcoef(m.ps_score_prom, m.ps_score_skin)[0, 1]
axB.scatter(m.ps_score_prom, m.ps_score_skin, s=34, color='#BDBDBD', edgecolor='#888', linewidth=0.4, zorder=2)
lim = [min(m.ps_score_prom.min(), m.ps_score_skin.min()) - 3, 2]
axB.plot(lim, lim, ls='--', color='#999', lw=0.9, zorder=1)
# highlight lead HES1 + positive-control SMAD3. MEF2C shown as a demoted lead — NOT a
# negative control on the KO axis (it is a strong significant KO hit, is_lead=True; it only
# serves as a negative control on the activity/trajectory/spatial axes).
hl = {'HES1': fs.LEAD_BLUE, 'SMAD3': fs.CTRL_COLOR, 'MEF2C': '#BDBDBD'}
for g, c in hl.items():
    row = m[m.gene == g]
    if len(row):
        axB.scatter(row.ps_score_prom, row.ps_score_skin, s=95, color=c, edgecolor='black',
                    linewidth=0.8, zorder=4)
        dx = 1.5 if g != 'HES1' else 1.5
        axB.annotate(g, (row.ps_score_prom.iloc[0], row.ps_score_skin.iloc[0]),
                     xytext=(6, -2), textcoords='offset points', fontsize=8.5,
                     fontweight='bold', color=c if g != 'MEF2C' else '#555')
axB.set_xlabel('KO perturbation score — generic promoter prior')
axB.set_ylabel('KO perturbation score — skin scATAC prior')
axB.set_title(f'KO effect is robust to GRN prior (r = {r:.2f})', fontsize=9.8)
axB.annotate('HES1: q=1.1e-10 (skin), q=8.7e-7 (promoter)\nsignificant on both priors',
             xy=(0.03, 0.06), xycoords='axes fraction', fontsize=7.4, style='italic', color='#333',
             va='bottom', ha='left')
fs.despine(axB)
fs.panel_label(axB, 'b')

fig.suptitle('In-silico knockout nominates HES1 as a druggable node, robustly across GRN priors',
             fontsize=11.5, fontweight='bold', y=1.02)
fs.save_fig(fig, OUT)
print('prior robustness r =', round(r, 3))
print('HES1 skin ps=%.2f q=%.1e | promoter ps=%.2f q=%.1e' % (
    skin.query("gene=='HES1'").ps_score.iloc[0], skin.query("gene=='HES1'").ps_qbh.iloc[0],
    prom.query("gene=='HES1'").ps_score.iloc[0], prom.query("gene=='HES1'").ps_qbh.iloc[0]))

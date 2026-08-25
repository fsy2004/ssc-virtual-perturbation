"""Figure 4: SCENIC+ reconstruction of HES1- and RBPJ-centred eRegulons
from skin snMultiome data (GSE312129).
Panels show eRegulon sizes, HES1 targets, cell-type activity and the HES1/RBPJ network.
"""
import os, sys
import numpy as np, pandas as pd, networkx as nx
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
ARCH = os.path.join(PROJECT, 'server_archive')
META = os.path.join(ARCH, 'tier1', 'scenicplus', 'scplus', 'eRegulon_metadata__direct_e_regulon_metadata.csv')
AUC = os.path.join(ARCH, 'scenicplus_figdata', 'direct_gene_based_AUC_by_celltype.csv')
OUT = os.path.join(PROJECT, '03_results_figures', 'figures', 'fig_scenic_engine2')

d = pd.read_csv(META)
nt = d.groupby('TF')['Gene'].nunique()
impcol = 'importance_TF2G' if 'importance_TF2G' in d.columns else 'importance_x_abs_rho'

fig = plt.figure(figsize=(13.6, 9.2))
gs = GridSpec(2, 3, height_ratios=[1.0, 1.15], width_ratios=[1.05, 0.9, 1.05], wspace=0.5, hspace=0.42, figure=fig)

# ===== A: cross-engine concordance =====
axA = fig.add_subplot(gs[0, 0])
NOTCH = {'HES1', 'RBPJ', 'HEY1', 'HEYL', 'HES4'}; POS = {'SMAD3'}; NEG = {'MEF2C'}
show = [t for t in ['SMAD3', 'MEF2C', 'HES1', 'RBPJ', 'FOSB', 'TWIST1', 'SNAI2', 'RUNX1', 'STAT3'] if t in nt.index]
vals = nt.loc[show].sort_values()
tcol = lambda t: fs.LEAD_BLUE if t in NOTCH else (fs.CTRL_COLOR if t in POS else ('#BDBDBD' if t in NEG else '#8C8C8C'))
y = np.arange(len(vals)); cols = [tcol(t) for t in vals.index]
axA.hlines(y, 0, vals.values, color=cols, lw=2.2); axA.scatter(vals.values, y, color=cols, s=52, zorder=3, edgecolor='white', linewidth=0.6)
axA.set_yticks(y); axA.set_yticklabels(vals.index, fontsize=8.5)
for tk, t in zip(axA.get_yticklabels(), vals.index):
    if t in NOTCH: tk.set_fontweight('bold'); tk.set_color(fs.LEAD_BLUE)
    elif t in POS: tk.set_color(fs.CTRL_COLOR)
for yy, v in zip(y, vals.values): axA.text(v + 6, yy, str(int(v)), va='center', ha='left', fontsize=7.6, color='#333')
axA.set_xlabel('eRegulon target genes (SCENIC+)'); axA.set_title('Reconstructed eRegulons', fontsize=9.8)
axA.set_xlim(0, max(vals.values) * 1.18)
axA.legend(handles=[Line2D([0],[0], marker='o', color='w', markerfacecolor=fs.LEAD_BLUE, ms=8, label='Notch factor'),
                    Line2D([0],[0], marker='o', color='w', markerfacecolor=fs.CTRL_COLOR, ms=8, label='SMAD3 (fibrotic ctrl)'),
                    Line2D([0],[0], marker='o', color='w', markerfacecolor='#8C8C8C', ms=8, label='other mainline TF')], loc='lower right', fontsize=6.8)
fs.despine(axA); fs.panel_label(axA, 'a')

# ===== B: HES1 eRegulon targets (Notch-feedback highlighted) =====
axB = fig.add_subplot(gs[0, 1])
htop = d[d.TF == 'HES1'].sort_values(impcol, ascending=False).drop_duplicates('Gene').head(12).iloc[::-1]
NOTCHTGT = {'NRARP', 'HES5', 'HES1', 'HEY1', 'HEYL', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'DTX1'}
y = np.arange(len(htop)); cols = [fs.LEAD_BLUE if g in NOTCHTGT else '#8C8C8C' for g in htop.Gene]
axB.hlines(y, 0, htop[impcol], color=cols, lw=2.2); axB.scatter(htop[impcol], y, color=cols, s=48, zorder=3, edgecolor='white', linewidth=0.6)
axB.set_yticks(y); axB.set_yticklabels(htop.Gene, fontsize=8)
for tk, g in zip(axB.get_yticklabels(), htop.Gene):
    if g in NOTCHTGT: tk.set_fontweight('bold'); tk.set_color(fs.LEAD_BLUE)
axB.set_xlabel('TF→gene importance'); axB.set_title('HES1 eRegulon targets', fontsize=9.8)
axB.annotate('Notch-feedback\nNRARP, HES5', xy=(0.96, 0.06), xycoords='axes fraction', ha='right', fontsize=7.2, style='italic', color='#2166AC', va='bottom')
axB.set_xlim(0, htop[impcol].max() * 1.15); fs.despine(axB); fs.panel_label(axB, 'b')

# ===== C: eRegulon activity by cell type =====
axC = fig.add_subplot(gs[0, 2])
a = pd.read_csv(AUC, index_col=0); cts = [c for c in a.columns if c != 'n_eRegulon_targets']
rows = {'HES1': 'HES1_direct_+/-_(53g)', 'RBPJ': 'RBPJ_direct_+/+_(39g)', 'SMAD3': 'SMAD3_direct_+/+_(145g)'}
mark = {'HES1': (fs.LEAD_BLUE, 'o'), 'RBPJ': ('#4C72B0', 's'), 'SMAD3': (fs.CTRL_COLOR, '^')}
order = sorted(cts, key=lambda c: -a.loc[rows['HES1'], c]); xpos = np.arange(len(order))
for lab, rn in rows.items():
    if rn in a.index:
        axC.plot(xpos, [a.loc[rn, c] for c in order], marker=mark[lab][1], color=mark[lab][0], lw=1.6, ms=7, markeredgecolor='white', markeredgewidth=0.6, label=lab)
axC.set_xticks(xpos); axC.set_xticklabels(order, rotation=35, ha='right', fontsize=7.8)
axC.set_ylabel('mean eRegulon activity (AUCell)'); axC.set_title('eRegulon activity by cell type', fontsize=9.8)
axC.legend(fontsize=7.6, loc='upper right')
axC.annotate('active in fibroblasts among\nother lineages (broad Notch)', xy=(0.03, 0.04), xycoords='axes fraction', fontsize=7.0, style='italic', color='#666', va='bottom')
fs.despine(axC); fs.panel_label(axC, 'c')

# ===== D: reconstructed HES1 / RBPJ Notch eRegulon network =====
axD = fig.add_subplot(gs[1, :])
TFS = ['HES1', 'RBPJ']; TOPN = 13
G = nx.DiGraph(); pos = {}; HUB = {'HES1': np.array([-1.6, 0.0]), 'RBPJ': np.array([1.6, 0.0])}
for tf in TFS:
    G.add_node(tf, kind='tf'); pos[tf] = HUB[tf]
    top = d[d.TF == tf].sort_values(impcol, ascending=False).drop_duplicates('Gene').head(TOPN)
    side = -1 if tf == 'HES1' else 1; n = len(top)
    for i, (_, r) in enumerate(top.iterrows()):
        g = r['Gene']; G.add_node(g, kind='notch_tgt' if g in NOTCHTGT else 'tgt'); G.add_edge(tf, g, w=float(r[impcol]))
        ang = np.pi * (0.5 - (i + 0.5) / n); pos[g] = HUB[tf] + np.array([side * 1.15 * np.cos(ang), 1.15 * np.sin(ang)])
ews = [G[u][v]['w'] for u, v in G.edges()]; wmax = max(ews)
nx.draw_networkx_edges(ax=axD, G=G, pos=pos, width=[0.5 + 2.2 * w / wmax for w in ews], edge_color='#CFCFCF', alpha=0.75, arrows=True, arrowsize=7, node_size=430)
tf_n = [n for n in G if G.nodes[n]['kind'] == 'tf']; nt_n = [n for n in G if G.nodes[n]['kind'] == 'notch_tgt']; tg_n = [n for n in G if G.nodes[n]['kind'] == 'tgt']
nx.draw_networkx_nodes(ax=axD, G=G, pos=pos, nodelist=tf_n, node_size=1150, node_color=fs.LEAD_BLUE, edgecolors='white', linewidths=1.4)
nx.draw_networkx_nodes(ax=axD, G=G, pos=pos, nodelist=nt_n, node_size=430, node_color='#2166AC', edgecolors='white', linewidths=0.9)
nx.draw_networkx_nodes(ax=axD, G=G, pos=pos, nodelist=tg_n, node_size=320, node_color='#D0D0D0', edgecolors='white', linewidths=0.9)
nx.draw_networkx_labels(ax=axD, G=G, pos=pos, labels={n: n for n in tf_n}, font_size=10, font_weight='bold', font_color='white')
nx.draw_networkx_labels(ax=axD, G=G, pos=pos, labels={n: n for n in nt_n + tg_n}, font_size=6.6, font_color='#222')
axD.set_title('Reconstructed HES1 / RBPJ Notch eRegulon network', fontsize=9.8); axD.axis('off')
axD.legend(handles=[Line2D([0],[0], marker='o', color='w', markerfacecolor=fs.LEAD_BLUE, ms=11, label='Notch TF'),
                    Line2D([0],[0], marker='o', color='w', markerfacecolor='#2166AC', ms=8, label='Notch-feedback target'),
                    Line2D([0],[0], marker='o', color='w', markerfacecolor='#D0D0D0', ms=7, label='other target')], loc='lower center', ncol=3, fontsize=7.4, frameon=False, bbox_to_anchor=(0.5, -0.02))
fs.panel_label(axD, 'd')

fig.suptitle('Engine 2 — a multiome eGRN recovers HES1 and the core Notch factor RBPJ',
             fontsize=11, fontweight='bold', y=1.0)
fs.save_fig(fig, OUT)
print('HES1:', int(nt['HES1']), 'RBPJ:', int(nt['RBPJ']), 'SMAD3:', int(nt['SMAD3']), '| network nodes:', G.number_of_nodes())

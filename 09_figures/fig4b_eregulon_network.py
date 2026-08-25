"""Fig 4 (advanced) — the HES1 / RBPJ Notch eRegulon module reconstructed by SCENIC+ (Engine 2).
Network of the two Notch transcription factors and their top target genes; canonical Notch-feedback
targets (NRARP, HES5) highlighted. Data: SCENIC+ direct eRegulon metadata (verified). Vector PDF.
"""
import os, sys
import numpy as np, pandas as pd, networkx as nx
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import figstyle as fs

META = r'C:\Users\fsy\Desktop\SSc_virtualKO_MR\server_archive\tier1\scenicplus\scplus\eRegulon_metadata__direct_e_regulon_metadata.csv'
OUT = r'C:\Users\fsy\Desktop\SSc_virtualKO_MR\03_results_figures\figures\Fig4b_eregulon_network'
TFS = ['HES1', 'RBPJ']
NOTCH_TGT = {'NRARP', 'HES5', 'HEY1', 'HEYL', 'DTX1', 'NOTCH1', 'NOTCH2', 'NOTCH3'}
TOPN = 14

d = pd.read_csv(META)
impcol = 'importance_TF2G' if 'importance_TF2G' in d.columns else 'importance_x_abs_rho'
G = nx.DiGraph(); pos = {}
HUB = {'HES1': np.array([-1.35, 0.0]), 'RBPJ': np.array([1.35, 0.0])}
for tf in TFS:
    G.add_node(tf, kind='tf'); pos[tf] = HUB[tf]
    top = d[d.TF == tf].sort_values(impcol, ascending=False).drop_duplicates('Gene').head(TOPN)
    genes = list(top['Gene'])
    side = -1 if tf == 'HES1' else 1
    n = len(genes)
    for i, (_, r) in enumerate(top.iterrows()):
        g = r['Gene']
        G.add_node(g, kind='notch_tgt' if g in NOTCH_TGT else 'tgt')
        G.add_edge(tf, g, w=float(r[impcol]))
        # fan the targets on the TF's side (semicircle facing outward)
        ang = np.pi * (0.5 - (i + 0.5) / n)          # -90..+90
        rad = 1.35
        pos[g] = HUB[tf] + np.array([side * rad * np.cos(ang), rad * np.sin(ang)])

fig, ax = plt.subplots(figsize=(8.4, 6.6))
ews = [G[u][v]['w'] for u, v in G.edges()]
wmax = max(ews)
nx.draw_networkx_edges(ax=ax, G=G, pos=pos, width=[0.6 + 2.6 * w / wmax for w in ews],
                       edge_color='#C9C9C9', alpha=0.7, arrows=True, arrowsize=8, node_size=650)
tf_nodes = [n for n in G if G.nodes[n]['kind'] == 'tf']
nt_nodes = [n for n in G if G.nodes[n]['kind'] == 'notch_tgt']
tg_nodes = [n for n in G if G.nodes[n]['kind'] == 'tgt']
nx.draw_networkx_nodes(ax=ax, G=G, pos=pos, nodelist=tf_nodes, node_size=1500, node_color=fs.LEAD_BLUE, edgecolors='white', linewidths=1.5)
nx.draw_networkx_nodes(ax=ax, G=G, pos=pos, nodelist=nt_nodes, node_size=560, node_color='#2166AC', edgecolors='white', linewidths=1.0)
nx.draw_networkx_nodes(ax=ax, G=G, pos=pos, nodelist=tg_nodes, node_size=420, node_color='#BDBDBD', edgecolors='white', linewidths=1.0)
nx.draw_networkx_labels(ax=ax, G=G, pos=pos, labels={n: n for n in tf_nodes}, font_size=11, font_weight='bold', font_color='white')
nx.draw_networkx_labels(ax=ax, G=G, pos=pos, labels={n: n for n in nt_nodes + tg_nodes}, font_size=7.5,
                        font_weight='normal', font_color='#222')
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0], marker='o', color='w', markerfacecolor=fs.LEAD_BLUE, markersize=13, label='Notch TF'),
                   Line2D([0],[0], marker='o', color='w', markerfacecolor='#2166AC', markersize=9, label='Notch-feedback target'),
                   Line2D([0],[0], marker='o', color='w', markerfacecolor='#BDBDBD', markersize=8, label='other target')],
          loc='lower center', ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.08))
ax.set_title('SCENIC+ Notch eRegulon module (HES1, RBPJ) — top targets', fontsize=11, fontweight='bold')
ax.axis('off')
fs.save_fig(fig, OUT)
print('nodes:', G.number_of_nodes(), 'edges:', G.number_of_edges(), '| Notch-feedback targets shown:', [n for n in nt_nodes])

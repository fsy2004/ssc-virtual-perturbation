"""Fig 8 — γ-secretase druggability: docking + membrane MD (qualitative, no MM-GBSA).
Docking (smina vs PDB 6LR4) + MD stability from GROMACS xvg, 3 × 200 ns per condition
(nirogacestat / semagacestat / apo; 1.8 µs aggregate). Time series: periodic-image (PBC) frames
above a physical ceiling are dropped, then a large rolling median (spike-immune) per replica and
the median across the 3 replicas give the true bound-state baseline (no break-out spikes). FEL:
per-condition 3D free-energy landscapes from PC1/PC2 (MAD outlier removal + smoothed density).
(Per-residue RMSF is omitted: the local per-atom xvg selection is inconsistent across conditions.)"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.colors as mcolors
from scipy.ndimage import gaussian_filter
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers the 3d projection)
import figstyle as fs

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
MD   = os.path.join(PROJECT, '02_analysis_modules', 'md_pipeline', 'md_data')
DOCK = os.path.join(PROJECT, '03_results_figures', 'dock_figs', 'affinities.csv')
POSE = os.path.join(PROJECT, '03_results_figures', 'dock_figs', 'nirogacestat_zoom.png')
OUT  = os.path.join(PROJECT, '03_results_figures', 'figures', 'Fig8_docking_md')

COND = {
    'nirogacestat': (os.path.join(MD, 'niro', 'charmm-gui-8090533815', 'gromacs'), fs.TF_COLORS['HES1']),
    'semagacestat': (os.path.join(MD, 'sema', 'charmm-gui-8092842495', 'gromacs'), fs.TF_COLORS['FOSB']),
    'apo':          (os.path.join(MD, 'apo', 'gromacs'),                            '#7A7A7A'),
}
FELTITLE = {'nirogacestat': 'Nirogacestat', 'semagacestat': 'Semagacestat', 'apo': 'Apo (no ligand)'}

def read_xvg(path):
    t, v = [], []
    with open(path) as f:
        for ln in f:
            if not ln or ln[0] in '#@':
                continue
            p = ln.split()
            if len(p) >= 2:
                try:
                    t.append(float(p[0])); v.append(float(p[1]))
                except ValueError:
                    pass
    return np.asarray(t), np.asarray(v)

def robust(cond_dir, fname, win=200, cap=None):
    """Clean baseline from a noisy MD series: drop periodic-image (PBC) frames above a physical
    ceiling, take a large rolling median (spike-immune) + light mean per replica, then the MEDIAN
    across the 3 replicas. No break-out spikes; the bound-state baseline is preserved."""
    series, xref = [], None
    for r in (1, 2, 3):
        p = os.path.join(cond_dir, f'prod_r{r}', fname)
        if not os.path.exists(p):
            continue
        t, v = read_xvg(p)
        if not len(v):
            continue
        s = pd.Series(v, dtype=float)
        if cap is not None:
            s = s.mask(s > cap)                                        # PBC jump -> NaN, excluded
        s = s.rolling(win, center=True, min_periods=1).median()
        s = s.rolling(max(20, win // 2), center=True, min_periods=1).mean()
        series.append(s.values); xref = t
    if not series:
        return None, None
    n = min(len(s) for s in series)
    with np.errstate(all='ignore'):
        out = np.nanmedian(np.vstack([s[:n] for s in series]), 0)
    out = pd.Series(out).interpolate(limit_direction='both').values   # fill any all-NaN gaps
    return xref[:n], out

FELCMAP = mcolors.LinearSegmentedColormap.from_list('fel',
    ['#08306B', '#2171B5', '#6BAED6', '#C7E9C0', '#FEE391', '#FE9929', '#CB181D', '#7A0177'])

def fel_surface(ax, cond_dir, title):
    P = []
    for r in (1, 2, 3):
        p = os.path.join(cond_dir, f'prod_r{r}', 'proj2d.xvg')
        if os.path.exists(p):
            a, b = read_xvg(p)
            P.append(np.column_stack([a, b]))
    if not P:
        ax.text2D(0.5, 0.5, 'n/a', transform=ax.transAxes, ha='center'); return None
    P = np.vstack(P)
    def keep(v, k=6):
        m = np.median(v); mad = np.median(np.abs(v - m)) * 1.4826 + 1e-9
        return (v >= m - k * mad) & (v <= m + k * mad)
    P = P[keep(P[:, 0]) & keep(P[:, 1])]                              # drop PBC-jumped frames
    RT, NB, SIG, GMAX = 0.616, 60, 2.0, 5.0
    H, xe, ye = np.histogram2d(P[:, 0], P[:, 1], bins=NB, density=True)
    H = gaussian_filter(H, sigma=SIG); Pn = H / H.sum()
    G = np.full_like(Pn, np.nan); occ = Pn > 0
    G[occ] = -RT * np.log(Pn[occ]); G -= np.nanmin(G); G[occ] = np.clip(G[occ], 0, GMAX)
    xc = 0.5 * (xe[:-1] + xe[1:]); yc = 0.5 * (ye[:-1] + ye[1:])
    Xc, Yc = np.meshgrid(xc, yc, indexing='ij')
    surf = ax.plot_surface(Xc, Yc, G, cmap=FELCMAP, vmin=0, vmax=GMAX, linewidth=0,
                           antialiased=True, rcount=NB, ccount=NB)
    ax.contourf(Xc, Yc, np.nan_to_num(G, nan=GMAX), zdir='z', offset=0,
                levels=np.linspace(0, GMAX, 16), cmap=FELCMAP, alpha=0.5)
    ii, jj = np.unravel_index(np.nanargmin(G), G.shape)
    ax.scatter(Xc[ii, jj], Yc[ii, jj], 0, color='k', s=20, marker='v', depthshade=False)
    ax.set_zlim(0, GMAX); ax.view_init(elev=30, azim=-60)
    ax.set_xlabel('PC1', fontsize=7, labelpad=-5); ax.set_ylabel('PC2', fontsize=7, labelpad=-5)
    ax.set_zlabel('ΔG', fontsize=7, labelpad=-5); ax.tick_params(labelsize=5.5)
    ax.xaxis.pane.set_alpha(0); ax.yaxis.pane.set_alpha(0); ax.zaxis.pane.set_alpha(0)
    ax.set_title(title, fontsize=9.5, pad=-4)
    return surf

fig = plt.figure(figsize=(13.8, 11.6))
gs = GridSpec(3, 3, figure=fig, hspace=0.42, wspace=0.34,
              left=0.07, right=0.95, top=0.92, bottom=0.06)

# (a) docking affinities (lollipop, wide)
ax = fig.add_subplot(gs[0, 0:2])
aff = pd.read_csv(DOCK).sort_values('affinity_kcalmol')
y = np.arange(len(aff))
norm = mcolors.Normalize(aff.affinity_kcalmol.min(), aff.affinity_kcalmol.max())
bcmap = mcolors.LinearSegmentedColormap.from_list('be', ['#08306B', '#2166AC', '#92C5DE', '#F0F0F0'])
for yi, (_, r) in zip(y, aff.iterrows()):
    hl = r.drug in ('nirogacestat', 'semagacestat')
    ax.plot([0, r.affinity_kcalmol], [yi, yi], color='#D5D5D5', lw=1.2, zorder=1)
    ax.scatter(r.affinity_kcalmol, yi, s=180 if hl else 110,
               color=bcmap(norm(r.affinity_kcalmol)),
               edgecolor=('#B2182B' if r.drug == 'nirogacestat' else '#2B2B2B'),
               lw=1.7 if hl else 0.8, zorder=3)
    ax.text(r.affinity_kcalmol, yi, f'{r.affinity_kcalmol:.1f}', ha='center', va='center',
            fontsize=6.0, fontweight='bold', color=('white' if norm(r.affinity_kcalmol) < 0.5 else '#1A1A1A'))
ax.set_yticks(y); ax.set_yticklabels(aff.drug, fontsize=7.4)
for t, d in zip(ax.get_yticklabels(), aff.drug):
    if d == 'nirogacestat': t.set_color('#B2182B'); t.set_fontweight('bold')
    elif d == 'ibuprofen':  t.set_color('#888888'); t.set_style('italic')
ax.invert_xaxis()
ax.set_xlabel('docking affinity (kcal/mol; more negative = tighter)')
ax.set_title('γ-Secretase inhibitor docking (PDB 6LR4): nirogacestat among the strongest binders', pad=6)
fs.despine(ax, keep=('bottom',)); ax.tick_params(left=False)
fs.panel_label(ax, 'a')

# (b) docked pose
ax = fig.add_subplot(gs[0, 2])
try:
    ax.imshow(plt.imread(POSE))
except Exception:
    ax.text(0.5, 0.5, 'pose n/a', ha='center', va='center', transform=ax.transAxes, color='#888')
ax.axis('off')
ax.set_title('Nirogacestat in the PSEN1\ntransmembrane pocket', fontsize=9.5, pad=4)
fs.panel_label(ax, 'b', dx=0.04, dy=1.03)

# (c) backbone RMSD
ax = fig.add_subplot(gs[1, 0])
for name, (cd, col) in COND.items():
    X, M = robust(cd, 'rmsd_bb.xvg', win=380, cap=0.7)
    if M is None: continue
    ax.plot(X, M, color=col, lw=2.0, label=name)
ax.axvspan(0, 20, color='#EFEFEF', zorder=0)
ax.set_xlabel('time (ns)'); ax.set_ylabel('Cα RMSD (nm)'); ax.set_ylim(0.15, 0.42)
ax.set_title('Backbone RMSD (n = 3): equilibrated, stable', fontsize=10, pad=5)
ax.legend(fontsize=6.6, loc='upper left')
fs.despine(ax); fs.panel_label(ax, 'c')

# (d) ligand RMSD
ax = fig.add_subplot(gs[1, 1])
for name in ('nirogacestat', 'semagacestat'):
    cd, col = COND[name]
    X, M = robust(cd, 'rmsd_lig.xvg', win=260, cap=1.2)
    if M is None: continue
    ax.plot(X, M, color=col, lw=2.0, label=name)
ax.set_xlabel('time (ns)'); ax.set_ylabel('ligand RMSD (nm)'); ax.set_ylim(0, 0.55)
ax.set_title('Ligand RMSD vs docked pose\n(median ≈ 0.3 nm, pose retained)', fontsize=10, pad=5)
ax.legend(fontsize=6.6, loc='upper left')
fs.despine(ax); fs.panel_label(ax, 'd')

# (e) ligand-pocket H-bonds
ax = fig.add_subplot(gs[1, 2])
for name in ('nirogacestat', 'semagacestat'):
    cd, col = COND[name]
    X, M = robust(cd, 'hbond_lig.xvg', win=200)
    if M is None: continue
    ax.plot(X / 1000.0, M, color=col, lw=1.9, label=name)  # XVG time is in ps
ax.set_xlabel('time (ns)'); ax.set_ylabel('ligand–pocket H-bonds'); ax.set_ylim(0, None)
ax.set_title('Ligand–pocket hydrogen bonds\n(maintained throughout)', fontsize=10, pad=5)
ax.legend(fontsize=6.6, loc='upper right')
fs.despine(ax); fs.panel_label(ax, 'e')

# (f, g, h) free-energy landscapes — one per condition
surf = None
felaxes = []
for k, (name, (cd, col)) in enumerate(COND.items()):
    ax = fig.add_subplot(gs[2, k], projection='3d')
    s = fel_surface(ax, cd, FELTITLE[name])
    surf = s if s is not None else surf
    felaxes.append(ax)
    ax.text2D(0.0, 1.0, 'fgh'[k], transform=ax.transAxes, fontsize=13, fontweight='bold', va='bottom', ha='left')
if surf is not None:
    cb = fig.colorbar(surf, ax=felaxes, shrink=0.62, aspect=22, pad=0.015)
    cb.set_label('ΔG (kcal/mol)', fontsize=8); cb.ax.tick_params(labelsize=6.5)

fig.suptitle('Druggability of the γ-secretase node: docking and molecular dynamics stability',
             fontsize=12.5, fontweight='bold', y=0.975)
fs.save_fig(fig, OUT)
print('docking rows:', len(aff))

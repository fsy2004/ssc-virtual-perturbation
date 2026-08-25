"""Shared publication style for SSc virtual-KO figures.
Nature/Cell-grade matplotlib aesthetics: clean, minimal, bold labels, panel borders,
curated low-saturation professional palette, consistent sizing.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

# ---- global rcParams ----
mpl.rcParams.update({
    'figure.dpi': 120,
    'savefig.dpi': 400,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9,
    'axes.titlesize': 10.5,
    'axes.labelsize': 9.5,
    'axes.titleweight': 'bold',
    'axes.labelweight': 'bold',
    'axes.linewidth': 1.0,
    'axes.edgecolor': '#2B2B2B',
    'xtick.labelsize': 8.5,
    'ytick.labelsize': 8.5,
    'xtick.major.width': 0.9,
    'ytick.major.width': 0.9,
    'xtick.minor.width': 0.6,
    'ytick.minor.width': 0.6,
    'xtick.major.size': 4.0,
    'ytick.major.size': 4.0,
    'xtick.color': '#2B2B2B',
    'ytick.color': '#2B2B2B',
    'legend.fontsize': 8,
    'legend.frameon': False,
    'legend.handlelength': 1.4,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,   # editable text in vector output
})

# ---- curated professional palette ----
# leads in strong colours; MEF2C (demoted) muted grey; positive controls green
TF_COLORS = {
    'HES1':  '#2166AC',   # Notch lead — steel blue
    'FOSB':  '#D6604D',   # AP-1 lead — muted orange-red
    'CEBPD': '#92C5DE',   # weak — light blue
    'MEF2C': '#BDBDBD',   # demoted — light grey
    'SMAD3': '#4DAC26',   # canonical fibrotic control — mid green
    'CEBPB': '#4DAC26',
    'KLF4':  '#4DAC26',
    'FOS':   '#4DAC26',
    'STAT3': '#4DAC26',
}
CTRL_COLOR   = '#4DAC26'    # positive-control TFs
LEAD_BLUE    = '#2166AC'    # primary lead
LEAD_RED     = '#D6604D'    # secondary lead
GROUP        = {'SSc': '#B2182B', 'HC': '#2166AC', 'Ctrl': '#2166AC', 'control': '#2166AC'}
# Diverging blue-white-red for heatmaps
SEQ_BLUE_RED = ['#053061', '#2166AC', '#92C5DE', '#F7F7F7', '#FDDBC7', '#D6604D', '#67001F']
# Continuous expression / activity colormap (perceptually uniform, colorblind-safe)
EXPR_CMAP    = 'viridis'

# ---- helper functions ----

def despine(ax, keep=('left', 'bottom')):
    """Remove top/right spines; optionally remove all but specified sides."""
    for s in ('top', 'right', 'left', 'bottom'):
        ax.spines[s].set_visible(s in keep)

def panel_label(ax, letter, size=13, **kwargs):
    """Bold panel label at a UNIFORM offset from each axes' top-left corner.
    Uses points offset (not axes fractions) so labels align on one line across a row
    and in a column regardless of y-axis labels/colourbars. Legacy dx/dy kwargs are
    accepted and ignored to keep the alignment uniform across all figures."""
    ax.annotate(letter.upper(), xy=(0.0, 1.0), xycoords='axes fraction',
                xytext=(-26, 8), textcoords='offset points',
                fontsize=size, fontweight='bold', va='bottom', ha='left',
                annotation_clip=False, zorder=20)

def star(p):
    """Return significance stars string."""
    if p < 0.001:  return '***'
    if p < 0.01:   return '**'
    if p < 0.05:   return '*'
    return 'ns'

def add_stat_annotation(ax, x, y_top, p, fontsize=9):
    """Place a significance star above a violin/bar at position x."""
    ax.text(x, y_top, star(p), ha='center', va='bottom',
            fontsize=fontsize, fontweight='bold', color='#222')

def save_fig(fig, stem):
    """Save PNG + PDF with consistent settings."""
    fig.savefig(f'{stem}.png', dpi=400, bbox_inches='tight', facecolor='white')
    fig.savefig(f'{stem}.pdf', bbox_inches='tight', facecolor='white')
    print(f'saved {stem}.png / .pdf')

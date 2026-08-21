"""Compact (40 x 50 mm) alternatives to alt3_sample_grouped_bars.

Same information: per-sample mean variant coverage (4 bins summed), colored by
assay. Instead of one labelled bar per sample, each sample is one point.
"""
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, FuncFormatter, NullFormatter
from pathlib import Path

REPO = Path('/Users/mkh/GitHub/mor_dms_analysis')
OUT_DIR = REPO / 'figures' / 'variant_counts' / 'plots'
MM = 1 / 25.4
W_MM, H_MM = 40, 40

ASSAY_ORDER = ['damgo_drc', 'morphine_drc', 'fentanyl_drc', 'multidrug', 'surface']
ASSAY_DISPLAY = {'damgo_drc': 'DAMGO', 'morphine_drc': 'Morphine',
                 'fentanyl_drc': 'Fentanyl', 'multidrug': 'Multiligand',
                 'surface': 'Surface'}
COLORS = {'damgo_drc': '#1f77b4', 'morphine_drc': '#ff7f0e',
          'fentanyl_drc': '#2ca02c', 'multidrug': '#d62728',
          'surface': '#9467bd'}

STYLE = {'font.family': 'Helvetica', 'font.size': 6, 'text.color': 'black',
         'axes.labelcolor': 'black', 'axes.edgecolor': 'black',
         'xtick.color': 'black', 'ytick.color': 'black',
         'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
         'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5,
         'patch.linewidth': 0.5, 'lines.linewidth': 0.5, 'pdf.fonttype': 42}

# ---- per-sample coverage: sum the 4 bin means for each sample ----
tbl = pd.read_csv(OUT_DIR / 'per_bin_coverage_summary.csv')
g = (tbl.groupby(['assay', 'date', 'replicate_display', 'sample', 'conc'], dropna=False)
        .agg(cov=('mean_coverage', 'sum'), n_bins=('bin', 'count'))
        .reset_index())
g = g[g['assay'].isin(ASSAY_ORDER)]
print(f"{len(g)} samples, coverage {g['cov'].min():.0f}-{g['cov'].max():.0f}x")

XLIM = (330, 1900)
COV_TICKS = [500, 1000]
# Log can't reach 0, so the log strip is anchored a full decade below the data;
# the symlog variant keeps a true 0 tick via a linear region below LINTHRESH.
STRIP_XLIM = (1, 2000)
STRIP_TICKS = [1, 10, 100, 1000]
SYMLOG_XLIM = (0, 2000)
SYMLOG_TICKS = [0, 100, 1000]
LINTHRESH = 100
tick_fmt = FuncFormatter(lambda v, _: f'{v:g}')


def finish(fig, ax, name):
    # matplotlib labels log minor ticks by default; at 40 mm they overlap
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_minor_formatter(NullFormatter())
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_linewidth(0.5)
    ax.tick_params(labelsize=6, colors='black', width=0.5, pad=1.5, length=2)
    ax.tick_params(which='minor', length=1, width=0.5)
    fig.set_size_inches(W_MM * MM, H_MM * MM)
    fig.savefig(OUT_DIR / f'{name}.pdf')
    fig.savefig(OUT_DIR / f'{name}.png', dpi=600)
    plt.close(fig)
    print(f'wrote: {name}.{{pdf,png}}')


# ---------- A: one strip per assay, labels inside the panel ----------
def draw_strip(name, scale):
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots()
        fig.subplots_adjust(left=0.06, right=0.99, top=1.0, bottom=0.155)
        rng = np.random.default_rng(0)
        for i, assay in enumerate(ASSAY_ORDER):
            sub = g[g['assay'] == assay]
            y = i + rng.uniform(-0.13, 0.13, len(sub))
            ax.scatter(sub['cov'], y, c=COLORS[assay], s=1.6, alpha=0.75,
                       edgecolors='none', clip_on=False, zorder=3)
            med = sub['cov'].median()
            ax.plot([med, med], [i - 0.22, i + 0.22], color='black', lw=0.5, zorder=4)
            ax.text(0.01, i - 0.30, f'{ASSAY_DISPLAY[assay]} ({len(sub)})',
                    transform=ax.get_yaxis_transform(),
                    fontsize=6, color='black', ha='left', va='bottom')
        if scale == 'symlog':
            ax.set_xscale('symlog', linthresh=LINTHRESH, linscale=0.35)
            ax.set_xlim(*SYMLOG_XLIM)
            ax.set_xticks(SYMLOG_TICKS)
        else:
            ax.set_xscale('log')
            ax.set_xlim(*STRIP_XLIM)
            ax.set_xticks(STRIP_TICKS)
        # tight top/bottom padding: just enough for the first label and last dot row
        ax.set_ylim(len(ASSAY_ORDER) - 0.72, -0.62)
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
        ax.xaxis.set_major_formatter(tick_fmt)
        ax.set_xlabel('Mean sequencing coverage (×)', fontsize=6, color='black',
                      labelpad=1)
        finish(fig, ax, name)


draw_strip('compact_coverage_strip', 'log')
draw_strip('compact_coverage_strip_symlog', 'symlog')

# ---------- B: rank-ordered coverage curve ----------
with plt.rc_context(STYLE):
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0.26, right=0.96, top=0.99, bottom=0.16)
    s = g.sort_values('cov').reset_index(drop=True)
    ax.scatter(np.arange(1, len(s) + 1), s['cov'],
               c=[COLORS[a] for a in s['assay']], s=1.6, alpha=0.85,
               edgecolors='none', zorder=3)
    ax.set_yscale('log')
    ax.set_ylim(*XLIM)
    ax.set_xlim(-3, len(s) + 4)
    ax.set_yticks(COV_TICKS)
    ax.yaxis.set_major_formatter(tick_fmt)
    ax.set_xticks([1, 50, 100, len(s)])
    ax.set_xlabel('Sample (ranked)', fontsize=6, color='black', labelpad=1)
    ax.set_ylabel('Mean coverage (×)', fontsize=6, color='black', labelpad=1)
    finish(fig, ax, 'compact_coverage_rank')

# ---------- C: stacked histogram of log10 coverage ----------
with plt.rc_context(STYLE):
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0.20, right=0.99, top=0.99, bottom=0.16)
    bins = np.logspace(np.log10(300), np.log10(2000), 16)
    ax.hist([g.loc[g['assay'] == a, 'cov'] for a in ASSAY_ORDER], bins=bins,
            stacked=True, color=[COLORS[a] for a in ASSAY_ORDER],
            edgecolor='none')
    ax.set_xscale('log')
    ax.set_xlim(300, 2000)
    ax.set_xticks(COV_TICKS)
    ax.xaxis.set_major_formatter(tick_fmt)
    ax.set_xlabel('Mean coverage (×)', fontsize=6, color='black', labelpad=1)
    ax.set_ylabel('Samples', fontsize=6, color='black', labelpad=1)
    finish(fig, ax, 'compact_coverage_hist')

print('done')

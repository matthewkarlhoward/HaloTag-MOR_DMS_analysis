"""Cell-count vs gDNA recovery + per-sample variant coverage plots."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

REPO = Path('/Users/mkh/GitHub/mor_dms_analysis')
COUNTS_ROOT = REPO / 'figures' / 'variant_counts'
OUT_DIR = COUNTS_ROOT / 'plots'
OUT_DIR.mkdir(parents=True, exist_ok=True)

PG_PER_CELL = 6
MM = 1 / 25.4

# ---- Load spreadsheet ----
df = pd.read_excel(REPO / 'figures' / '20260522_sort_fastq_gDNA_info.xlsx', sheet_name='samples')
df = df.rename(columns={c: c.strip().replace('\xa0','') for c in df.columns})
df = df.rename(columns={'concentration (log , M)': 'conc',
                        'tube#': 'tube',
                        'ug dna': 'ug_dna',
                        'gDNA (ng/uL)': 'gdna_conc'})

def assay_for(row):
    d = int(row['date']); samp = row['sample']
    if d < 20240000:
        if samp == 'Surface expression': return 'surface'
        if samp == 'DAMGO': return 'damgo_drc'
        return 'other'
    if 20251017 <= d <= 20251019: return 'morphine_drc'
    if 20251009 <= d <= 20251011: return 'fentanyl_drc'
    if 20250805 <= d <= 20250816: return 'multidrug'
    return 'other'

df['assay'] = df.apply(assay_for, axis=1)

def sample_label(row):
    r, s = int(row['replicate']), row['sample']
    c = row['conc']
    s_short = {'Mitragynine Pseudoindoxyl':'MP','Surface expression':'Surface'}.get(s, s)
    if pd.isna(c): return f'R{r}_{s_short}'
    return f'R{r}_{s_short}_{c:g}M'

df['sample_id'] = df.apply(sample_label, axis=1)

ASSAY_ORDER = ['damgo_drc','morphine_drc','fentanyl_drc','multidrug','surface','other']
ASSAY_DISPLAY = {'damgo_drc':'damgo_drc', 'morphine_drc':'morphine_drc',
                 'fentanyl_drc':'fentanyl_drc', 'multidrug':'multiligand',
                 'surface':'surface', 'other':'other'}
COLORS = {
    'damgo_drc':    '#1f77b4',
    'morphine_drc': '#ff7f0e',
    'fentanyl_drc': '#2ca02c',
    'multidrug':    '#d62728',
    'surface':      '#9467bd',
    'other':        '#7f7f7f',
}

STYLE = {'font.family': 'Helvetica', 'font.size': 6, 'text.color': 'black',
         'axes.labelcolor': 'black', 'axes.edgecolor': 'black',
         'xtick.color': 'black', 'ytick.color': 'black',
         'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
         'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5,
         'patch.linewidth': 0.5, 'lines.linewidth': 0.5}

# ---------- Plot 1: cells vs ug DNA ----------
from scipy.stats import pearsonr
N_OUTLIERS = 2
threshold = df['cells'].nlargest(N_OUTLIERS).min()
df_plot = df[df['cells'] < threshold]
x_max = df_plot['cells'].max() * 1.05
r, p = pearsonr(df_plot['cells'], df_plot['ug_dna'])

with plt.rc_context(STYLE):
    fig, ax = plt.subplots(figsize=(40*MM, 40*MM))
    for assay in ASSAY_ORDER:
        sub = df_plot[df_plot['assay'] == assay]
        full = df[df['assay'] == assay]
        if len(full) == 0: continue
        ax.scatter(sub['cells'], sub['ug_dna'], c=COLORS[assay],
                   label=f'{ASSAY_DISPLAY[assay]} (n={len(full)})',
                   alpha=0.6, s=2.5, edgecolors='none')

    xs = np.linspace(0, x_max, 200)
    ax.plot(xs, xs * PG_PER_CELL / 1e6, 'k--', lw=0.5, alpha=0.6,
            label=f'1:1 ({PG_PER_CELL} pg/cell)')

    from matplotlib.ticker import FuncFormatter
    fmt_M = FuncFormatter(lambda x, _: f'{x/1e6:g}')
    ax.xaxis.set_major_formatter(fmt_M)

    ax.set_xlabel('Cells sorted (× 10⁶)', fontsize=6, color='black')
    ax.set_ylabel('gDNA recovered (μg)', fontsize=6, color='black')

    sec = ax.secondary_yaxis('right',
                             functions=(lambda y: y * 1e6 / PG_PER_CELL,
                                        lambda ge: ge * PG_PER_CELL / 1e6))
    sec.set_ylabel('Genome equivalents (× 10⁶)', fontsize=6, color='black')
    sec.yaxis.set_major_formatter(fmt_M)
    sec.tick_params(labelsize=6, colors='black', width=0.5, pad=1, length=2)
    sec.spines['right'].set_linewidth(0.5)

    ax.set_xlim(0, x_max); ax.set_ylim(bottom=0)
    ax.tick_params(labelsize=6, colors='black', width=0.5, pad=1, length=2)
    for spine in ax.spines.values(): spine.set_linewidth(0.5)

    fig.tight_layout(pad=0.3)
    fig.savefig(OUT_DIR / 'cells_vs_ug_dna.pdf')
    fig.savefig(OUT_DIR / 'cells_vs_ug_dna.png', dpi=600)
    plt.close(fig)
print('wrote: cells_vs_ug_dna.{pdf,png}')

# ---------- Compute per-bin variant coverage (one record per bin file) ----------
import re
def display_rep(row):
    """For DAMGO/FSK multidrug controls the spreadsheet rep is the gel-day
    number (1–4) but the filename keeps a per-drug rep that extends to 7.
    Use the filename rep for display so the controls are correctly distinguished."""
    if row['sample'] not in ('DAMGO', 'FSK'):
        return int(row['replicate'])
    fn = row['variant_counts_filename']
    m = re.search(r'_rep_?(\d+)[_.]', fn)
    return int(m.group(1)) if m else int(row['replicate'])

records = []
for _, row in df.iterrows():
    assay = row['assay']
    path = COUNTS_ROOT / assay / row['variant_counts_filename']
    if not path.exists():
        alt = COUNTS_ROOT / 'other_unpublished' / row['variant_counts_filename']
        if alt.exists(): path = alt
        else: print(f'MISSING: {path}'); continue
    counts = pd.read_csv(path, usecols=['count'])['count'].values.astype(np.int64)
    records.append({'sample_id': row['sample_id'], 'assay': assay,
                    'date': int(row['date']), 'replicate': int(row['replicate']),
                    'replicate_display': display_rep(row),
                    'sample': row['sample'], 'conc': row['conc'], 'bin': str(row['bin']),
                    'variant_counts_filename': row['variant_counts_filename'],
                    'mean_coverage': float(counts.mean()),
                    'median_coverage': float(np.median(counts)),
                    'n_variants': int(len(counts)),
                    'counts': counts})

print(f'computed {len(records)} bins')

tbl = pd.DataFrame([{k:v for k,v in r.items() if k!='counts'} for r in records])
tbl = tbl.sort_values(['assay','date','replicate','sample','conc','bin']).reset_index(drop=True)
tbl.to_csv(OUT_DIR / 'per_bin_coverage_summary.csv', index=False)
print('wrote: per_bin_coverage_summary.csv')

# ---------- Plot 2: per-sample mean coverage bars (4 bins summed → 1 bar per sample) ----------
SHORT = {'Mitragynine Pseudoindoxyl':'MP','Surface expression':'Surface'}
def sample_bar_label(r):
    parts = [f"R{r['replicate_display']}", SHORT.get(r['sample'], r['sample'])]
    if pd.notna(r['conc']): parts.append(f"{r['conc']:g}M")
    return ' | '.join(parts)

sample_groups = {}
for r in records:
    sid = (r['date'], r['replicate'], r['sample'], r['conc'] if pd.notna(r['conc']) else None)
    sample_groups.setdefault(sid, []).append(r)

grouped_records = []
for sid, rs in sample_groups.items():
    if not rs: continue
    summed = np.zeros_like(rs[0]['counts'])
    for r in rs: summed = summed + r['counts']
    grouped_records.append({'sample_id': rs[0]['sample_id'], 'assay': rs[0]['assay'],
                            'date': rs[0]['date'], 'replicate': rs[0]['replicate'],
                            'replicate_display': rs[0]['replicate_display'],
                            'sample': rs[0]['sample'], 'conc': rs[0]['conc'],
                            'counts': summed, 'n_bins': len(rs)})

g_sort_key = lambda r: (ASSAY_ORDER.index(r['assay']),
                        r['date'], r['replicate_display'], r['sample'],
                        (r['conc'] if pd.notna(r['conc']) else 99))
grouped_sorted = sorted(grouped_records, key=g_sort_key)
n = len(grouped_sorted)
ncols = 3
chunk = int(np.ceil(n / ncols))
panels = [grouped_sorted[i*chunk:(i+1)*chunk] for i in range(ncols)]
print(f'sample-grouped: {n} bars, {ncols} cols x ~{chunk} bars')

from matplotlib.patches import Patch
handles = [Patch(facecolor=COLORS[a], alpha=0.9, edgecolor='none', label=ASSAY_DISPLAY[a])
           for a in ASSAY_ORDER if any(r['assay']==a for r in grouped_sorted)]

x_max = max(float(np.mean(r['counts'])) for r in grouped_sorted) * 1.1
COL_W_MM = 40
COL_H_MM = 140
mid = ncols // 2
with plt.rc_context(STYLE):
    fig, axes = plt.subplots(1, ncols, figsize=(COL_W_MM*ncols*MM, COL_H_MM*MM))
    if ncols == 1: axes = [axes]
    for i, (ax, panel) in enumerate(zip(axes, panels)):
        positions = np.arange(len(panel))
        means = [float(np.mean(r['counts'])) for r in panel]
        colors = [COLORS[r['assay']] for r in panel]
        labels = [sample_bar_label(r) for r in panel]
        ax.barh(positions, means, color=colors, edgecolor='none', height=0.9)
        ax.set_xscale('log')
        ax.set_xlim(1, x_max)
        if i == mid:
            ax.set_xlabel('Mean coverage (4 bins summed, log)', fontsize=6, color='black')
        ax.set_yticks(positions)
        ax.set_yticklabels(labels, fontsize=6, color='black')
        ax.invert_yaxis()
        ax.set_ylim(len(panel)-0.5, -0.5)
        ax.tick_params(axis='x', labelsize=6, colors='black', width=0.5, pad=2)
        ax.tick_params(axis='y', length=1.5, width=0.3, pad=1)
        for sp in ('top','right'): ax.spines[sp].set_visible(False)
        for sp in ('left','bottom'): ax.spines[sp].set_linewidth(0.5)
    fig.tight_layout(pad=0.3)
    fig.savefig(OUT_DIR / 'alt3_sample_grouped_bars.pdf')
    fig.savefig(OUT_DIR / 'alt3_sample_grouped_bars.png', dpi=400)
    plt.close(fig)
print('wrote: alt3_sample_grouped_bars.{pdf,png}')

print('done')

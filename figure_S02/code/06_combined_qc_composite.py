#!/usr/bin/env python3
"""
Combined QC composite: dose-response + multi-drug saturating effects in one PDF.

Top block (2 rows x 3 cols):
    Row 1: dose-response Score (Morphine | Fentanyl | DAMGO)
    Row 2: dose-response SE    (Morphine | Fentanyl | DAMGO)

Bottom block (2 rows x 1 col):
    Row 3: multi-drug Score (Surface + 16 ligands; MCAM excluded)
    Row 4: multi-drug SE

Left/right edges of all rows align; Score rows share y-axis range/ticks;
SE rows share y-axis range/ticks.

Output: DMS_QC/06_combined_qc_composite.pdf  (160 mm x ~100 mm)
"""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

MM = 1.0 / 25.4

mpl.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 6,
    "text.color": "black",
    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "axes.titlesize": 6,
    "axes.labelsize": 6,
    "xtick.color": "black", "ytick.color": "black",
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
})

HERE = Path(__file__).resolve().parent
CSV = HERE.parent / "dms_scores" / "composite_dms_scores.csv"
df = pd.read_csv(CSV, low_memory=False)

# ===== Dose-response series =====
FENTANYL = [
    (f"MOR_Fentanyl_{v}M", f"{v}", False)
    for v in ["-12", "-11", "-10", "-9", "-8", "-7", "-6", "-5"]
]
MORPHINE = [
    (f"MOR_Morphine_{v}M", f"{float(v)/10:g}", False)
    for v in ["-115", "-105", "-95", "-85", "-75", "-65", "-55", "-45"]
]
DAMGO_MAP = {5: -9, 4: -8, 3: -7, 2: -6, 1: -5}
DAMGO = (
    [("DAMGO_mor_dms_6", "no lig.", True)]
    + [(f"DAMGO_mor_dms_{i}", f"{DAMGO_MAP[i]}", False) for i in [5, 4, 3, 2, 1]]
)
DRC_SERIES = [("Morphine", MORPHINE), ("Fentanyl", FENTANYL), ("DAMGO", DAMGO)]

# ===== Multi-drug series =====
DRUGS = [
    "Surface",
    "FSK", "Naltrexone", "Naloxone", "Nalbuphine", "Buprenorphine",
    "Butorphanol", "MP", "TRV130", "PZM21", "Methadone",
    "C6guano", "Morphine", "SR17018", "Fentanyl", "Carfentanil", "DAMGO",
]
DISPLAY = {"FSK": "No Ligand", "C6guano": "C6-Guano", "SR17018": "SR-17018"}
XLABELS_MD = [DISPLAY.get(d, d) for d in DRUGS]

GREY = "#888888"
GREY_DARK = "#222222"
SUBSAMPLE_N = 800

def fetch_drc(prefix, col):
    a = df[f"{prefix}_rescaled_{col}"].to_numpy()
    return a[~np.isnan(a)]

def fetch_md(drug, col):
    a = df[f"{drug}_{col}"].to_numpy()
    return a[~np.isnan(a)]

def _violin_jitter(d, max_width, rng):
    if len(d) < 2 or np.allclose(d.std(), 0):
        return np.zeros(len(d))
    kde = gaussian_kde(d)
    dens = kde(d)
    norm = dens / dens.max()
    return rng.uniform(-1, 1, size=len(d)) * norm * max_width

def strip_panel(ax, data_list, positions, special_idx=None,
                ylim=None, ylabel=None, jitter=0.38, rng_seed=0,
                point_size=0.7, point_alpha=0.55, connect_medians=False):
    rng = np.random.default_rng(rng_seed)
    for i, (x, d) in enumerate(zip(positions, data_list)):
        if len(d) > SUBSAMPLE_N:
            d_plot = rng.choice(d, size=SUBSAMPLE_N, replace=False)
        else:
            d_plot = d
        c = GREY_DARK if (special_idx is not None and i == special_idx) else GREY
        xs = x + _violin_jitter(d_plot, jitter, rng)
        ax.scatter(xs, d_plot, s=point_size, c=c, alpha=point_alpha,
                   linewidths=0, rasterized=True)
    meds = [np.median(d) for d in data_list]
    ax.plot(positions, meds, color="black", lw=0.8, marker="o", ms=2.4,
            linestyle=("-" if connect_medians else ""), zorder=4)
    for x, d in zip(positions, data_list):
        q25, q75 = np.percentile(d, [25, 75])
        ax.plot([x, x], [q25, q75], color="black", lw=0.7, zorder=3)
    ax.axhline(0, color="black", lw=0.3, ls="--", alpha=0.4)
    if ylim is not None:
        ax.set_ylim(*ylim)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

SCORE_LIM = (-1.25, 0.75)
SCORE_TICKS = [-1.0, -0.5, 0.0, 0.5]
SE_LIM = (0.0, 0.4)
SE_TICKS = [0.0, 0.2, 0.4]

# ===== Figure layout =====
fig = plt.figure(figsize=(120 * MM, 100 * MM))
outer = fig.add_gridspec(
    2, 1, height_ratios=[1, 1], hspace=0.55,
    left=0.10, right=0.99, top=0.96, bottom=0.16,
)
top_gs = outer[0].subgridspec(
    2, 3, width_ratios=[len(s[1]) for s in DRC_SERIES],
    wspace=0.10, hspace=0.20,
)
bot_gs = outer[1].subgridspec(2, 1, hspace=0.20)

# ----- Top block: dose-response -----
for c, (drug, series) in enumerate(DRC_SERIES):
    scores, ses, positions, ticklabels, no_lig_idx = [], [], [], [], None
    for i, (prefix, label, is_nl) in enumerate(series):
        scores.append(fetch_drc(prefix, "score"))
        ses.append(fetch_drc(prefix, "sd"))
        positions.append(i); ticklabels.append(label)
        if is_nl:
            no_lig_idx = i
    xlim = (-0.6, max(positions) + 0.6)

    ax_s = fig.add_subplot(top_gs[0, c])
    strip_panel(ax_s, scores, positions, special_idx=no_lig_idx,
                ylim=SCORE_LIM, ylabel=("Score" if c == 0 else None),
                connect_medians=True)
    ax_s.set_yticks(SCORE_TICKS)
    if c > 0:
        ax_s.set_yticklabels([])
    ax_s.set_xticks(positions); ax_s.tick_params(axis="x", labelbottom=False)
    ax_s.set_xlim(*xlim)
    ax_s.set_title(drug, loc="center", pad=4)

    ax_e = fig.add_subplot(top_gs[1, c], sharex=ax_s)
    strip_panel(ax_e, ses, positions, special_idx=no_lig_idx,
                ylim=SE_LIM, ylabel=("SE" if c == 0 else None))
    ax_e.set_yticks(SE_TICKS)
    if c > 0:
        ax_e.set_yticklabels([])
    ax_e.set_xticks(positions); ax_e.set_xticklabels(ticklabels)
    ax_e.set_xlim(*xlim)
    ax_e.set_xlabel("log[ligand] (M)")

# ----- Bottom block: multi-drug -----
md_positions = list(range(len(DRUGS)))
md_scores = [fetch_md(d, "effect") for d in DRUGS]
md_ses    = [fetch_md(d, "effect_se") for d in DRUGS]
surf_idx = 0

ax_md_s = fig.add_subplot(bot_gs[0, 0])
strip_panel(ax_md_s, md_scores, md_positions, special_idx=surf_idx,
            ylim=SCORE_LIM, ylabel="Score")
ax_md_s.set_yticks(SCORE_TICKS)
ax_md_s.set_xticks(md_positions); ax_md_s.tick_params(axis="x", labelbottom=False)
ax_md_s.set_xlim(-0.6, len(DRUGS) - 0.4)

ax_md_e = fig.add_subplot(bot_gs[1, 0], sharex=ax_md_s)
strip_panel(ax_md_e, md_ses, md_positions, special_idx=surf_idx,
            ylim=SE_LIM, ylabel="SE")
ax_md_e.set_yticks(SE_TICKS)
ax_md_e.set_xticks(md_positions)
ax_md_e.set_xticklabels(XLABELS_MD, rotation=45, ha="right",
                        rotation_mode="anchor")
ax_md_e.set_xlim(-0.6, len(DRUGS) - 0.4)

# Dashed vertical line separating Surface from ligand panel (multi-drug block)
for a in (ax_md_s, ax_md_e):
    a.axvline(0.5, color="black", lw=0.4, ls="--", alpha=0.7, zorder=0)

out = HERE / "06_combined_qc_composite.pdf"
fig.savefig(out, dpi=600)
print(f"Wrote {out}")

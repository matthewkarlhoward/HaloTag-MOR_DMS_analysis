#!/usr/bin/env python3
"""
Multi-drug saturating-effect QC composite.

Row 1: effect score per drug (jittered point cloud, KDE-shaped jitter)
Row 2: effect SE per drug (same)
All 17 conditions on one axis (Surface + 16 drugs; MCAM excluded).

Output: DMS_QC/05_multidrug_composite.pdf  (160 mm x 50 mm)
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

# Surface first (baseline), then ligands ordered antagonist -> partial -> full agonist.
# MCAM excluded. FSK is the no-ligand control.
DRUGS = [
    "Surface",
    "FSK", "Naltrexone", "Naloxone", "Nalbuphine", "Buprenorphine",
    "Butorphanol", "MP", "TRV130", "PZM21", "Methadone",
    "C6guano", "Morphine", "SR17018", "Fentanyl", "Carfentanil", "DAMGO",
]
# Display labels for the x-axis (column name -> pretty label)
DISPLAY = {
    "FSK":     "No Ligand",
    "C6guano": "C6-Guano",
    "SR17018": "SR-17018",
}
XLABELS = [DISPLAY.get(d, d) for d in DRUGS]

GREY = "#888888"
GREY_DARK = "#222222"
SURFACE_COL = GREY_DARK  # darker grey to mark the Surface baseline

SUBSAMPLE_N = 800

def fetch(drug, col):
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
                point_size=0.7, point_alpha=0.55):
    rng = np.random.default_rng(rng_seed)
    for i, (x, d) in enumerate(zip(positions, data_list)):
        if len(d) > SUBSAMPLE_N:
            d_plot = rng.choice(d, size=SUBSAMPLE_N, replace=False)
        else:
            d_plot = d
        c = SURFACE_COL if (special_idx is not None and i == special_idx) else GREY
        xs = x + _violin_jitter(d_plot, jitter, rng)
        ax.scatter(xs, d_plot, s=point_size, c=c, alpha=point_alpha,
                   linewidths=0, rasterized=True)
    meds = [np.median(d) for d in data_list]
    ax.plot(positions, meds, color="black", lw=0.8, marker="o", ms=2.4,
            linestyle="", zorder=4)
    for x, d in zip(positions, data_list):
        q25, q75 = np.percentile(d, [25, 75])
        ax.plot([x, x], [q25, q75], color="black", lw=0.7, zorder=3)
    ax.axhline(0, color="black", lw=0.3, ls="--", alpha=0.4)
    if ylim is not None:
        ax.set_ylim(*ylim)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

# Pull data
scores = [fetch(d, "effect") for d in DRUGS]
ses    = [fetch(d, "effect_se") for d in DRUGS]
positions = list(range(len(DRUGS)))
surf_idx = 0  # Surface is the first column

# ---- Layout 160 x 50 mm ----
fig = plt.figure(figsize=(160 * MM, 50 * MM))
gs = fig.add_gridspec(
    2, 1, height_ratios=[1, 1],
    hspace=0.20,
    left=0.06, right=0.99, top=0.96, bottom=0.30,
)

SCORE_LIM = (-1.25, 0.75)
SCORE_TICKS = [-1.0, -0.5, 0.0, 0.5]
SE_LIM = (0.0, 0.4)

ax_top = fig.add_subplot(gs[0, 0])
strip_panel(ax_top, scores, positions, special_idx=surf_idx,
            ylim=SCORE_LIM, ylabel="Score")
ax_top.set_yticks(SCORE_TICKS)
ax_top.set_xticks(positions)
ax_top.tick_params(axis="x", labelbottom=False)
ax_top.set_xlim(-0.6, len(DRUGS) - 0.4)

ax_bot = fig.add_subplot(gs[1, 0], sharex=ax_top)
strip_panel(ax_bot, ses, positions, special_idx=surf_idx,
            ylim=SE_LIM, ylabel="SE")
ax_bot.set_xticks(positions)
ax_bot.set_xticklabels(XLABELS, rotation=45, ha="right",
                       rotation_mode="anchor")
ax_bot.set_xlim(-0.6, len(DRUGS) - 0.4)

# Dashed vertical line separating Surface (column 0) from ligand conditions
for a in (ax_top, ax_bot):
    a.axvline(0.5, color="black", lw=0.4, ls="--", alpha=0.7, zorder=0)

out = HERE / "05_multidrug_composite.pdf"
fig.savefig(out, dpi=600)
print(f"Wrote {out}")

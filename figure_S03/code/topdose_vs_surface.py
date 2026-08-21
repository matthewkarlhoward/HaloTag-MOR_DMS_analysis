#!/usr/bin/env python3
"""
Top-dose (highest concentration) fitness score vs surface-expression score,
one column per ligand.

Top dose per ligand (rescaled single-concentration score):
    Morphine  MOR_Morphine_-45M_rescaled_score   (-4.5 log[M], ~31.6 uM top)
    Fentanyl  MOR_Fentanyl_-5M_rescaled_score     (-5 log[M], 10 uM)
    DAMGO     DAMGO_mor_dms_1_rescaled_score       (index 1 = -5 log[M], 10 uM)

Rows:    variant-level (top) and position-level mean, n>=5 (bottom)
x-axis:  Surface_effect (surface-expression score)

Filter:  missense only.  Spearman r (n) annotated per panel.

Output: plots/scatter/ec50_emax/topdose_vs_surface.{pdf,png}
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[3]
OUT  = ROOT / "plots" / "scatter" / "ec50_emax"

TOPDOSE = {
    "Morphine": ("MOR_Morphine_-45M_rescaled_score", "31.6 uM"),
    "Fentanyl": ("MOR_Fentanyl_-5M_rescaled_score",  "10 uM"),
    "DAMGO":    ("DAMGO_mor_dms_1_rescaled_score",    "10 uM"),
}
LIGANDS = ["Morphine", "Fentanyl", "DAMGO"]
MIN_VARIANTS = 5

plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 6,
    "axes.labelsize": 6, "axes.titlesize": 6,
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "text.color": "black",
    "axes.edgecolor": "black", "axes.labelcolor": "black",
    "axes.titlecolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "lines.linewidth": 0.5, "patch.linewidth": 0.5,
})

PT_COL = "#3b3b3b"

dms = pd.read_csv(ROOT / "dms_scores/composite_dms_scores.csv", low_memory=False)
mis = dms[dms.type == "missense"].copy()

def annotate_r(ax, x, y):
    ok = x.notna() & y.notna()
    if ok.sum() < 3:
        return
    rho, _ = spearmanr(x[ok], y[ok])
    ax.text(0.04, 0.95, f"$r_s$ = {rho:+.2f}\nn = {ok.sum():,}",
            transform=ax.transAxes, va="top", ha="left", fontsize=6,
            color="black")

MM = 1 / 25.4
fig, axes = plt.subplots(2, 3, figsize=(165 * MM, 110 * MM),
                         sharex=True, sharey=True)

for j, drug in enumerate(LIGANDS):
    col, dose = TOPDOSE[drug]

    # Row 0: variant-level
    ax = axes[0, j]
    ax.scatter(mis.Surface_effect, mis[col], s=3, color=PT_COL,
               edgecolors="none", alpha=0.35, zorder=2, rasterized=True)
    annotate_r(ax, mis.Surface_effect, mis[col])
    ax.set_title(f"{drug} ({dose})", color="black")
    if j == 0:
        ax.set_ylabel("top-dose score\n(per variant)", color="black")

    # Row 1: position-level mean (n>=5)
    pos = (mis.groupby("position")
              .agg(score=(col, "mean"), surface=("Surface_effect", "mean"),
                   n=(col, "count"))
              .query("n >= @MIN_VARIANTS").reset_index())
    ax = axes[1, j]
    ax.scatter(pos.surface, pos.score, s=6, color=PT_COL,
               edgecolors="none", alpha=0.55, zorder=2)
    annotate_r(ax, pos.surface, pos.score)
    if j == 0:
        ax.set_ylabel("mean top-dose score\n(per position)", color="black")
    ax.set_xlabel("surface expression score", color="black")

for ax in axes.ravel():
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.5)

fig.tight_layout()
out_pdf = OUT / "topdose_vs_surface.pdf"
out_png = OUT / "topdose_vs_surface.png"
fig.savefig(out_pdf, dpi=600)
fig.savefig(out_png, dpi=600)
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")
plt.close()

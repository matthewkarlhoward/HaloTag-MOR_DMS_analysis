#!/usr/bin/env python3
"""
Variant-level Morphine vs Fentanyl Emax scatter (40 x 40 mm), Activity-scaled.

One point = one missense variant whose top plateau (the emax parameter) is
anchored in BOTH ligands -- curve_type sigmoid, flat, or no_baseline.
  - Source: robust 3-parameter refit (refit_3param_robust_{drug}.csv)
  - Emax is put on the Activity scale per drug: that drug's sigmoid-synonymous
    mean Emin -> 0 and mean Emax -> 1 (1 = WT-like efficacy, 0 = dead baseline).
  - x = Fentanyl Emax (Activity),  y = Morphine Emax (Activity)

Points are colored by local 2-D density (gaussian KDE, magma, dense-on-top); a
grey dotted y = x identity line is drawn for reference. Both axes tick every 0.5
and label only the whole numbers.

Split out from variant_morphine_vs_fentanyl.py (the EC50 half is its sibling
variant_morphine_vs_fentanyl_ec50.py).

Output: variant_morphine_vs_fentanyl_emax.pdf
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter
from scipy.stats import gaussian_kde
from pathlib import Path

HERE = Path(__file__).resolve().parent

plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 6,
    "mathtext.default": "regular",
    "axes.labelsize": 6, "axes.titlesize": 6,
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "text.color": "black", "axes.edgecolor": "black",
    "axes.labelcolor": "black", "axes.titlecolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "lines.linewidth": 0.5, "patch.linewidth": 0.5,
})
MM = 1 / 25.4

def setup(ax, xlab, ylab):
    ax.set_xlabel(xlab, labelpad=1)
    ax.set_ylabel(ylab, labelpad=1)
    ax.tick_params(axis="both", pad=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.5)

def diag(ax, x, y):
    lo = float(min(x.min(), y.min()))
    hi = float(max(x.max(), y.max()))
    ax.plot([lo, hi], [lo, hi], color="grey", lw=0.5, ls=":", zorder=1)

def scatter_panel(x, y, xlab, ylab, out_pdf, tick_step=None):
    x, y = np.asarray(x, float), np.asarray(y, float)
    fig, ax = plt.subplots(figsize=(40 * MM, 40 * MM))
    diag(ax, x, y)
    # Color by local 2-D density (gaussian KDE); draw dense points on top.
    z = gaussian_kde(np.vstack([x, y]))(np.vstack([x, y]))
    o = z.argsort()
    ax.scatter(x[o], y[o], c=z[o], s=1, cmap="magma",
               edgecolors="none", alpha=0.85, zorder=2)
    setup(ax, xlab, ylab)
    if tick_step:
        # Tick every tick_step on both axes; label only every other tick (at
        # multiples of 2*tick_step) -> whole numbers for Emax (step 0.5).
        lab = 2 * tick_step
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_locator(MultipleLocator(tick_step))
            axis.set_major_formatter(FuncFormatter(
                lambda v, _, s=lab: f"{int(round(v))}".replace("-", "−")
                if abs(v / s - round(v / s)) < 1e-9 else ""))
    fig.subplots_adjust(left=0.205, bottom=0.175, right=0.97, top=0.98)
    fig.savefig(out_pdf, dpi=600)
    plt.close(fig)
    print(f"Saved: {out_pdf}  (n = {len(x)})")

# ── Emax: robust refit missense in both ligands, Activity-scaled ──────────────
mor = pd.read_csv(HERE / "refit_3param_robust_morphine.csv")
fen = pd.read_csv(HERE / "refit_3param_robust_fentanyl.csv")
# Keep curve types whose TOP plateau (the emax parameter) is anchored: sigmoid,
# flat, and no_baseline (only its bottom baseline is unreached). Adding
# "no_plateau" would gain ~1000 points but its emax is extrapolated; "boundary"
# fits are pinned at the parameter bound and stay excluded.
EMAX_OK = {"sigmoid", "flat", "no_baseline"}

def _activity_miss(df, prefix):
    """Missense Emax on the Activity scale (sigmoid-syn Emin -> 0, Emax -> 1)."""
    syn = df[(df.type == "synonymous") & (df.new_curve_type == "sigmoid")]
    e0, e1 = syn.fitted_emin.mean(), syn.fitted_emax.mean()
    sub = df[df.fit_ok & (df.type == "missense") & df.new_curve_type.isin(EMAX_OK)].copy()
    sub[f"{prefix}_emax"] = (sub.fitted_emax - e0) / (e1 - e0)
    return sub[["hgvs", f"{prefix}_emax"]]

emax = _activity_miss(mor, "Mor").merge(_activity_miss(fen, "Fen"), on="hgvs")
scatter_panel(emax.Fen_emax, emax.Mor_emax,
              r"E$_{max}$ Fentanyl", r"E$_{max}$ Morphine",
              HERE / "variant_morphine_vs_fentanyl_emax.pdf", tick_step=0.5)

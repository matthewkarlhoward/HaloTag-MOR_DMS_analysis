#!/usr/bin/env python3
"""
Variant-level Morphine vs Fentanyl EC50 scatter with Q126 highlighted (40 x 40 mm).

Same panel and same filtering as variant_morphine_vs_fentanyl_ec50.py -- one
point = one missense variant with a clean sigmoid fit in BOTH ligands -- with
the position-126 variants drawn on top in green so their location in the cloud
is obvious.

Q126 sits well to the RIGHT of the identity line: fentanyl EC50 is shifted ~1.6
log units weaker (median -6.52 vs -8.10 overall) while morphine barely moves
(-7.70 vs -8.00), i.e. fentanyl-selective potency loss.

Highlighted points carry a black outline over a thin white halo, which keeps
them separable where they sit on top of the dark background cloud.

Two outputs:
  variant_morphine_vs_fentanyl_ec50_q126.pdf        single "Q126" group label
  variant_morphine_vs_fentanyl_ec50_q126_plain.pdf  points only, label in caption

Backend is forced to Agg and the size set explicitly: the macosx backend rounds
figsize to 2 dp, which breaks the exact 40 mm panel width.
"""
import matplotlib
matplotlib.use("Agg")

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
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
    "pdf.fonttype": 42,
})
MM = 1 / 25.4

POSITION = 126
# Highlight colour must clear two hurdles: it cannot collide with magma (black /
# purple / orange / yellow) or the overlay vanishes into the density cloud, and
# it cannot be blue or red, which mean "morphine bias" / "fentanyl bias" in the
# neighbouring panels. That leaves the green-teal wedge. Override from argv[1].
HILITE = sys.argv[1] if len(sys.argv) > 1 else "#00A651"
OUT_SUFFIX = sys.argv[2] if len(sys.argv) > 2 else ""

# Halo profile: (stroke width pt, alpha), widest/faintest drawn first. A handful
# of steps bands visibly under magnification, so the default ramps many thin
# strokes at low alpha instead -- cumulative opacity 1-(1-a)^k rises smoothly
# toward the marker and the rings disappear.
def _ramp(width, n=20, alpha=0.13):
    return [(w, alpha) for w in np.linspace(width, 0.2, n)]


# Widths are stroke diameters centred on the marker outline, so the halo reaches
# width/2 beyond the point: the 1.75 pt default radiates ~0.31 mm at 40 mm.
HALO_PROFILES = {
    "hard":   [(1.2, 1.00)],
    "soft":   [(2.6, 0.12), (2.1, 0.20), (1.6, 0.34), (1.1, 0.60)],
    "smooth": _ramp(1.75),
    "smooth_wide": _ramp(2.4, n=26, alpha=0.11),
}
HALO = HALO_PROFILES[os.environ.get("HALO_PROFILE", "smooth")]


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


def scatter_panel(df, hit, xlab, ylab, out_pdf, tick_step=1, annotate=True):
    x, y = df.Fen.to_numpy(float), df.Mor.to_numpy(float)

    fig, ax = plt.subplots()
    fig.set_size_inches(40 * MM, 40 * MM)
    diag(ax, x, y)

    # Background: color by local 2-D density (gaussian KDE), dense points on top.
    z = gaussian_kde(np.vstack([x, y]))(np.vstack([x, y]))
    o = z.argsort()
    ax.scatter(x[o], y[o], c=z[o], s=1, cmap="magma",
               edgecolors="none", alpha=0.85, zorder=2)

    # Overlay: every variant at the highlighted position. Behind each marker sit
    # nested white strokes, widest and faintest first, so the halo fades out
    # instead of ending on a hard rim. A true blur would need an agg filter,
    # which rasterises; stacked strokes keep the panel vector.
    ov = ax.scatter(hit.Fen, hit.Mor, s=4, facecolors=HILITE,
                    edgecolors="black", linewidths=0.25, zorder=3)
    ov.set_path_effects(
        [pe.withStroke(linewidth=w, foreground="white", alpha=a)
         for w, a in HALO] + [pe.Normal()])

    setup(ax, xlab, ylab)
    # Tick every tick_step on both axes; label only every other tick.
    lab = 2 * tick_step
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(MultipleLocator(tick_step))
        axis.set_major_formatter(FuncFormatter(
            lambda v, _, s=lab: f"{int(round(v))}".replace("-", "−")
            if abs(v / s - round(v / s)) < 1e-9 else ""))

    wt = hit.wildtype.iloc[0].upper()
    if annotate:
        # Short leader from the group centroid up-right into open space; a label
        # parked in a far corner would drag its line across the whole cloud.
        ax.annotate(f"{wt}{POSITION}",
                    xy=(hit.Fen.median(), hit.Mor.median()),
                    xytext=(-6.05, -6.30), textcoords="data",
                    size=6, color=HILITE, ha="left", va="bottom", zorder=4,
                    arrowprops=dict(arrowstyle="-", color=HILITE, lw=0.5,
                                    shrinkA=1, shrinkB=2))

    fig.subplots_adjust(left=0.205, bottom=0.175, right=0.97, top=0.98)
    fig.savefig(out_pdf, dpi=600)
    plt.close(fig)
    print(f"Saved: {out_pdf.name}  (n = {len(x)}, {wt}{POSITION} n = {len(hit)})")


# ── EC50: robust refit, sigmoid missense in both ligands ─────────────────────
mor = pd.read_csv(HERE / "refit_3param_robust_morphine.csv")
fen = pd.read_csv(HERE / "refit_3param_robust_fentanyl.csv")


def _refit_miss(df, prefix, keep_meta=False):
    sub = df[df.fit_ok & (df.type == "missense") &
             (df.new_curve_type == "sigmoid")]
    cols = ["hgvs", "fitted_ec50_logM"]
    if keep_meta:
        cols += ["position", "wildtype", "mutation"]
    return sub[cols].rename(columns={"fitted_ec50_logM": prefix})


ec50 = _refit_miss(mor, "Mor", keep_meta=True).merge(
    _refit_miss(fen, "Fen"), on="hgvs")
hit = ec50[ec50.position == POSITION]
if hit.empty:
    raise SystemExit(f"No variants pass the EC50 filter at position {POSITION}")

for annotate in (True, False):
    stem = "variant_morphine_vs_fentanyl_ec50_q126"
    scatter_panel(ec50, hit,
                  r"EC$_{50}$ Fentanyl", r"EC$_{50}$ Morphine",
                  HERE / f"{stem}{'' if annotate else '_plain'}{OUT_SUFFIX}.pdf",
                  tick_step=1, annotate=annotate)

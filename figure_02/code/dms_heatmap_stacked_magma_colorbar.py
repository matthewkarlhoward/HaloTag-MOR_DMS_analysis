#!/usr/bin/env python3
"""
Standalone colour scale for dms_heatmap_stacked_magma.pdf, as a stripped-down
PDF to drop into Illustrator.

Same ramp and same limits as the "range" scaling in dms_stacked_heatmap_magma.R:
reversed magma_fastyellow (LOF = bright yellow, GOF = black) over the empirical
2-98% range of the plotted effects, recomputed here from the same rows the
heatmap draws (20 substitutions, positions 65-355, all 16 conditions).

The bar is drawn as 256 abutting filled rectangles, not an image or a smooth
gradient, so it stays fully vector and stays editable after Illustrator import.

Outputs:
  dms_heatmap_stacked_magma_colorbar_h.pdf   horizontal, 6.5 x 2 mm bar
  dms_heatmap_stacked_magma_colorbar_v.pdf   vertical,   2 x 6.5 mm bar
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")                      # exact figure size (see memory: macosx snaps figsize)
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[2]
OUT  = ROOT / "plots" / "heatmaps"
MM   = 1 / 25.4

plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 6,
    "text.color": "black", "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "xtick.labelcolor": "black", "ytick.labelcolor": "black",
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "axes.linewidth": 0.5, "lines.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    # TrueType, not the matplotlib default Type 3: keeps the labels as editable
    # Helvetica text in Illustrator instead of unstylable outlines
    "pdf.fonttype": 42,
})

TXT      = 6
LW       = 0.5      # pt
BAR_LONG = 6.5      # mm
BAR_THIN = 2.0      # mm
NSTEP    = 256      # rectangles across the bar

# ---- ramp: magma resampled at t^0.55, then reversed (LOF -> yellow) ----------
_magma = plt.get_cmap("magma")
CMAP   = LinearSegmentedColormap.from_list(
    "magma_fastyellow", _magma(np.linspace(0, 1, 256) ** 0.55)).reversed()

# ---- limits: 2-98% of exactly the effects the heatmap plots -----------------
CONDITIONS = ["FSK", "Naltrexone", "Naloxone", "Nalbuphine", "Buprenorphine",
              "Butorphanol", "MP", "TRV130", "PZM21", "Methadone", "C6guano",
              "Fentanyl", "Carfentanil", "Morphine", "SR17018", "DAMGO"]
AA = list("ACDEFGHIKLMNPQRSTVWY")

dms = pd.read_csv(ROOT / "dms_scores" / "composite_dms_scores.csv")
sub = dms[dms["mutation"].isin(AA) & dms["position"].between(65, 355)]
vals = np.concatenate([sub[f"{c}_effect"].to_numpy() for c in CONDITIONS])
vals = vals[~np.isnan(vals)]
LO, HI = (float(np.quantile(vals, 0.02)), float(np.quantile(vals, 0.98)))
# Ticks are labelled with what they mean, not with numbers: the two ends of the
# window and the synonymous baseline at 0.
TICKS  = [LO, 0.0, HI]
LABELS = ["LOF", "Neutral", "GOF"]


def draw(vertical):
    """One bar + frame + ticks, sized so the bar is exactly BAR_LONG x BAR_THIN."""
    # padding leaves room for the tick labels and the LOF/GOF anchors, which
    # sit outside the axes; too little and they clip at the canvas edge.
    pad_l, pad_r = (1.5, 8.0) if vertical else (7.0, 7.0)
    pad_b, pad_t = (2.5, 2.5) if vertical else (4.0, 2.0)
    bw, bh = (BAR_THIN, BAR_LONG) if vertical else (BAR_LONG, BAR_THIN)

    fig = plt.figure(figsize=((bw + pad_l + pad_r) * MM, (bh + pad_b + pad_t) * MM))
    fig.patch.set_alpha(0)                       # clear background, no white box
    ax = fig.add_axes([pad_l / (bw + pad_l + pad_r), pad_b / (bh + pad_b + pad_t),
                       bw / (bw + pad_l + pad_r), bh / (bh + pad_b + pad_t)])
    ax.patch.set_alpha(0)

    # bar as NSTEP abutting rectangles (vector, no image, no gradient mesh)
    edges = np.linspace(0.0, 1.0, NSTEP + 1)
    for lo_e, hi_e in zip(edges[:-1], edges[1:]):
        c = CMAP((lo_e + hi_e) / 2)
        xy, w, h = ((0.0, lo_e), 1.0, hi_e - lo_e) if vertical \
              else ((lo_e, 0.0), hi_e - lo_e, 1.0)
        ax.add_patch(Rectangle(xy, w, h, facecolor=c, edgecolor="none",
                               linewidth=0, antialiased=False))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    tpos = [(t - LO) / (HI - LO) for t in TICKS]

    if vertical:
        # ticks and labels on the RIGHT of the bar
        ax.set_xticks([])
        ax.set_yticks(tpos)
        ax.set_yticklabels(LABELS)
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.tick_params(axis="y", length=2, width=LW, pad=1, labelsize=TXT)
    else:
        # "Neutral" is the only label that fits under the bar; the two ends are
        # set outside the bar instead, where they cannot collide with it
        ax.set_yticks([])
        ax.set_xticks([tpos[1]])
        ax.set_xticklabels([LABELS[1]])
        ax.tick_params(axis="x", length=2, width=LW, pad=1, labelsize=TXT)
        # offset in points, not axes fraction: a fraction of a 6.5 mm bar is
        # too small a gap and the words touch the frame
        ax.annotate("LOF", (0, 0.5), xycoords="axes fraction", xytext=(-3, 0),
                    textcoords="offset points", ha="right", va="center", fontsize=TXT)
        ax.annotate("GOF", (1, 0.5), xycoords="axes fraction", xytext=(3, 0),
                    textcoords="offset points", ha="left", va="center", fontsize=TXT)

    for s in ax.spines.values():
        s.set_linewidth(LW)
        s.set_visible(True)
    return fig


for vertical, tag in ((False, "h"), (True, "v")):
    f = OUT / f"dms_heatmap_stacked_magma_colorbar_{tag}.pdf"
    fig = draw(vertical)
    # transparent=True: no white canvas rectangle to delete in Illustrator
    fig.savefig(f, format="pdf", transparent=True)   # no bbox_inches: exact mm size
    plt.close(fig)
    print(f"Saved: {f.name}   limits [{LO:.3f}, {HI:.3f}]")

#!/usr/bin/env python3
"""
Efficacy-order reversal shown as dose-response curves: WT vs an A119L double.

Two panels, four ligands each (DAMGO / PZM21 / Nalbuphine / Naloxone), coloured
by efficacy class. Signal-down BRET (baseline 1, deeper = more activation), so
curve depth = efficacy.

  WT           : DAMGO > PZM21 > Nalbuphine > Naloxone  (textbook order)
  A119L+V175N  : DAMGO > Nalbuphine > Naloxone > PZM21  (PZM21 falls to last;
                 naloxone, a WT antagonist, is now a mid-rank agonist)

Curves are the F-test-gated fits with mean +/- SEM points.
Source: figures/doubles_merged_points.csv
Outputs: plots/doubles/rank_example_curves.pdf / .png
"""
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from drc_fit import sigmoid, fit_drc

mpl.rcParams.update({
    "font.family": "Helvetica", "font.weight": "normal", "font.size": 6,
    "axes.titlesize": 6, "axes.labelsize": 6, "xtick.labelsize": 6,
    "ytick.labelsize": 6, "axes.titleweight": "normal",
    "axes.labelweight": "normal", "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "text.color": "black",
    "axes.edgecolor": "black", "axes.labelcolor": "black",
    "xtick.color": "black", "ytick.color": "black", "mathtext.default": "regular",
})

MM = 1 / 25.4
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "plots" / "doubles"

LIG_COL = {"DAMGO": "#1a1a1a", "PZM21": "#c0392b",
           "Nalbuphine": "#2e6db4", "Naloxone": "#9159a8"}
# 2x2 mutant cycle (row-major): WT, A119L / K100N, A119L+K100N
PANELS = [("WT", "WT"), ("119L", "A119L"),
          ("100N", "K100N"), ("119L_100N", "A119L + K100N")]
XMIN, XMAX = -12.0, -4.0            # data clipped at -12 (drops the -13 baseline)
XGRID = np.linspace(XMIN, XMAX, 300)
XPAD = 0.4                          # axis padding so edge points aren't clipped


def build():
    df = pd.read_csv(ROOT / "figures" / "doubles_merged_points.csv")
    df = df[df.logM >= XMIN]                             # clip at -12
    fig, axes = plt.subplots(2, 2, sharex=False, sharey=True,
                             figsize=(91 * MM, 97 * MM))   # ~square panels
    fig.subplots_adjust(left=0.12, right=0.99, top=0.965, bottom=0.085,
                        wspace=0.10, hspace=0.30)
    for ax in axes.flat:
        ax.set_box_aspect(1)                              # force square boxes

    for idx, (tok, title) in enumerate(PANELS):
        ax = axes.flat[idx]
        row, col_i = idx // 2, idx % 2
        for lig, col in LIG_COL.items():
            s = df[(df.ligand == lig) & (df.variant == f"{lig}_{tok}")]
            if tok in ("WT", "119L"):                    # drop the -4 point here
                s = s[s.logM <= -4.5]
            if s.empty:
                continue
            f = fit_drc(s.logM.to_numpy(), s.response.to_numpy())
            xlo, xhi = float(s.logM.min()), float(s.logM.max())   # trim to data
            if f["ok"]:
                if f["responsive"]:
                    g = np.linspace(xlo, xhi, 300)
                    ax.plot(g, sigmoid(g, *f["params"]), color=col,
                            lw=1.0, zorder=3)
                else:
                    ax.plot([xlo, xhi], [f["level"]] * 2, color=col,
                            lw=1.0, zorder=2)
            p = s.groupby("logM")["response"].agg(["mean", "sem"]).reset_index()
            ax.errorbar(p["logM"], p["mean"], yerr=p["sem"], fmt="o", ms=1.8,
                        mfc=col, mec=col, ecolor=col, elinewidth=0.5, capsize=0,
                        lw=0, zorder=4, alpha=0.9)
        ax.axhline(1.0, color="black", lw=0.4, ls=(0, (2, 2)), zorder=1)
        ax.set_xlim(XMIN - XPAD, XMAX + XPAD)
        ax.set_ylim(0.5, 1.12)
        ax.set_yticks([0.5, 0.75, 1.0])
        xt = np.arange(XMIN, XMAX + 1, 1)
        ax.set_xticks(xt)
        ax.set_xticklabels([f"{int(x)}" if x % 2 == 0 else "" for x in xt])
        ax.set_title(title, pad=3)
        ax.set_xlabel("log[ligand]", labelpad=2)         # on every panel
        if col_i == 0:
            ax.set_ylabel("Norm. BRET", labelpad=2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    handles = [Line2D([], [], color=c, lw=1.0, label=l) for l, c in LIG_COL.items()]
    axes.flat[0].legend(handles=handles, loc="lower left", frameon=False,
                        fontsize=5, handlelength=1.1, labelspacing=0.28,
                        handletextpad=0.4, borderaxespad=0.3)

    fig.savefig(OUT / "rank_example_curves.pdf")
    fig.savefig(OUT / "rank_example_curves.png", dpi=300)
    plt.close(fig)
    print("wrote rank_example_curves.pdf / .png")


if __name__ == "__main__":
    build()

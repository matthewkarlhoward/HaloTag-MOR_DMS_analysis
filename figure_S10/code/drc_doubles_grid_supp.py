#!/usr/bin/env python3
"""
Supplement: every double-mutant dose-response curve on one page.

6 partner rows (K100D, K100N, V175N, R278D, I280A, I280K) x 4 ligand columns.
Each panel overlays WT / A119L / partner single / A119L+partner double, as
mean +/- SEM points with the F-test-gated fit (drc_fit.fit_drc): a significant
curve is drawn as a solid sigmoid, a non-response as a flat dashed line.

Colours match Fig 8: WT black, A119L blue, single red, double purple.
House style: black text size 6 Helvetica, 0.5pt lines. Exact 165 x 180 mm.

Source: figures/doubles_merged_points.csv
Outputs: plots/doubles/supp_doubles_drc_grid.pdf / .png
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

LIGANDS = ["DAMGO", "PZM21", "Nalbuphine", "Naloxone"]
PARTNERS = [("100D", "K100D"), ("100N", "K100N"), ("175N", "V175N"),
            ("278D", "R278D"), ("280A", "I280A"), ("280K", "I280K")]

C_WT, C_A119L, C_SINGLE, C_DOUBLE = "#1a1a1a", "#2e6db4", "#c0392b", "#9159a8"
XMIN, XMAX = -13.0, -4.0
XGRID = np.linspace(XMIN, XMAX, 300)


def build():
    df = pd.read_csv(ROOT / "figures" / "doubles_merged_points.csv")

    def sub(lig, tok):
        return df[(df.ligand == lig) & (df.variant == f"{lig}_{tok}")]

    nrow, ncol = len(PARTNERS), len(LIGANDS)
    fig, axes = plt.subplots(nrow, ncol, sharex=True, sharey=True,
                             figsize=(165 * MM, 140 * MM))
    fig.subplots_adjust(left=0.095, right=0.985, top=0.945, bottom=0.075,
                        wspace=0.20, hspace=0.16)

    for i, (tok, plabel) in enumerate(PARTNERS):
        series = [("WT", C_WT), ("119L", C_A119L),
                  (tok, C_SINGLE), (f"119L_{tok}", C_DOUBLE)]
        for j, lig in enumerate(LIGANDS):
            ax = axes[i, j]
            for var, col in series:
                s = sub(lig, var)
                if s.empty:
                    continue
                f = fit_drc(s.logM.to_numpy(), s.response.to_numpy())
                if f["ok"]:
                    if f["responsive"]:
                        ax.plot(XGRID, sigmoid(XGRID, *f["params"]), color=col,
                                lw=0.9, ls="-", zorder=3)
                    else:                              # no response -> flat, solid
                        ax.plot([XMIN, XMAX], [f["level"]] * 2, color=col,
                                lw=0.9, ls="-", zorder=2)
                prof = s.groupby("logM")["response"].agg(["mean", "sem"]).reset_index()
                ax.errorbar(prof["logM"], prof["mean"], yerr=prof["sem"], fmt="o",
                            ms=1.4, mfc=col, mec=col, ecolor=col, elinewidth=0.5,
                            capsize=0, lw=0, zorder=4, alpha=0.9)
            ax.axhline(1.0, color="black", lw=0.4, ls=(0, (2, 2)), zorder=1)
            ax.set_xlim(XMIN, XMAX)
            ax.set_ylim(0.5, 1.2)
            ax.set_yticks([0.5, 0.75, 1.0])
            xt = np.arange(XMIN, XMAX + 1, 1)
            ax.set_xticks(xt)
            ax.set_xticklabels([f"{int(x)}" if x % 2 == 0 else "" for x in xt])
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if i == 0:
                ax.set_title(lig, pad=3)
            if j == 0:
                ax.set_ylabel(f"{plabel}\nNorm. BRET", labelpad=2)
            if i == nrow - 1:
                ax.set_xlabel("log[agonist] (M)", labelpad=2)

    handles = [
        Line2D([], [], color=C_WT, lw=0.9, label="WT"),
        Line2D([], [], color=C_A119L, lw=0.9, label="A119L"),
        Line2D([], [], color=C_SINGLE, lw=0.9, label="partner single"),
        Line2D([], [], color=C_DOUBLE, lw=0.9, label="A119L + partner"),
        Line2D([], [], color="black", marker="o", ms=1.4, lw=0, label="mean ± SEM"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.995),
               frameon=False, fontsize=6, ncol=5, handlelength=1.3,
               columnspacing=1.4, handletextpad=0.4)

    fig.savefig(OUT / "supp_doubles_drc_grid.pdf")
    fig.savefig(OUT / "supp_doubles_drc_grid.png", dpi=300)
    plt.close(fig)
    print("wrote supp_doubles_drc_grid.pdf / .png")


if __name__ == "__main__":
    build()

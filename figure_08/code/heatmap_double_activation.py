#!/usr/bin/env python3
"""
Standalone heatmap: A119L+mutant activation (% of WT DAMGO), no cell numbers.

Rows: WT reference (top, divider) + the six A119L+partner doubles.
Cols: DAMGO / PZM21 / Nalbuphine / Naloxone. Value = activation window
(-Span_used, F-test-gated) as % of the WT DAMGO window. Purple = the double-
mutant bar colour. Final file is exactly 40 x 50 mm.

Source: figures/doubles_merged_points.csv
Outputs: plots/doubles/heatmap_double_activation.pdf / .png
"""
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from drc_fit import fit_drc

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
DOUBLE_CMAP = LinearSegmentedColormap.from_list("purple", ["#ffffff", "#9159a8"])


def build():
    df = pd.read_csv(ROOT / "figures" / "doubles_merged_points.csv")

    def W(lig, var):
        s = df[(df.ligand == lig) & (df.variant == f"{lig}_{var}")]
        return -fit_drc(s.logM.to_numpy(), s.response.to_numpy())["span_used"] \
            if not s.empty else np.nan

    ref = W("DAMGO", "WT")
    wt = [100 * W(l, "WT") / ref for l in LIGANDS]
    dbl = [[100 * W(l, f"119L_{t}") / ref for l in LIGANDS] for t, _ in PARTNERS]
    M = np.array([wt] + dbl)
    ylabels = ["WT"] + [lab for _, lab in PARTNERS]
    nr = M.shape[0]

    fig = plt.figure(figsize=(40 * MM, 50 * MM))
    axh = fig.add_axes([0.23, 0.30, 0.50, 0.58])
    im = axh.imshow(M, cmap=DOUBLE_CMAP, vmin=0, vmax=110, aspect="auto")
    axh.set_xticks(range(len(LIGANDS)))
    axh.set_xticklabels(LIGANDS, rotation=45, ha="right")
    axh.set_yticks(range(nr))
    axh.set_yticklabels(ylabels)
    axh.set_title("A119L + mutant\n(% DAMGO Emax)", pad=3)
    axh.set_xticks(np.arange(-.5, len(LIGANDS), 1), minor=True)
    axh.set_yticks(np.arange(-.5, nr, 1), minor=True)
    axh.grid(which="minor", color="white", linewidth=0.5)
    axh.tick_params(which="minor", length=0)
    axh.axhline(0.5, color="black", lw=0.5, zorder=6)   # WT / mutant divider
    for sp in axh.spines.values():
        sp.set_visible(True)
        sp.set_edgecolor("black")
        sp.set_linewidth(0.5)
        sp.set_zorder(6)

    pos = axh.get_position()
    cax = fig.add_axes([pos.x1 + 0.025, pos.y0, 0.045, pos.height])
    cb = fig.colorbar(im, cax=cax, ticks=[0, 50, 100])
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(width=0.5)

    fig.savefig(OUT / "heatmap_double_activation.pdf", transparent=True)
    fig.savefig(OUT / "heatmap_double_activation.png", dpi=300, transparent=True)
    plt.close(fig)
    print("wrote heatmap_double_activation.pdf / .png")


if __name__ == "__main__":
    build()

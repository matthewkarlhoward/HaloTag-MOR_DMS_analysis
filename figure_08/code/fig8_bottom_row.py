#!/usr/bin/env python3
"""
Figure 8 bottom row: K100N activation bars (left) + A119L-rescue heatmap (right).

Efficacy metric (both panels): the activation window of a variant = -Span_used
(0 if the fit is not a significant dose-response), expressed as a PERCENT of the
WT DAMGO window (WT DAMGO = 100% activation). This is the fraction of the full-
agonist response each variant reaches.

Left  : % activation of the K100N series (WT / A119L / K100N / A119L+K100N) per
        ligand, mean +/- SEM.
Right : the DOUBLE-mutant activation -- window of A119L+X in the same % units --
        for every partner X and ligand. This is exactly the purple bar on the
        left, generalised to all partners; purple heatmap matches the purple bar.
        Colour mnemonic: A119L (blue) + partner (red) = double (purple).

Source: figures/doubles_merged_points.csv
Outputs: plots/doubles/fig8_bottom_row.pdf / .png
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

LIGANDS = ["DAMGO", "PZM21", "Nalbuphine", "Naloxone"]     # bar + heatmap order

# ---- colours: reuse the figure's black / blue / red / purple ----
# WT = black, A119L = blue, partner single = red, double = purple (blue + red)
C_WT, C_A119L, C_SINGLE, C_DOUBLE = "#1a1a1a", "#2e6db4", "#c0392b", "#9159a8"
K100N = [("WT", "WT", C_WT), ("119L", "A119L", C_A119L),
         ("100N", "K100N", C_SINGLE), ("119L_100N", "A119L+K100N", C_DOUBLE)]
# white -> the exact double-mutant bar purple (light enough that black in-cell
# numbers stay legible on the darkest cells)
DOUBLE_CMAP = LinearSegmentedColormap.from_list("purple", ["#ffffff", C_DOUBLE])

PARTNERS = [("100D", "K100D"), ("100N", "K100N"), ("175N", "V175N"),
            ("278D", "R278D"), ("280A", "I280A"), ("280K", "I280K")]


def fit_of(df, lig, tok):
    s = df[(df.ligand == lig) & (df.variant == f"{lig}_{tok}")]
    return fit_drc(s.logM.to_numpy(), s.response.to_numpy()) if not s.empty else None


def window(df, lig, tok):
    f = fit_of(df, lig, tok)
    return -f["span_used"] if f else np.nan


def per_rep_span_sem(df, lig, tok):
    s = df[(df.ligand == lig) & (df.variant == f"{lig}_{tok}")]
    per = []
    for _, rep in s.groupby(["experiment", "rep"]):
        r = fit_drc(rep.logM.to_numpy(), rep.response.to_numpy())
        if r["ok"] and r["responsive"]:
            per.append(-r["span"])
    return np.std(per, ddof=1) / np.sqrt(len(per)) if len(per) > 1 else np.nan


def build():
    df = pd.read_csv(ROOT / "figures" / "doubles_merged_points.csv")
    ref = window(df, "DAMGO", "WT")            # 100% activation reference

    fig = plt.figure(figsize=(165 * MM, 40 * MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.62, 1.28],
                          wspace=0.30, left=0.078, right=0.93,
                          top=0.82, bottom=0.19)

    # ---------- left: K100N % activation bars ----------
    axb = fig.add_subplot(gs[0, 0])
    n = len(K100N)
    bw = 0.19
    offs = (np.arange(n) - (n - 1) / 2) * bw
    xg = np.arange(len(LIGANDS))
    for k, (tok, lab, col) in enumerate(K100N):
        hs, es = [], []
        for lig in LIGANDS:
            f = fit_of(df, lig, tok)
            resp = bool(f and f["responsive"])
            w = (-f["span_used"]) if f else np.nan
            hs.append(100 * w / ref)
            es.append(100 * per_rep_span_sem(df, lig, tok) / ref if resp else np.nan)
        axb.bar(xg + offs[k], hs, bw, yerr=es, color=col, edgecolor="black",
                linewidth=0.5, error_kw=dict(elinewidth=0.5, capsize=1.2,
                capthick=0.5), zorder=2, label=lab)
    axb.set_xticks(xg)
    axb.set_xticklabels(LIGANDS)
    axb.set_ylabel("Activation\n(% of WT DAMGO)", labelpad=2)
    axb.set_ylim(0, 118)
    axb.set_yticks([0, 25, 50, 75, 100])
    for s in ("top", "right"):
        axb.spines[s].set_visible(False)
    axb.legend(frameon=False, fontsize=6, ncol=4, loc="lower center",
               bbox_to_anchor=(0.5, 1.0), handlelength=0.9, handleheight=1.0,
               columnspacing=0.7, handletextpad=0.3)

    # ---------- right: activation heatmap; WT reference row on top ----------
    axh = fig.add_subplot(gs[0, 1])
    wt_row = [100 * window(df, lig, "WT") / ref for lig in LIGANDS]
    dbl = [[100 * window(df, lig, f"119L_{tok}") / ref for lig in LIGANDS]
           for tok, _ in PARTNERS]
    M = np.array([wt_row] + dbl)                      # row 0 = WT reference
    ylabels = ["WT"] + [lab for _, lab in PARTNERS]
    nr = M.shape[0]
    im = axh.imshow(M, cmap=DOUBLE_CMAP, vmin=0, vmax=110, aspect="auto")
    axh.set_xticks(range(len(LIGANDS)))
    axh.set_xticklabels(LIGANDS)                     # horizontal, bar order
    axh.set_yticks(range(nr))
    axh.set_yticklabels(ylabels)
    axh.set_title("A119L + mutant activation (% of WT DAMGO)", pad=4)
    for i in range(nr):
        for j in range(M.shape[1]):
            v = M[i, j]
            v = 0.0 if abs(v) < 0.5 else v            # avoid "-0"
            axh.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6,
                     color="black")
    axh.set_xticks(np.arange(-.5, len(LIGANDS), 1), minor=True)
    axh.set_yticks(np.arange(-.5, nr, 1), minor=True)
    axh.grid(which="minor", color="white", linewidth=0.5)
    axh.tick_params(which="minor", length=0)
    axh.plot([-.5, len(LIGANDS) - .5], [0.5, 0.5], color="black",
             lw=0.5, zorder=6)                        # WT / mutant divider
    for sp in axh.spines.values():                   # 0.5pt black border
        sp.set_visible(True)
        sp.set_edgecolor("black")
        sp.set_linewidth(0.5)
        sp.set_zorder(5)

    # colourbar placed tight against the heatmap
    pos = axh.get_position()
    cax = fig.add_axes([pos.x1 + 0.008, pos.y0, 0.013, pos.height])
    cb = fig.colorbar(im, cax=cax, ticks=[0, 50, 100])
    cb.set_label("%", labelpad=1)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(width=0.5)

    fig.savefig(OUT / "fig8_bottom_row.pdf")
    fig.savefig(OUT / "fig8_bottom_row.png", dpi=300)
    plt.close(fig)
    print("wrote fig8_bottom_row.pdf / .png")


if __name__ == "__main__":
    build()

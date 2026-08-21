#!/usr/bin/env python3
"""
Supplemental figures for the A119L double-mutant experiments.

  supp_doubles_series1  K100N / K100D / V175N
  supp_doubles_series2  R278D / I280A / I280K
      Span bars over the dose-response curves (WT / A119L / single / A119L+single)
      per ligand, grouped by position family so the severity pairs sit together.

  supp_doubles_rescue_all
      The rescue matrix of the main panel, extended to all six substitutions.
      I280A is included for completeness but is not a rescue test: it does not
      appreciably disrupt partial agonist signalling on its own.

Source: figures/doubles_merged_points.csv
Outputs: plots/doubles/supp_doubles_series1.pdf, supp_doubles_series2.pdf,
         supp_doubles_rescue_all.pdf  (+ .png)
"""
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from bar_100n_175_span_with_curves import COL, LIGANDS, MM, OUT, ROOT, draw_block
from fig_doubles_summary import draw_heatmap, rescue_matrix

PAGES = [
    ("supp_doubles_series1",
     [("K100N series", "100N", "K100N"),
      ("K100D series", "100D", "K100D"),
      ("V175N series", "175N", "V175N")]),
    ("supp_doubles_series2",
     [("R278D series", "278D", "R278D"),
      ("I280A series", "280A", "I280A"),
      ("I280K series", "280K", "I280K")]),
]

ALL_SUBS = [("100N", "K100N"), ("100D", "K100D"), ("175N", "V175N"),
            ("278D", "R278D"), ("280A", "I280A"), ("280K", "I280K")]


def series_legend(fig, anchor=(0.805, 0.5)):
    handles = [
        Patch(facecolor=COL["WT"], edgecolor="black", lw=0.5, label="WT"),
        Patch(facecolor=COL["A119L"], edgecolor="black", lw=0.5, label="A119L"),
        Patch(facecolor=COL["single"], edgecolor="black", lw=0.5, label="single mutant"),
        Patch(facecolor=COL["double"], edgecolor="black", lw=0.5, label="A119L + single"),
        Line2D([], [], color="black", lw=0.9, ls=(0, (2, 1.5)),
               label="no response (F-test n.s.)"),
        Line2D([], [], color="black", marker="o", ms=1.5, lw=0, label="mean ± SEM"),
    ]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=anchor,
               frameon=False, fontsize=5.6, handlelength=1.3, handleheight=1.2)


def build_series(df, name, blocks):
    nblk = len(blocks)
    H = 22 + 84 * nblk
    fig = plt.figure(figsize=(150 * MM, H * MM))
    gs = fig.add_gridspec(2 * nblk, len(LIGANDS),
                          height_ratios=[1.0, 1.15] * nblk,
                          hspace=0.55, wspace=0.22,
                          left=0.09, right=0.80,
                          top=1 - 8 / H, bottom=12 / H)
    for b, (title, tok, label) in enumerate(blocks):
        draw_block(fig, gs, 2 * b, 2 * b + 1, df, title, tok, label)
    series_legend(fig)
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=300)
    plt.close(fig)
    print(f"wrote {name}.pdf / .png")


def build_rescue_all(df):
    frac, single_ok, dbl_ok = rescue_matrix(df, ALL_SUBS)
    fig = plt.figure()
    fig.set_size_inches(100 * MM, 56 * MM)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.075], wspace=0.06,
                          left=0.155, right=0.855, top=0.85, bottom=0.20)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])
    draw_heatmap(ax, cax, frac, single_ok, dbl_ok, ALL_SUBS)
    ax.set_title("Signaling in the A119L background, relative to A119L alone",
                 loc="left", fontsize=7, pad=3)
    fig.legend(handles=[
        Line2D([], [], marker="o", ms=1.8, mfc="none", mec="black", mew=0.5, lw=0,
               label="single mutant alone still signals"),
    ], loc="lower center", bbox_to_anchor=(0.5, 0.005), frameon=False, fontsize=5.6,
        handlelength=1.3)
    fig.savefig(OUT / "supp_doubles_rescue_all.pdf")
    fig.savefig(OUT / "supp_doubles_rescue_all.png", dpi=300)
    plt.close(fig)
    print("wrote supp_doubles_rescue_all.pdf / .png")


def build():
    df = pd.read_csv(ROOT / "figures" / "doubles_merged_points.csv")
    for name, blocks in PAGES:
        build_series(df, name, blocks)
    build_rescue_all(df)


if __name__ == "__main__":
    build()

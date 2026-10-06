#!/usr/bin/env python3
"""
Supplementary normalised DRC grids for the K100 / V175 / R278 TruPath plates:
every ligand each plate actually ran, not just the three shown in the main
figure. Same normalisation and F-test fitting as drc_singles_grid.py.

  singles_supp_icl.pdf    20260125 ICL plate, read 2 -- 6 ligands
                          rows: K100D/K100N and V175E/V175N against that
                          plate's own WT
  singles_supp_r278.pdf   20260121 plate -- 6 ligands, WT vs R278D. An earlier,
                          independent R278D run (no PZM21); the main figure's
                          R278D comes from the 20260427 plate.

Outputs: plots/singles_trupath/singles_supp_{icl,r278}.{pdf,png}
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from drc_singles_grid import (MM, REPO, XMIN, XMAX, XPAD_LO, XPAD_HI,
                              C_WT, C_M1, C_M2, draw)

OUT = Path(__file__).resolve().parent

FIGURES = [
    ("singles_supp_icl", "ICL",
     ["DAMGO", "Morphine", "Fentanyl", "Butorphanol", "PZM21", "Nalbuphine"],
     [("K100", [("WT", C_WT, "WT"), ("K100D", C_M1, "K100D"),
                ("K100N", C_M2, "K100N")]),
      ("V175", [("WT", C_WT, "WT"), ("V175E", C_M1, "V175E"),
                ("V175N", C_M2, "V175N")])]),
    ("singles_supp_r278", "R278_run1",
     ["DAMGO", "Morphine", "Fentanyl", "Buprenorphine", "MP", "Nalbuphine"],
     [("R278", [("WT", C_WT, "WT"), ("R278D", C_M1, "R278D")])]),
]


def build_one(stem, source, ligands, rows, df):
    ncol, nrow = len(ligands), len(rows)
    pw, ph, gx, gy = 25 * MM, 21 * MM, 5 * MM, 6.5 * MM
    left, bottom, top_pad = 13 * MM, 11 * MM, 6 * MM
    fig_w = left + ncol * pw + (ncol - 1) * gx + 4 * MM
    fig_h = bottom + nrow * ph + (nrow - 1) * gy + top_pad

    fig = plt.figure()
    fig.set_size_inches(fig_w, fig_h)
    params = []
    for i, (row_label, variants) in enumerate(rows):
        for j, lig in enumerate(ligands):
            x0 = left + j * (pw + gx)
            y0 = fig_h - top_pad - (i + 1) * ph - i * gy
            ax = fig.add_axes([x0 / fig_w, y0 / fig_h, pw / fig_w, ph / fig_h])
            for variant, color, label in variants:
                sub = df[(df.source == source) & (df.ligand == lig)
                         & (df.variant == variant)]
                fit = draw(ax, sub, color)
                if fit is not None:
                    params.append(dict(
                        source=source, row=row_label, ligand=lig,
                        variant=label, responsive=fit["responsive"],
                        logec50=fit["params"][2] if fit["params"] else np.nan,
                        emax=fit["params"][1] if fit["params"] else np.nan,
                        window=-fit["span_used"], r2=fit["r2"],
                        pval=fit["pval"]))
            ax.axhline(1.0, color="black", lw=0.4, ls=(0, (1.5, 1.5)), zorder=1)
            ax.set_xlim(XPAD_LO, XPAD_HI)
            ax.set_ylim(0.45, 1.2)
            ax.set_yticks([0.5, 0.75, 1.0])
            ax.set_xticks(np.arange(-12, -3, 2))
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            if i == 0:
                ax.set_title(lig, pad=2.5)
            if j == 0:
                ax.set_ylabel("Norm. BRET", labelpad=1.5)
            else:
                ax.set_yticklabels([])
            if i == nrow - 1:
                ax.set_xlabel("log[ligand], M", labelpad=1.5)
            if j == 0:
                handles = [Line2D([], [], color=c, marker="o", ms=1.9, lw=0.8,
                                  label=lab) for _, c, lab in variants]
                ax.legend(handles=handles, loc="lower left", frameon=False,
                          handlelength=1.1, handletextpad=0.4, borderpad=0.0,
                          labelspacing=0.15, fontsize=5.5)

    fig.savefig(OUT / f"{stem}.pdf", transparent=True)
    fig.savefig(OUT / f"{stem}.png", dpi=400, facecolor="white")
    plt.close(fig)
    print(f"-> {(OUT / stem).relative_to(REPO)}.pdf")
    return params


def main():
    df = pd.read_csv(REPO / "figures" / "singles_norm_points.csv")
    df = df[(df.logM >= XMIN) & (df.logM <= XMAX)]
    allp = []
    for stem, source, ligands, rows in FIGURES:
        allp += build_one(stem, source, ligands, rows, df)
    tab = pd.DataFrame(allp)
    tab.to_csv(OUT / "singles_supp_params.csv", index=False)
    print(tab.round(3).to_string(index=False))


if __name__ == "__main__":
    main()

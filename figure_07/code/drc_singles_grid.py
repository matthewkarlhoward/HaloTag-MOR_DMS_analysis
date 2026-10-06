#!/usr/bin/env python3
"""
K100 / V175 / R278 dose-response panels on the NORMALISED scale.

Same data as panels d-f of the G-protein figure (which plot raw 515/410 BRET),
re-drawn after the per-curve normalisation used for the WT / A119L / double
mutants: every curve divided by its own fitted no-drug plateau, so baseline =
1.0 and agonist response reads downward (signal-down assay).

  row 1  K100D / K100N   20260427 plate  (the plate panel d plots; the same
                         plate behind figures/doubles1_points.csv)
  row 2  V175E / V175N   20260125 ICL plate, read 2
  row 3  R278D           20260427 plate

Curves are re-fit on the normalised points with the doubles F-test rule
(plots/doubles/drc_fit.py): a sigmoid is drawn only where it beats a flat line
at p < 0.05, otherwise the curve is drawn flat at the data mean.

Inputs : figures/singles_norm_points.csv   (normalize_singles.py)
Outputs: plots/singles_trupath/singles_norm_grid.{pdf,png}
"""
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "plots" / "doubles"))
from drc_fit import fit_drc, sigmoid            # noqa: E402

mpl.rcParams.update({
    "font.family": "Helvetica", "font.weight": "normal", "font.size": 6,
    "axes.titlesize": 6, "axes.labelsize": 6, "xtick.labelsize": 6,
    "ytick.labelsize": 6, "axes.titleweight": "normal",
    "axes.labelweight": "normal", "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "text.color": "black",
    "axes.edgecolor": "black", "axes.labelcolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "mathtext.default": "regular", "pdf.fonttype": 42, "ps.fonttype": 42,
})

MM = 1 / 25.4
OUT = Path(__file__).resolve().parent

LIGANDS = ["DAMGO", "PZM21", "Nalbuphine"]
C_WT, C_M1, C_M2 = "#000000", "#a4262c", "#3b7fc4"

# (row label, source, [(variant, colour, legend label), ...])
ROWS = [
    ("K100", "doubles_raw", [("WT", C_WT, "WT"),
                             ("100D", C_M1, "K100D"),
                             ("100N", C_M2, "K100N")]),
    ("V175", "ICL",         [("WT", C_WT, "WT"),
                             ("V175E", C_M1, "V175E"),
                             ("V175N", C_M2, "V175N")]),
    ("R278", "doubles_raw", [("WT", C_WT, "WT"),
                             ("278D", C_M1, "R278D")]),
]

XMIN, XMAX = -12.0, -4.0          # data / fitted-curve range
XPAD_LO, XPAD_HI = -12.5, -3.5    # axis range: half a log of air each side so
                                  # the first and last points are not clipped
XGRID = np.linspace(XMIN, XMAX, 400)


def draw(ax, sub, color):
    """Mean +/- SEM points and the F-tested fit for one curve."""
    if sub is None or sub.empty:
        return None
    prof = sub.groupby("logM")["norm"].agg(["mean", "sem"]).reset_index()
    fit = fit_drc(sub.logM.to_numpy(), sub.norm.to_numpy())
    if fit["ok"] and fit["responsive"]:
        ax.plot(XGRID, sigmoid(XGRID, *fit["params"]), color=color, lw=0.8,
                zorder=3, solid_capstyle="round")
    else:
        ax.plot([XMIN, XMAX], [fit["level"]] * 2, color=color, lw=0.8,
                zorder=3, solid_capstyle="round")
    ax.errorbar(prof["logM"], prof["mean"], yerr=prof["sem"], fmt="o", ms=1.9,
                mfc=color, mec=color, ecolor=color, elinewidth=0.5, capsize=0,
                lw=0, zorder=4)
    return fit


def build():
    df = pd.read_csv(REPO / "figures" / "singles_norm_points.csv")
    df = df[(df.logM >= XMIN) & (df.logM <= XMAX)]

    # exact 105 x 110 mm canvas; panel size falls out of the margins so the
    # grid always fills it (see reference_matplotlib_exact_mm_figsize)
    fig_w, fig_h = 105 * MM, 110 * MM
    left, right, bottom, top_pad = 12 * MM, 3 * MM, 10 * MM, 5 * MM
    gx, gy = 5 * MM, 7 * MM
    ncol, nrow = len(LIGANDS), len(ROWS)
    pw = (fig_w - left - right - (ncol - 1) * gx) / ncol
    ph = (fig_h - bottom - top_pad - (nrow - 1) * gy) / nrow

    fig = plt.figure()
    fig.set_size_inches(fig_w, fig_h)
    rows_out = []

    for i, (row_label, source, variants) in enumerate(ROWS):
        for j, lig in enumerate(LIGANDS):
            x0 = left + j * (pw + gx)
            y0 = fig_h - top_pad - (i + 1) * ph - i * gy
            ax = fig.add_axes([x0 / fig_w, y0 / fig_h, pw / fig_w, ph / fig_h])
            for variant, color, label in variants:
                sub = df[(df.source == source) & (df.ligand == lig)
                         & (df.variant == variant)]
                fit = draw(ax, sub, color)
                if fit is not None:
                    rows_out.append(dict(
                        row=row_label, source=source, ligand=lig,
                        variant=label, responsive=fit["responsive"],
                        logec50=fit["params"][2] if fit["params"] else np.nan,
                        emax=fit["params"][1] if fit["params"] else np.nan,
                        window=-fit["span_used"], r2=fit["r2"],
                        pval=fit["pval"]))
            ax.axhline(1.0, color="black", lw=0.4, ls=(0, (1.5, 1.5)),
                       zorder=1)
            ax.set_xlim(XPAD_LO, XPAD_HI)
            ax.set_ylim(0.45, 1.15)
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
            if i == len(ROWS) - 1:
                ax.set_xlabel(f"log[{lig}], M", labelpad=1.5)
            handles = [Line2D([], [], color=c, marker="o", ms=1.9, lw=0.8,
                              label=lab) for _, c, lab in variants]
            if j == 0:
                ax.legend(handles=handles, loc="lower left", frameon=False,
                          handlelength=1.1, handletextpad=0.4,
                          borderpad=0.0, labelspacing=0.15, fontsize=5.5)

    fig.savefig(OUT / "singles_norm_grid.pdf", transparent=True)
    fig.savefig(OUT / "singles_norm_grid.png", dpi=400, facecolor="white")
    plt.close(fig)

    tab = pd.DataFrame(rows_out)
    tab.to_csv(OUT / "singles_norm_grid_params.csv", index=False)
    print(tab.round(3).to_string(index=False))
    print(f"-> {(OUT / 'singles_norm_grid.pdf').relative_to(REPO)}")


if __name__ == "__main__":
    build()

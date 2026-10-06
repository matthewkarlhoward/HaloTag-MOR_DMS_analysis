#!/usr/bin/env python3
"""Plate maps for the 20260502 TruPath run, so the row structure is visible.

Prints a text grid and writes a two-panel figure:
  donor 410 nm  -- the luminescence available in each well
  BRET ratio    -- 515/410, the readout (LOWER = more receptor activation)

Layout recovered earlier and validated against the .prism value-for-value:
  rows    DAMGO A/F/K, buprenorphine B/C/G/H/L/M, SR-17018 D/E/I/J/N/O, P = blanks
  columns dose pairs, 1-2 = 1e-4 M down to 21-22 = 1e-14 M, 23-24 = vehicle

The point of looking: ligand assignment is confounded with a strong row gradient in
donor signal (row O ~6,900, row G ~24,300), and SR-17018 drew the weak rows.
"""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
WELLS = HERE / "sr17018_trupath_wells.csv"
READ = "READ2"
ROWS = list("ABCDEFGHIJKLMNOP")
DOSE_LABEL = ["-4", "-5", "-6", "-7", "-8", "-9", "-10", "-11",
              "-12", "-13", "-14", "veh"]
LIG_COLOR = {"DAMGO": "#000000", "Buprenorphine": "#2E75B6",
             "SR17018": "#9E2A2B", "blank": "#999999"}

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 6, "axes.labelsize": 6, "xtick.labelsize": 5, "ytick.labelsize": 5,
    "text.color": "black", "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 1.5, "ytick.major.size": 1.5,
})


def load():
    donor = np.full((16, 24), np.nan)
    bret = np.full((16, 24), np.nan)
    lig = {}
    with open(WELLS) as fh:
        for r in csv.DictReader(fh):
            if r["read"] != READ:
                continue
            i, j = ROWS.index(r["plate_row"]), int(r["plate_col"]) - 1
            donor[i, j] = float(r["donor_410"])
            if r["bret_ratio"]:
                bret[i, j] = float(r["bret_ratio"])
            lig[r["plate_row"]] = r["ligand"] or "blank"
    return donor, bret, lig


def text_grid(donor, lig):
    print("donor 410 nm, in thousands (READ2). '.' = blank row\n")
    head = "      " + "".join(f"{c:>3d}" for c in range(1, 25))
    print(head)
    print("      " + "".join(f"{DOSE_LABEL[(c - 1) // 2]:>3s}" if c % 2 else "   "
                             for c in range(1, 25)))
    for i, R in enumerate(ROWS):
        cells = ""
        for j in range(24):
            v = donor[i, j]
            cells += "  ." if np.isnan(v) or v < 300 else f"{v / 1000:3.0f}"
        print(f"  {R} {cells}   {lig.get(R, '')}")
    print("\n  (dose label sits over the first column of each pair)")


def panel(ax, M, title, cmap, lig, invert=False):
    data = M.copy()
    data[15, :] = np.nan                       # blank row, do not scale to it
    finite = data[np.isfinite(data)]
    vmin, vmax = np.percentile(finite, 2), np.percentile(finite, 98)
    im = ax.imshow(data, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation="nearest")
    ax.set_title(title, fontsize=6, pad=3)
    ax.set_xticks(np.arange(0, 24, 2))
    ax.set_xticklabels(DOSE_LABEL, rotation=90)
    ax.set_yticks(range(16))
    ax.set_yticklabels(ROWS)
    for i, R in enumerate(ROWS):
        ax.get_yticklabels()[i].set_color(LIG_COLOR.get(lig.get(R, "blank"), "#999999"))
    # separate the dose pairs
    for c in range(2, 24, 2):
        ax.axvline(c - 0.5, color="white", linewidth=0.4)
    ax.tick_params(length=1.5)
    cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=4.5, length=1.5, width=0.4)
    return im


def main():
    donor, bret, lig = load()
    text_grid(donor, lig)

    fig = plt.figure()
    fig.set_size_inches(170 / 25.4, 62 / 25.4)
    ax1 = fig.add_axes([0.055, 0.17, 0.40, 0.70])
    ax2 = fig.add_axes([0.555, 0.17, 0.40, 0.70])
    panel(ax1, donor, "donor 410 nm (luminescence per well)", "Blues", lig)
    panel(ax2, bret, "BRET ratio 515/410 (lower = more activated)", "Purples_r", lig)
    for ax in (ax1, ax2):
        ax.set_xlabel("log[ligand] (M)", fontsize=6, labelpad=1)
    fig.text(0.055, 0.015,
             "row labels coloured by ligand: DAMGO black, buprenorphine blue, "
             "SR-17018 red, P = cell-free blanks", fontsize=5)
    fig.savefig(HERE / "plate_map_20260502.pdf", transparent=True)
    fig.savefig(HERE / "plate_map_20260502.png", dpi=600, facecolor="white")
    plt.close(fig)
    print("\nwrote plate_map_20260502.pdf / .png")


if __name__ == "__main__":
    main()

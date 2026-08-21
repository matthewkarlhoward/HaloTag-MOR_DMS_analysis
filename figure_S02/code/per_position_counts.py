#!/usr/bin/env python3
"""
Per-position total counts for the baseline sample (baseline_1 + baseline_2).

Output: lib_generation/per_position_counts.pdf  (60 mm wide x 30 mm tall)
"""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

MM = 1.0 / 25.4  # mm -> inches

mpl.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 6,
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.edgecolor": "black",
    "axes.titlesize": 6,
    "axes.labelsize": 6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.linewidth": 0.4,
    "xtick.major.width": 0.4,
    "ytick.major.width": 0.4,
    "xtick.major.size": 2,
    "ytick.major.size": 2,
})

HERE = Path(__file__).resolve().parent
B1 = pd.read_csv(HERE / "baseline_1.csv")
B2 = pd.read_csv(HERE / "baseline_2.csv")

baseline = B1[["pos"]].copy()
baseline["count"] = B1["count"].values + B2["count"].values
bl_by_pos = baseline.groupby("pos")["count"].sum()

fig = plt.figure(figsize=(60 * MM, 30 * MM))
ax = fig.add_axes([0.16, 0.27, 0.80, 0.66])  # left, bottom, width, height (figure fraction)

ax.bar(bl_by_pos.index, bl_by_pos.values, color="#888888", width=1.0, linewidth=0)

pos_min, pos_max = int(bl_by_pos.index.min()), int(bl_by_pos.index.max())
grid_start = ((pos_min + 49) // 50) * 50
for x in range(grid_start, pos_max + 1, 50):
    ax.axvline(x, color="black", lw=0.3, ls="--", alpha=0.5)

ax.set_yscale("log")
ax.set_ylim(1, 2000)
ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())
ax.set_yticks([1, 10, 100, 1000])
ax.set_xlim(0, pos_max + 0.5)
ax.set_xlabel("Position")
ax.set_ylabel("Total Count")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

out = HERE / "per_position_counts.pdf"
fig.savefig(out, dpi=600)  # no bbox_inches="tight" — preserve exact 60x30 mm
print(f"median: {bl_by_pos.median():.0f}   min: {bl_by_pos.min()}   max: {bl_by_pos.max()}")
print(f"Wrote {out}")

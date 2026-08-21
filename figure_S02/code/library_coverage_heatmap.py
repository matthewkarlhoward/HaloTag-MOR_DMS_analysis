#!/usr/bin/env python3
"""
Library coverage heatmap from baseline counts.

- Sums counts across baseline_1.csv and baseline_2.csv (per-variant).
- Plots pos (x) x variant (y), colored by summed count.
- Reports % of designed variants observed (count>0), overall and missense-only.

Output: lib_generation/library_coverage_heatmap.pdf
"""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm

mpl.rcParams.update({
    "font.family": "Helvetica", "font.size": 6, "font.weight": "normal",
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
})

HERE = Path(__file__).resolve().parent
B1 = pd.read_csv(HERE / "baseline_1.csv")
B2 = pd.read_csv(HERE / "baseline_2.csv")
MV = pd.read_csv(HERE / "mor_variants.csv")

# Merge on the variant identity columns and sum counts
key_cols = ["pos", "mutation_type", "name", "mutation", "length", "hgvs"]
b = B1[key_cols].copy()
b["count"] = B1["count"].values + B2["count"].values

# Build a "row" label for the y-axis
# - M / S: amino acid letter
# - D: del1 / del2 / del3
# - I: insG / insGG / insGGG (Gs added — single-letter glycine per length)
def variant_row(r):
    mt, mut, ln = r["mutation_type"], r["mutation"], int(r["length"])
    if mt in ("M", "S"):
        return mut
    if mt == "D":
        return f"del{ln}"
    if mt == "I":
        return "ins" + ("G" * ln)
    return None

b["row"] = b.apply(variant_row, axis=1)

# Row order: insertions, deletions, then AAs grouped by property
indel_order = ["insGGG", "insGG", "insG", "del3", "del2", "del1"]
aa_order = list("PRKHDEFWYNQCSTILVMGA")  # property-grouped (matches plots/heatmaps style)
row_order = indel_order + aa_order
n_rows = len(row_order)

pos_min, pos_max = int(b["pos"].min()), int(b["pos"].max())
positions = np.arange(pos_min, pos_max + 1)
n_cols = len(positions)

# Pivot to row x position matrix (summed counts)
mat = (
    b.pivot_table(index="row", columns="pos", values="count", aggfunc="sum")
    .reindex(index=row_order, columns=positions)
    .to_numpy(dtype=float)
)

# Mark cells that were NOT in the design (NaN) vs. designed but unobserved (0)
designed_mask = (
    MV.assign(row=MV.apply(variant_row, axis=1))
    .pivot_table(index="row", columns="pos", values="count", aggfunc="size")
    .reindex(index=row_order, columns=positions)
    .notna()
    .to_numpy()
)
plot_mat = np.where(designed_mask, mat, np.nan)
# Designed-but-zero -> very small positive so LogNorm still shows them as low end
plot_mat_for_log = np.where(designed_mask & (plot_mat == 0), 0.5, plot_mat)

# ---- Coverage stats ----
n_designed = len(MV)
observed = b["count"] > 0
n_observed = int(observed.sum())
pct_overall = 100 * n_observed / n_designed

mis_designed = (MV["mutation_type"] == "M").sum()
mis_observed = int(((b["mutation_type"] == "M") & observed).sum())
pct_missense = 100 * mis_observed / mis_designed

print(f"Designed variants total:    {n_designed}")
print(f"Observed variants (count>0): {n_observed}  ({pct_overall:.2f}%)")
print(f"Designed missense:          {mis_designed}")
print(f"Observed missense:          {mis_observed}  ({pct_missense:.2f}%)")

# ---- Plot ----
cmap = LinearSegmentedColormap.from_list(
    "coverage", ["#f5f5f5", "#cfe1ee", "#3b8db5", "#0a3d62"], N=256
)
cmap.set_bad("#dddddd")  # not-in-design cells

vmax = np.nanpercentile(plot_mat_for_log, 99)
vmax = max(vmax, 2.0)

fig, ax = plt.subplots(figsize=(11, 3.6))
im = ax.imshow(
    plot_mat_for_log,
    aspect="auto",
    cmap=cmap,
    norm=LogNorm(vmin=0.5, vmax=vmax),
    interpolation="nearest",
    extent=[pos_min - 0.5, pos_max + 0.5, n_rows - 0.5, -0.5],
)

# Y ticks
ax.set_yticks(np.arange(n_rows))
ax.set_yticklabels(row_order, fontsize=5)

# X ticks every 20
xt = np.arange(((pos_min + 19) // 20) * 20, pos_max + 1, 20)
ax.set_xticks(xt)
ax.set_xticklabels(xt, fontsize=5)
ax.set_xlabel("Residue position", fontsize=7)
ax.set_ylabel("Variant", fontsize=7)

# Separator between indels and AAs
ax.axhline(len(indel_order) - 0.5, color="black", lw=0.6)

ax.set_title(
    f"MOR library coverage (baseline_1 + baseline_2)   "
    f"observed: {n_observed}/{n_designed} ({pct_overall:.1f}%)   "
    f"missense: {mis_observed}/{mis_designed} ({pct_missense:.1f}%)",
    fontsize=7, loc="left",
)

cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
cbar.set_label("Summed read count (log)", fontsize=6)
cbar.ax.tick_params(labelsize=5)

# Legend for grey = not in design
from matplotlib.patches import Patch
ax.legend(
    handles=[Patch(facecolor="#dddddd", edgecolor="black", lw=0.4, label="not in design")],
    loc="upper right", bbox_to_anchor=(1.0, 1.18), frameon=False, fontsize=5,
)

plt.tight_layout()
out = HERE / "library_coverage_heatmap.pdf"
plt.savefig(out, dpi=300, bbox_inches="tight")
print(f"\nWrote {out}")

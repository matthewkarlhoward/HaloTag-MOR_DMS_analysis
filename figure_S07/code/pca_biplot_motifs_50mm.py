#!/usr/bin/env python3
"""
50x50 mm PC1 x PC2 position scatter; canonical GPCR motif residues highlighted by
colour (no per-residue text labels; motifs identified via the legend).
All text size 6; 0.5 pt lines. Reads precomputed position scores; no PCA recompute.

Output: plots/pca/pca_biplot_motifs_50mm.pdf
"""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "plots" / "pca" / "position"
MM = 1/25.4

# House style: black text, Helvetica, size 6, 0.5 pt lines
plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 6,
    "mathtext.default": "regular",
    "text.color": "black",
    "axes.labelsize": 6, "axes.labelcolor": "black", "axes.edgecolor": "black",
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "xtick.color": "black", "ytick.color": "black",
    "legend.fontsize": 6,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 1.5, "ytick.major.size": 1.5,
})

pos = pd.read_csv(OUT / "pca_position_scores.csv")

MOTIFS = {
    "DRY (3x49-51)":    [166, 167, 168],
    "CWxP (6x47-50)":   [294, 295, 296, 297],
    "NPxxY (7x49-53)":  [334, 335, 336, 337, 338],
    "PIF toggle":       [157, 246, 291],
    "Na+ site (D2x50)": [116],
}
MOTIF_COL = {
    "DRY (3x49-51)":    "#9b1c1d",
    "CWxP (6x47-50)":   "#0e6d39",
    "NPxxY (7x49-53)":  "#3779b9",
    "PIF toggle":       "#e66101",
    "Na+ site (D2x50)": "#6f08a3",
}
pos_to_motif = {p: m for m, ps in MOTIFS.items() for p in ps}

fig, ax = plt.subplots(figsize=(45*MM, 45*MM), constrained_layout=True)
fig.set_size_inches(45*MM, 45*MM)  # exact: figsize= arg rounds to 2 dp

ax.scatter(pos.PC1, pos.PC2, s=2.5, color="#a8a8a8", alpha=0.6,
           edgecolors="none", zorder=1)

for p, m in pos_to_motif.items():
    sub = pos[pos.position == p]
    if len(sub) == 0:
        continue
    r = sub.iloc[0]
    ax.scatter(r.PC1, r.PC2, s=11, color=MOTIF_COL[m],
               edgecolors="black", linewidths=0.3, zorder=4)

ax.axhline(0, color="grey", lw=0.5, ls="--")
ax.axvline(0, color="grey", lw=0.5, ls="--")

ax.set_xlim(pos.PC1.min() - 1.0, pos.PC1.max() + 2.0)
ax.set_ylim(pos.PC2.min() - 0.6, pos.PC2.max() + 0.8)

ax.set_xlabel("PC1 (66.8%)")
ax.set_ylabel("PC2 (15.6%)")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Legend: motif name only (BW position numbers -> caption); superscript Na+
SHORT = {"DRY (3x49-51)": "DRY", "CWxP (6x47-50)": "CWxP",
         "NPxxY (7x49-53)": "NPxxY", "PIF toggle": "PIF",
         "Na+ site (D2x50)": r"Na$^+$"}
handles = [mpatches.Patch(color=c, label=SHORT[m]) for m, c in MOTIF_COL.items()]
ax.legend(handles=handles, frameon=False, loc="lower left",
          handlelength=1.0, handletextpad=0.4, labelspacing=0.3,
          borderpad=0.2, borderaxespad=0.2)

fig.savefig(OUT / "pca_biplot_motifs_50mm.pdf", dpi=600)
plt.close()
print("Saved: pca_biplot_motifs_50mm.pdf")

#!/usr/bin/env python3
"""
50x50 mm TRUE biplot: position scores (dots) + drug loadings (arrows) on PC1 x PC2.
 - grey dots = all positions; colored dots = canonical GPCR motif residues (labeled)
 - arrows = drug loadings, colored by efficacy class, scaled by one isotropic factor k
Reads precomputed scores + loadings; does not recompute PCA.

Output: plots/pca/pca_biplot_motifs_loadings_50mm.pdf
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from adjustText import adjust_text

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "plots" / "pca" / "position"
MM = 1/25.4

plt.rcParams.update({
    "font.family": "Helvetica",
    "axes.linewidth": 0.4,
    "xtick.major.width": 0.4, "ytick.major.width": 0.4,
    "xtick.major.size": 1.5, "ytick.major.size": 1.5,
})

# ── Load precomputed scores + loadings ───────────────────────────────────────
pos = pd.read_csv(OUT / "pca_position_scores.csv")
load = pd.read_csv(OUT / "pca_loadings_table.csv", index_col=0)
gpcrdb = pd.read_csv(ROOT / "annotations/GPCRdb_OPRM1_table.csv").rename(
    columns={"pos": "position"})
wt_lookup = dict(zip(gpcrdb.position, gpcrdb.wt_aa))

# ── Canonical GPCR motifs (Na+ adj & TM6 ionic-lock dropped per request) ─────
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

# ── Drug efficacy-class palette (matches pca_analysis.py) ────────────────────
CLS_COL = {"strong": "#231f20", "intermediate": "#9b1c1d",
           "weak": "#3779b9", "antagonist": "#6f08a3"}
CLS_ORDER = ["strong", "intermediate", "weak", "antagonist"]

# Isotropic loading scale: one factor for BOTH axes so arrow angles are preserved.
score_maxabs = np.abs(pos[["PC1", "PC2"]].to_numpy()).max(axis=0)
load_maxabs = np.abs(load[["PC1", "PC2"]].to_numpy()).max(axis=0)
k = 0.80 * float(np.min(score_maxabs / load_maxabs))
print(f"loading scale factor k = {k:.3f}")

# ── Figure ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(50*MM, 50*MM), constrained_layout=True)
fig.set_size_inches(50*MM, 50*MM)  # exact: figsize= arg rounds to 2 dp

# All positions as grey background
ax.scatter(pos.PC1, pos.PC2, s=2.5, color="#a8a8a8", alpha=0.6,
           edgecolors="none", zorder=1)

# Drug loading arrows (under the motif dots); labels collected for declutter
drug_texts = []
for d, row in load.iterrows():
    col = CLS_COL.get(row["class"], "#333333")
    x, y = row["PC1"]*k, row["PC2"]*k
    ax.annotate("", xy=(x, y), xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", lw=0.5, color=col,
                                mutation_scale=3, shrinkA=0, shrinkB=0),
                zorder=3)
    drug_texts.append(ax.text(x, y, d, fontsize=2.4, color=col,
                              fontweight="bold", ha="center", va="center",
                              zorder=6))

# Motif residues colored + labeled (on top)
for p, m in pos_to_motif.items():
    sub = pos[pos.position == p]
    if len(sub) == 0:
        continue
    r = sub.iloc[0]
    ax.scatter(r.PC1, r.PC2, s=11, color=MOTIF_COL[m],
               edgecolors="black", linewidths=0.3, zorder=4)
    ax.annotate(f"{wt_lookup.get(p, '?')}{p}", (r.PC1, r.PC2),
                fontsize=3, fontweight="bold", color=MOTIF_COL[m],
                xytext=(2, 1.5), textcoords="offset points",
                clip_on=False, zorder=5)

ax.axhline(0, color="grey", lw=0.3, ls="--")
ax.axvline(0, color="grey", lw=0.3, ls="--")

# Limits: include cloud + scaled arrow tips; pad right/top for edge labels
xs = list(pos.PC1) + list(load.PC1*k)
ys = list(pos.PC2) + list(load.PC2*k)
ax.set_xlim(min(xs) - 1.0, max(xs) + 2.0)
ax.set_ylim(min(ys) - 0.6, max(ys) + 0.8)

# Repel the drug labels apart (degenerate strong-agonist cluster) with leaders
try:
    adjust_text(drug_texts, ax=ax, expand=(1.3, 1.5),
                arrowprops=dict(arrowstyle="-", color="0.55", lw=0.2))
except Exception as e:
    print(f"adjust_text skipped: {e}")

ax.set_xlabel("PC1", fontsize=5)
ax.set_ylabel("PC2", fontsize=5)
ax.tick_params(labelsize=4)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Two legends, titled to disambiguate the shared palette
motif_handles = [mpatches.Patch(color=c, label=m) for m, c in MOTIF_COL.items()]
leg1 = ax.legend(handles=motif_handles, frameon=False, fontsize=3,
                 loc="lower right", title="Motif (position)", title_fontsize=3,
                 handlelength=1.0, handletextpad=0.4, labelspacing=0.2,
                 borderpad=0.2, borderaxespad=0.2)
ax.add_artist(leg1)
cls_handles = [Line2D([0], [0], color=CLS_COL[c], lw=1, label=c)
               for c in CLS_ORDER]
ax.legend(handles=cls_handles, frameon=False, fontsize=3,
          loc="upper left", title="Drug loading (efficacy)", title_fontsize=3,
          handlelength=1.2, handletextpad=0.4, labelspacing=0.2,
          borderpad=0.2, borderaxespad=0.2)

fig.savefig(OUT / "pca_biplot_motifs_loadings_50mm.pdf", dpi=600)
plt.close()
print("Saved: pca_biplot_motifs_loadings_50mm.pdf")

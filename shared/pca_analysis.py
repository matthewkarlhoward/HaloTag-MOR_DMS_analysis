#!/usr/bin/env python3
"""
PCA on per-position mean missense effects across 15 ligands (MCAM and FSK excluded).
Positions restricted to 65-355.

Each row = one position; each column = mean missense effect for one ligand.

Outputs (plots/pca/position/):
  scree_position.pdf
  pca_position_scatter.pdf   — PC1 vs PC2 scatter, colored by SSE
  pca_drug_loadings.pdf      — drug loading vectors in PC1/PC2 space
  pca_receptor_map.pdf       — PC1/PC2/PC3 along receptor sequence
  pca_loadings_table.csv
  pca_position_scores.csv
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "plots" / "pca" / "position"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 6,
    "axes.labelsize": 6, "axes.titlesize": 6,
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "lines.linewidth": 0.5,
})
MM = 1/25.4

# ── Drugs and metadata ──────────────────────────────────────────────────────
EXCLUDE = {"MCAM", "FSK"}
GROUPS = {
    "strong":       ["DAMGO","Carfentanil","Fentanyl","SR17018","Morphine","C6guano"],
    "intermediate": ["Methadone","TRV130","PZM21"],
    "weak":         ["MP","Buprenorphine","Nalbuphine","Butorphanol"],
    "antagonist":   ["Naloxone","Naltrexone"],
}
DRUG_CLS = {d: c for c, ds in GROUPS.items() for d in ds}
CLS_COL = {"strong":"#231f20","intermediate":"#9b1c1d",
           "weak":"#3779b9","antagonist":"#6f08a3"}

# ── Load data ────────────────────────────────────────────────────────────────
dms = pd.read_csv(ROOT / "dms_scores/composite_dms_scores.csv")
gpcrdb = pd.read_csv(ROOT / "annotations/GPCRdb_OPRM1_table.csv").rename(
    columns={"pos":"position"})

# Identify effect columns, excluding MCAM/FSK
effect_cols = [c for c in dms.columns
               if c.endswith("_effect")
               and not any(x in c for x in EXCLUDE)
               and "Surface" not in c
               and "pos_effect" not in c]
drugs = [c.replace("_effect", "") for c in effect_cols]
drugs = [d for d in drugs if d in DRUG_CLS]
effect_cols = [f"{d}_effect" for d in drugs]
print(f"Drugs used ({len(drugs)}): {drugs}")

POS_MIN, POS_MAX = 65, 355
mis = dms[(dms.type == "missense") &
          (dms.position >= POS_MIN) &
          (dms.position <= POS_MAX)]
print(f"Filtered missense rows (positions {POS_MIN}-{POS_MAX}): {len(mis)}")

# ── Per-position mean missense effect matrix ───────────────────────────────
pos_means = mis.groupby("position")[effect_cols].mean().dropna()
print(f"Positions with complete data across all {len(drugs)} drugs: {len(pos_means)}")

X = pos_means.values.astype(float)
# Standardize (z-score each drug column)
mu = X.mean(axis=0); sd = X.std(axis=0); sd[sd==0] = 1
Xz = (X - mu) / sd

# PCA via SVD
U, S, Vt = np.linalg.svd(Xz - Xz.mean(0), full_matrices=False)
scores = U * S
var_explained = (S**2) / (S**2).sum()
print(f"Variance explained (PC1-5): {var_explained[:5].round(3)}")

# Loadings: drug contribution to each PC
loadings = Vt.T * S
loadings_df = pd.DataFrame(loadings[:, :5],
                             index=drugs,
                             columns=[f"PC{i+1}" for i in range(5)])
loadings_df["class"] = [DRUG_CLS[d] for d in drugs]
loadings_df.to_csv(OUT / "pca_loadings_table.csv")

# Position scores
pos_means_out = pos_means.copy()
for i in range(5):
    pos_means_out[f"PC{i+1}"] = scores[:, i]
pos_means_out = pos_means_out.reset_index().merge(
    gpcrdb[["position","SSE"]], on="position", how="left")
pos_means_out.to_csv(OUT / "pca_position_scores.csv", index=False)

# ════════════════════════════════════════════════════════════════════════════
# FIGURES
# ════════════════════════════════════════════════════════════════════════════

# ── Scree plot ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(80*MM, 50*MM))
cum = np.cumsum(var_explained)
ax.bar(range(1, len(var_explained)+1), var_explained*100,
       color="#4393c3", edgecolor="black", linewidth=0.3)
ax.plot(range(1, len(var_explained)+1), cum*100,
        color="#b2182b", marker="o", markersize=3, lw=1)
ax.set_xlabel("PC")
ax.set_ylabel("% variance")
ax.set_title("Position-level PCA: variance explained",
             fontweight="bold", fontsize=7, loc="left")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
for i, v in enumerate(var_explained):
    if v > 0.01:
        ax.text(i+1, v*100 + 1, f"{v*100:.0f}%", ha="center", fontsize=5)
fig.tight_layout()
fig.savefig(OUT / "scree_position.pdf", dpi=600, bbox_inches="tight")
plt.close()
print("Saved: scree_position.pdf")

# ── PC1 vs PC2 scatter colored by SSE ───────────────────────────────────────
fig, ax = plt.subplots(figsize=(95*MM, 85*MM))

SSE_COL = {
    "TM1":"#e41a1c","TM2":"#377eb8","TM3":"#4daf4a",
    "TM4":"#984ea3","TM5":"#ff7f00","TM6":"#a65628",
    "TM7":"#f781bf","H8":"#999999",
    "ICL1":"#cccccc","ICL2":"#cccccc","ICL3":"#cccccc",
    "ECL1":"#eeeeee","ECL2":"#eeeeee","ECL3":"#eeeeee",
}

for sse in sorted(pos_means_out.SSE.dropna().unique()):
    sub = pos_means_out[pos_means_out.SSE == sse]
    ax.scatter(sub.PC1, sub.PC2, s=15,
               color=SSE_COL.get(sse, "#cccccc"),
               edgecolors="black", linewidths=0.2,
               alpha=0.85, label=sse)

KEY = {119:"A119", 148:"I148", 152:"N152", 280:"I280", 77:"Y77",
       341:"L341", 287:"V287", 175:"V175", 149:"D149", 151:"Y151",
       167:"R167", 297:"P297"}
for pos, lbl in KEY.items():
    sub = pos_means_out[pos_means_out.position == pos]
    if len(sub) == 0: continue
    r = sub.iloc[0]
    ax.annotate(lbl, (r.PC1, r.PC2), fontsize=5, fontweight="bold",
                xytext=(3, 3), textcoords="offset points")

ax.axhline(0, color="grey", lw=0.3, ls="--")
ax.axvline(0, color="grey", lw=0.3, ls="--")
ax.set_xlabel(f"PC1 ({var_explained[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({var_explained[1]*100:.1f}%)")
ax.set_title(f"Position-level PCA (n={len(pos_means_out)} positions)",
             fontweight="bold", fontsize=7, loc="left")
ax.legend(fontsize=4.5, frameon=False, loc="best", ncol=2)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "pca_position_scatter.pdf", dpi=600, bbox_inches="tight")
plt.close()
print("Saved: pca_position_scatter.pdf")

# ── Drug loadings in PC1/PC2 space ──────────────────────────────────────────
fig, ax = plt.subplots(figsize=(90*MM, 90*MM))
for d in drugs:
    x, y = loadings_df.loc[d, "PC1"], loadings_df.loc[d, "PC2"]
    ax.annotate("", xy=(x, y), xytext=(0, 0),
                arrowprops=dict(arrowstyle="->",
                                color=CLS_COL[DRUG_CLS[d]], lw=0.8))
    ax.text(x*1.05, y*1.05, d, fontsize=5, color=CLS_COL[DRUG_CLS[d]],
            ha="center" if abs(x) < 0.1 else ("left" if x > 0 else "right"))

ax.axhline(0, color="grey", lw=0.3, ls="--")
ax.axvline(0, color="grey", lw=0.3, ls="--")
lim = max(abs(loadings_df[["PC1","PC2"]].values).max() * 1.15, 1)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.set_aspect("equal")
ax.set_xlabel(f"PC1 loading ({var_explained[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 loading ({var_explained[1]*100:.1f}%)")
ax.set_title("Drug loadings on position-level PCA",
             fontweight="bold", fontsize=7, loc="left")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

handles = [mpatches.Patch(color=CLS_COL[c], label=c) for c in CLS_COL]
ax.legend(handles=handles, fontsize=5, frameon=False, loc="best")
fig.tight_layout()
fig.savefig(OUT / "pca_drug_loadings.pdf", dpi=600, bbox_inches="tight")
plt.close()
print("Saved: pca_drug_loadings.pdf")

# ── PC1/PC2/PC3 along receptor sequence ─────────────────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(170*MM, 70*MM), sharex=True)

sse_ranges = {"TM1":(67,100),"TM2":(105,130),"TM3":(145,170),"TM4":(174,200),
              "TM5":(218,260),"TM6":(267,300),"TM7":(312,340),"H8":(343,358)}
sse_colors_bg = {"TM1":"#fee0d2","TM2":"#deebf7","TM3":"#e5f5e0",
                 "TM4":"#f2e6ff","TM5":"#fff3e0","TM6":"#fce4ec",
                 "TM7":"#fff9c4","H8":"#f5f5f5"}

pos_sorted = pos_means_out.sort_values("position")
for ax, pc in zip(axes, ["PC1","PC2","PC3"]):
    for sse, (s, e) in sse_ranges.items():
        ax.axvspan(s, e, color=sse_colors_bg[sse], alpha=0.4, zorder=0)
    colors = ["#b2182b" if v > 0 else "#2166ac" for v in pos_sorted[pc]]
    ax.bar(pos_sorted.position, pos_sorted[pc], color=colors,
           edgecolor="none", width=1)
    ax.axhline(0, color="black", lw=0.3)
    ax.set_ylabel(pc, fontsize=6)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    if pc == "PC1":
        for p, lbl in KEY.items():
            r = pos_sorted[pos_sorted.position == p]
            if len(r):
                r = r.iloc[0]
                ax.text(r.position, r[pc] + np.sign(r[pc])*0.3,
                        lbl, fontsize=4, ha="center", fontweight="bold")
        for sse, (s, e) in sse_ranges.items():
            ax.text((s+e)/2, ax.get_ylim()[1]*0.92, sse, fontsize=4,
                    ha="center", color="grey")

axes[-1].set_xlabel("Receptor position")
axes[-1].set_xlim(POS_MIN, POS_MAX)
fig.suptitle("Position-level PC scores along receptor sequence",
             fontweight="bold", fontsize=7, y=0.99)
fig.tight_layout()
fig.savefig(OUT / "pca_receptor_map.pdf", dpi=600, bbox_inches="tight")
plt.close()
print("Saved: pca_receptor_map.pdf")

print("\nAll position-level PCA outputs saved to plots/pca/position/")

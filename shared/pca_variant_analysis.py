#!/usr/bin/env python3
"""
Variant-level PCA on per-VARIANT missense effects across 15 ligands (MCAM & FSK
excluded), positions 65-355. Companion to ../position/pca_analysis.py, which runs
the same recipe on per-POSITION mean effects.

Each row = one missense variant; each column = that variant's effect for one ligand.

PC signs are aligned to the position-level PCA (../position/pca_position_scores.csv)
so the two PCAs share a canonical orientation. Downstream scripts
(build_pc2_variant.py, pca_three_panel_variant*.py) re-align defensively, so the
absolute sign here does not affect any figure.

Outputs (plots/pca/variant/):
  pca_variant_scores.csv   (hgvs, position, wildtype, mutation, <15>_effect, PC1-3)
  scree_variant.pdf
  pca_variant_scatter.pdf
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "plots" / "pca" / "variant"
POS = ROOT / "plots" / "pca" / "position"
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

EXCLUDE = {"MCAM", "FSK"}
POS_MIN, POS_MAX = 65, 355

# ── Load data ────────────────────────────────────────────────────────────────
dms = pd.read_csv(ROOT / "dms_scores/composite_dms_scores.csv")

# 15 ligand effect columns, sorted by drug name (matches the position PCA / prior file)
effect_cols = [c for c in dms.columns
               if c.endswith("_effect")
               and not any(x in c for x in EXCLUDE)
               and "Surface" not in c
               and "pos_effect" not in c]
drugs = sorted(c[:-len("_effect")] for c in effect_cols)
effect_cols = [f"{d}_effect" for d in drugs]
print(f"Drugs used ({len(drugs)}): {drugs}")

mis = dms[(dms.type == "missense") &
          (dms.position >= POS_MIN) &
          (dms.position <= POS_MAX)].dropna(subset=effect_cols).copy()
print(f"Variants with complete data across {len(drugs)} drugs "
      f"(positions {POS_MIN}-{POS_MAX}): {len(mis)}")

# ── PCA via SVD (z-score each drug column, mirror position-level recipe) ──────
X = mis[effect_cols].values.astype(float)
mu = X.mean(axis=0); sd = X.std(axis=0); sd[sd == 0] = 1
Xz = (X - mu) / sd
U, S, Vt = np.linalg.svd(Xz - Xz.mean(0), full_matrices=False)
scores = U * S
var_explained = (S**2) / (S**2).sum()
print(f"Variance explained (PC1-5): {var_explained[:5].round(3)}")

for i in range(3):
    mis[f"PC{i+1}"] = scores[:, i]

# ── Align PC1-3 signs to the position-level PCA (canonical orientation) ───────
pos_path = POS / "pca_position_scores.csv"
if pos_path.exists():
    posPCA = pd.read_csv(pos_path).set_index("position")
    vm = mis.groupby("position")[["PC1", "PC2", "PC3"]].mean()
    for pc in ["PC1", "PC2", "PC3"]:
        common = vm.index.intersection(posPCA.index)
        r = np.corrcoef(posPCA.loc[common, pc].values, vm.loc[common, pc].values)[0, 1]
        if r < 0:
            mis[pc] = -mis[pc]
            print(f"  {pc}: flipped to align with position PCA")
        else:
            print(f"  {pc}: already aligned with position PCA")
else:
    print("WARN: position scores not found; PC signs left as raw SVD output")

# ── Write scores table ───────────────────────────────────────────────────────
cols = ["hgvs", "position", "wildtype", "mutation"] + effect_cols + ["PC1", "PC2", "PC3"]
mis[cols].to_csv(OUT / "pca_variant_scores.csv", index=False)
print(f"Wrote pca_variant_scores.csv ({len(mis)} variants)")

# ── Scree ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(80*MM, 50*MM))
ax.bar(range(1, len(var_explained)+1), var_explained*100,
       color="#4393c3", edgecolor="black", linewidth=0.3)
ax.plot(range(1, len(var_explained)+1), np.cumsum(var_explained)*100,
        color="#b2182b", marker="o", markersize=3, lw=1)
ax.set_xlabel("PC"); ax.set_ylabel("% variance")
ax.set_title("Variant-level PCA: variance explained",
             fontweight="bold", fontsize=7, loc="left")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
for i, v in enumerate(var_explained):
    if v > 0.01:
        ax.text(i+1, v*100 + 1, f"{v*100:.0f}%", ha="center", fontsize=5)
fig.tight_layout()
fig.savefig(OUT / "scree_variant.pdf", dpi=600, bbox_inches="tight")
plt.close()
print("Saved: scree_variant.pdf")

# ── PC1 vs PC2 variant scatter ───────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(85*MM, 85*MM))
ax.scatter(mis.PC1, mis.PC2, s=3, color="#4d4d4d", edgecolors="none", alpha=0.4)
ax.axhline(0, color="grey", lw=0.3, ls="--"); ax.axvline(0, color="grey", lw=0.3, ls="--")
ax.set_xlabel(f"PC1 ({var_explained[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({var_explained[1]*100:.1f}%)")
ax.set_title(f"Variant-level PCA (n={len(mis)} variants)",
             fontweight="bold", fontsize=7, loc="left")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "pca_variant_scatter.pdf", dpi=600, bbox_inches="tight")
plt.close()
print("Saved: pca_variant_scatter.pdf")
print("\nAll variant-level PCA outputs saved to plots/pca/variant/")

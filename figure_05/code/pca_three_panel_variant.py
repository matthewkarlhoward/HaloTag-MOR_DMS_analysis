#!/usr/bin/env python3
"""
Three-panel VARIANT-level PCA summary. One figure, exactly 130 x 50 mm, three
equal-size panels. All three panels come from ONE variant PCA fit, so PC2's
variance % (10.2%) is consistent across the figure (panel C no longer borrows the
position PCA's 15.6% loadings).

  A  variant PC1 vs PC2 biplot
  B  variant PC1 score vs that variant's surface-expression score
  C  variant PC2 loading vs TruPath Gi1 Emax (efficacy-class coloured)

Variant PCA: PC1 49.1%, PC2 10.2%. Scores + loadings are sign-aligned to the
position PCA (PC1 already matches; PC2 is flipped) so orientation is consistent.
Output: plots/pca/pca_three_panel_variant.pdf
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")  # exact figure size (macosx backend snaps figsize to 2 dp)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "plots" / "pca" / "variant"
MM = 1/25.4

plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 6,
    "mathtext.default": "regular", "text.color": "black",
    "axes.labelsize": 6, "axes.labelcolor": "black", "axes.edgecolor": "black",
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "xtick.color": "black", "ytick.color": "black", "legend.fontsize": 6,
    "axes.linewidth": 0.5, "lines.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
})

# ── Palettes ──────────────────────────────────────────────────────────────────
MOTIFS = {"DRY": [166,167,168], "CWxP": [294,295,296,297],
          "NPxxY": [334,335,336,337,338], "PIF": [157,246,291], "Na+": [116]}
MOTIF_COL = {"DRY":"#9b1c1d","CWxP":"#0e6d39","NPxxY":"#3779b9",
             "PIF":"#e66101","Na+":"#6f08a3"}
MOTIF_LBL = {"DRY":"DRY","CWxP":"CWxP","NPxxY":"NPxxY","PIF":"PIF","Na+":r"Na$^+$"}
pos_to_motif = {p: m for m, ps in MOTIFS.items() for p in ps}
CLS_COL = {"strong":"#231f20","intermediate":"#9b1c1d",
           "weak":"#3779b9","antagonist":"#6f08a3"}

# ── Data: ONE variant PCA fit -> scores (A/B) AND loadings (C) ────────────────
# Recomputing the variant PCA here (z-score per drug + SVD, same recipe as the
# position PCA) gives the variant loadings for panel C, so PC2's variance % is the
# SAME (10.2%) across all three panels instead of borrowing the position PCA's
# 15.6% loadings. Sign-aligned to the position PCA so orientation matches.
dms = pd.read_csv(ROOT / "dms_scores/composite_dms_scores.csv")
vfull = pd.read_csv(OUT / "pca_variant_scores.csv")
ecols = [c for c in vfull.columns if c.endswith("_effect")]
Xv = vfull[ecols].to_numpy(float)
keep = ~np.isnan(Xv).any(axis=1)
var = vfull.loc[keep, ["hgvs", "position"]].reset_index(drop=True)
Xz = Xv[keep]
mu = Xz.mean(0); sd = Xz.std(0); sd[sd == 0] = 1
Xz = (Xz - mu) / sd
U, S, Vt = np.linalg.svd(Xz - Xz.mean(0), full_matrices=False)
varpct = (S**2) / (S**2).sum() * 100                  # PC1 49.1%, PC2 10.2%
var["PC1"] = (U * S)[:, 0]; var["PC2"] = (U * S)[:, 1]
drugs = [c[:-7] for c in ecols]
ld = pd.DataFrame((Vt.T * S)[:, :2], index=drugs, columns=["PC1", "PC2"])

# Sign-align PC1/PC2 (scores + matching loadings) to the position-level PCA
posPCA = pd.read_csv(ROOT / "plots" / "pca" / "position" / "pca_position_scores.csv")[["position", "PC1", "PC2"]]
for col in ("PC1", "PC2"):
    vm = var.groupby("position")[col].mean().rename("v").reset_index()
    mm = posPCA[["position", col]].merge(vm, on="position").dropna()
    if np.corrcoef(mm[col], mm["v"])[0, 1] < 0:
        var[col] *= -1; ld[col] *= -1
PC1_PCT, PC2_PCT = varpct[0], varpct[1]

# ── Figure: 3 equal panels, exactly 130 x 50 mm ───────────────────────────────
fig, (axA, axB, axC) = plt.subplots(1, 3)
fig.set_size_inches(130*MM, 50*MM)
fig.subplots_adjust(left=0.07, right=0.99, bottom=0.18, top=0.95, wspace=0.40)

# ===== Panel A: variant PC1 vs PC2 biplot (motifs) ===========================
axA.scatter(var.PC1, var.PC2, s=3, color="black", alpha=0.25,
            edgecolors="none", zorder=1)
axA.axhline(0, color="grey", lw=0.5, ls="--")
axA.axvline(0, color="grey", lw=0.5, ls="--")
axA.set_xlim(var.PC1.min()-1.0, var.PC1.max()+2.0)
axA.set_ylim(var.PC2.min()-0.6, var.PC2.max()+0.8)
axA.set_xlabel(f"PC1 ({PC1_PCT:.1f}%)", labelpad=1)
axA.set_ylabel(f"PC2 ({PC2_PCT:.1f}%)", labelpad=1)
axA.tick_params(axis="both", pad=1)
axA.spines["top"].set_visible(False); axA.spines["right"].set_visible(False)

# ===== Panel B: variant PC1 vs surface (motifs) ==============================
dB = var.merge(dms[["hgvs", "Surface_effect"]], on="hgvs").rename(
    columns={"Surface_effect": "surface"}).dropna(subset=["PC1", "surface"])
rB, _ = pearsonr(dB.surface, dB.PC1); rhoB, _ = spearmanr(dB.surface, dB.PC1)
axB.axhline(0, color="grey", lw=0.5, ls="--")
axB.axvline(0, color="grey", lw=0.5, ls="--")
xs = np.linspace(dB.surface.min(), dB.surface.max(), 100)
mB, bB = np.polyfit(dB.surface, dB.PC1, 1)
axB.scatter(dB.surface, dB.PC1, s=3, color="black", alpha=0.25,
            edgecolors="none", zorder=1)
axB.plot(xs, mB*xs + bB, color="0.45", lw=0.8, zorder=2)
axB.text(0.03, 0.04, f"ρ = {rhoB:.2f}", transform=axB.transAxes,
         ha="left", va="bottom",
         bbox=dict(facecolor="white", edgecolor="none", pad=0.5))
axB.set_xlabel("Surface expression score", labelpad=1)
axB.set_ylabel(f"PC1 ({PC1_PCT:.1f}%)", labelpad=1)
axB.tick_params(axis="both", pad=1)
axB.spines["top"].set_visible(False); axB.spines["right"].set_visible(False)

# ===== Panel C: variant PC2 loading vs TruPath Emax (same variant PCA) ========
cls = pd.read_csv(ROOT / "plots" / "pca" / "position" / "pca_loadings_table.csv", index_col=0)["class"]  # drug->class
ldC = ld.reset_index().rename(columns={"index": "drug"})
ldC["class"] = ldC.drug.map(cls)
ldC["key"] = ldC.drug.str.lower()
eff = pd.read_excel(
    ROOT / "trupath_vs_camp_efficacy/WT_cAMP_vs_TRUPATH_All_Drugs.xlsx",
    sheet_name="cAMP vs TRUPATH (WT)").rename(columns={"TRUPATH_Gi1_Emax_pct_DAMGO": "emax"})
eff["key"] = eff.Ligand.str.lower()
dC = ldC.merge(eff[["key", "emax"]], on="key").dropna(subset=["PC2", "emax"])
xC, yC = dC.emax.to_numpy(), dC.PC2.to_numpy()
rhoC, _ = spearmanr(xC, yC)
axC.axhline(0, color="grey", lw=0.5, ls="--", zorder=0)
xs = np.linspace(xC.min(), xC.max(), 100)
mC, bC = np.polyfit(xC, yC, 1)
axC.plot(xs, mC*xs + bC, color="0.45", lw=0.5, zorder=1)
axC.scatter(xC, yC, s=12, c=dC["class"].map(CLS_COL).to_numpy(),
            edgecolors="none", zorder=3)
axC.text(0.97, 0.97, f"ρ = {rhoC:.2f}", transform=axC.transAxes,
         ha="right", va="top")
axC.set_xticks([25, 50, 75, 100])
axC.set_xlabel(r"TruPath Gi1 E$_{max}$", labelpad=1)
axC.set_ylabel(f"PC2 loading ({PC2_PCT:.1f}%)", labelpad=1)
axC.tick_params(axis="both", pad=1)
axC.spines["top"].set_visible(False); axC.spines["right"].set_visible(False)
hC = [mpatches.Patch(color=CLS_COL[c], label=c.capitalize())
      for c in ["strong", "intermediate", "weak", "antagonist"]]
axC.legend(handles=hC, frameon=False, loc="lower left", handlelength=1.0,
           handletextpad=0.4, labelspacing=0.2, borderpad=0.15, borderaxespad=0.25)

fig.savefig(OUT / "pca_three_panel_variant.pdf")
plt.close()
print(f"A/B variant (n={len(dB)}) | B PC1~surf r={rB:.2f} | C PC2~Emax rho={rhoC:.2f}")
print("Saved: pca_three_panel_variant.pdf")

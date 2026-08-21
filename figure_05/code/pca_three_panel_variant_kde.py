#!/usr/bin/env python3
"""
VARIANT-level three-panel PCA summary, panels A & B drawn as a 2D Gaussian-KDE of
the ~5.5k variants instead of points, in two styles:

  pca_three_panel_variant_contour.pdf  — grey->black filled contour bands (white bg)
  pca_three_panel_variant_density.pdf  — inferno density faded onto a white bg

Panel C (13 ligand loadings) stays as points. All from ONE variant PCA fit so PC2 =
10.2% throughout. One figure, exactly 130 x 50 mm, three equal panels.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")  # exact figure size
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import spearmanr, gaussian_kde

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
CLS_COL = {"strong":"#231f20","intermediate":"#9b1c1d",
           "weak":"#3779b9","antagonist":"#6f08a3"}

# ── Variant PCA (one fit -> scores for A/B, loadings for C) ────────────────────
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
varpct = (S**2) / (S**2).sum() * 100
var["PC1"] = (U*S)[:, 0]; var["PC2"] = (U*S)[:, 1]
drugs = [c[:-7] for c in ecols]
ld = pd.DataFrame((Vt.T*S)[:, :2], index=drugs, columns=["PC1", "PC2"])
posPCA = pd.read_csv(ROOT / "plots" / "pca" / "position" / "pca_position_scores.csv")[["position", "PC1", "PC2"]]
for col in ("PC1", "PC2"):
    vm = var.groupby("position")[col].mean().rename("v").reset_index()
    mm = posPCA[["position", col]].merge(vm, on="position").dropna()
    if np.corrcoef(mm[col], mm["v"])[0, 1] < 0:
        var[col] *= -1; ld[col] *= -1
PC1_PCT, PC2_PCT = varpct[0], varpct[1]

# Panel B data (variant PC1 vs per-variant surface)
dB = var.merge(dms[["hgvs", "Surface_effect"]], on="hgvs").rename(
    columns={"Surface_effect": "surface"}).dropna(subset=["PC1", "surface"])
rhoB = spearmanr(dB.surface, dB.PC1)[0]

# Panel C data (variant PC2 loading vs Emax)
cls = pd.read_csv(ROOT / "plots" / "pca" / "position" / "pca_loadings_table.csv", index_col=0)["class"]
ldC = ld.reset_index().rename(columns={"index": "drug"})
ldC["class"] = ldC.drug.map(cls); ldC["key"] = ldC.drug.str.lower()
eff = pd.read_excel(
    ROOT / "trupath_vs_camp_efficacy/WT_cAMP_vs_TRUPATH_All_Drugs.xlsx",
    sheet_name="cAMP vs TRUPATH (WT)").rename(columns={"cAMP Response": "emax"})
eff["key"] = eff.Ligand.str.lower()
dC = ldC.merge(eff[["key", "emax"]], on="key").dropna(subset=["PC2", "emax"])
xC, yC = dC.emax.to_numpy(), dC.PC2.to_numpy()
rhoC = spearmanr(xC, yC)[0]


def kde_grid(x, y, n=130, pad=0.06):
    xmin, xmax = x.min(), x.max(); ymin, ymax = y.min(), y.max()
    dx, dy = (xmax-xmin)*pad, (ymax-ymin)*pad
    xmin -= dx; xmax += dx; ymin -= dy; ymax += dy
    xx, yy = np.mgrid[xmin:xmax:n*1j, ymin:ymax:n*1j]
    zz = gaussian_kde(np.vstack([x, y]))(np.vstack([xx.ravel(), yy.ravel()]))
    return xx, yy, zz.reshape(xx.shape), (xmin, xmax, ymin, ymax)


def density_panel(ax, x, y, mode):
    xx, yy, zz, ext = kde_grid(x, y)
    if mode == "contour":                       # grey -> black filled bands, white bg
        levels = np.linspace(zz.max()*0.04, zz.max(), 8)
        colors = plt.cm.Greys(np.linspace(0.12, 0.92, len(levels)-1))
        ax.contourf(xx, yy, zz, levels=levels, colors=colors, zorder=1)
    else:                                       # inferno density faded onto white
        norm = zz / zz.max()
        rgba = mpl.cm.inferno(norm)
        rgba[..., 3] = np.clip(norm*3.5, 0, 1)
        ax.imshow(np.transpose(rgba, (1, 0, 2)), origin="lower", extent=ext,
                  aspect="auto", interpolation="bilinear", zorder=1)
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])


def build(mode):
    fig, (axA, axB, axC) = plt.subplots(1, 3)
    fig.set_size_inches(130*MM, 50*MM)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.18, top=0.95, wspace=0.40)

    # Panel A: variant PC1 vs PC2 density
    density_panel(axA, var.PC1.values, var.PC2.values, mode)
    axA.axhline(0, color="grey", lw=0.5, ls="--", alpha=0.6, zorder=2)
    axA.axvline(0, color="grey", lw=0.5, ls="--", alpha=0.6, zorder=2)
    axA.set_xlabel(f"PC1 ({PC1_PCT:.1f}%)", labelpad=1)
    axA.set_ylabel(f"PC2 ({PC2_PCT:.1f}%)", labelpad=1)
    axA.tick_params(axis="both", pad=1)
    axA.spines["top"].set_visible(False); axA.spines["right"].set_visible(False)

    # Panel B: variant PC1 vs surface density + trend + rho (bottom-left)
    density_panel(axB, dB.surface.values, dB.PC1.values, mode)
    axB.axhline(0, color="grey", lw=0.5, ls="--", alpha=0.6, zorder=2)
    axB.axvline(0, color="grey", lw=0.5, ls="--", alpha=0.6, zorder=2)
    xs = np.linspace(dB.surface.min(), dB.surface.max(), 100)
    mB, bB = np.polyfit(dB.surface, dB.PC1, 1)
    axB.plot(xs, mB*xs + bB, color="0.25", lw=0.8, zorder=3)
    axB.text(0.03, 0.04, f"ρ = {rhoB:.2f}", transform=axB.transAxes,
             ha="left", va="bottom",
             bbox=dict(facecolor="white", edgecolor="none", pad=0.5))
    axB.set_xlabel("Surface expression score", labelpad=1)
    axB.set_ylabel(f"PC1 ({PC1_PCT:.1f}%)", labelpad=1)
    axB.tick_params(axis="both", pad=1)
    axB.spines["top"].set_visible(False); axB.spines["right"].set_visible(False)

    # Panel C: variant PC2 loading vs Emax (points, unchanged)
    axC.axhline(0, color="grey", lw=0.5, ls="--", zorder=0)
    xsC = np.linspace(xC.min(), xC.max(), 100)
    mC, bC = np.polyfit(xC, yC, 1)
    axC.plot(xsC, mC*xsC + bC, color="0.45", lw=0.5, zorder=1)
    axC.scatter(xC, yC, s=12, c=dC["class"].map(CLS_COL).to_numpy(),
                edgecolors="none", zorder=3)
    if mode == "contour":      # rho sits just above the ligand class legend
        axC.text(0.04, 0.31, f"ρ = {rhoC:.2f}", transform=axC.transAxes,
                 ha="left", va="bottom")
    else:
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
               handletextpad=0.4, labelspacing=0.2, borderpad=0.15,
               borderaxespad=0.25)

    out = OUT / f"pca_three_panel_variant_{mode}.pdf"
    fig.savefig(out); plt.close()
    print(f"Saved: {out.name}")


for mode in ("contour", "density"):
    build(mode)

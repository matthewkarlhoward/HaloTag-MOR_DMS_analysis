#!/usr/bin/env python3
"""
Scatter: PC2 drug loading (y) vs TruPath Gi1 efficacy (x), one dot per drug.
Tests whether the PCA's PC2 axis tracks measured G-protein efficacy.

TruPath efficacy source: trupath_vs_camp_efficacy/WT_cAMP_vs_TRUPATH_All_Drugs.xlsx.
NOTE: that file's headers are SWAPPED (see camp_vs_trupath_combined.R) — the column
labeled 'cAMP Response' is actually TruPath Gi1 Emax (% DAMGO). We use that column.

13 drugs (Naltrexone & SR17018 have no TruPath value).
Output: plots/pca/pca_pc2_vs_trupath_efficacy.pdf
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")  # exact figure size (macosx backend snaps figsize to 2 dp)
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "plots" / "pca" / "position"
MM = 1/25.4

# Efficacy-class palette (matches radars / stacked heatmaps)
CLS_COL = {"strong": "#231f20", "intermediate": "#9b1c1d",
           "weak": "#3779b9", "antagonist": "#6f08a3"}

# House style: black text, Helvetica, size 6, 0.5 pt lines
plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 6,
    "mathtext.default": "regular",  # keep subscripts upright in Helvetica
    "text.color": "black",
    "axes.labelsize": 6, "axes.labelcolor": "black", "axes.edgecolor": "black",
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "xtick.color": "black", "ytick.color": "black",
    "legend.fontsize": 6,
    "axes.linewidth": 0.5, "lines.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
})

# ── Load PC2 loadings + TruPath efficacy ─────────────────────────────────────
load = pd.read_csv(OUT / "pca_loadings_table.csv", index_col=0)
load = load.reset_index().rename(columns={"index": "drug"})
load["key"] = load["drug"].str.lower()

eff = pd.read_excel(
    ROOT / "trupath_vs_camp_efficacy/WT_cAMP_vs_TRUPATH_All_Drugs.xlsx",
    sheet_name="cAMP vs TRUPATH (WT)")
# Undo the swapped headers: 'cAMP Response' is really TruPath Gi1 Emax (% DAMGO)
eff = eff.rename(columns={"cAMP Response": "trupath_emax"})
eff["key"] = eff["Ligand"].str.lower()

df = load.merge(eff[["key", "trupath_emax"]], on="key").dropna(
    subset=["PC2", "trupath_emax"])
print(f"Drugs in scatter (n={len(df)}): {sorted(df.drug)}")

x = df["trupath_emax"].to_numpy()
y = df["PC2"].to_numpy()
r, p = pearsonr(x, y)
rho, p_s = spearmanr(x, y)
print(f"Pearson r={r:.2f} p={p:.1e} | Spearman rho={rho:.2f} p={p_s:.1e}")

# ── Figure ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(layout="constrained")
fig.set_size_inches(33*MM, 33*MM)            # exact 33 x 33 mm
fig.get_layout_engine().set(w_pad=0.01, h_pad=0.01, wspace=0, hspace=0)

ax.axhline(0, color="grey", lw=0.5, ls="--", zorder=0)

# Linear trend guide
xs = np.linspace(x.min(), x.max(), 100)
m, b = np.polyfit(x, y, 1)
ax.plot(xs, m*xs + b, color="grey", lw=0.5, zorder=1)

# Points coloured by efficacy class (radar / heatmap palette)
ax.scatter(x, y, s=12, c=df["class"].map(CLS_COL).to_numpy(),
           edgecolors="none", zorder=3)

# Spearman rho only (label omitted)
ax.text(0.97, 0.97, f"ρ = {rho:.2f}",
        transform=ax.transAxes, ha="right", va="top", fontsize=6, color="black")

ax.set_xlabel(r"TruPath Gi1 E$_{max}$", labelpad=1)
ax.set_ylabel("PC2 loading (15.6%)", labelpad=1)
ax.set_xticks([25, 50, 75, 100])
ax.tick_params(axis="both", pad=1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.savefig(OUT / "pca_pc2_vs_trupath_efficacy.pdf", dpi=600)
plt.close()
print("Saved: pca_pc2_vs_trupath_efficacy.pdf")

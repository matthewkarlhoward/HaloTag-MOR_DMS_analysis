#!/usr/bin/env python3
"""
Surface expression vs morphine potency and efficacy, as a 2 x 2 grid.

           EC50                     Emax
  row 1    variant level            variant level        (density-coloured)
  row 2    position level           position level       (position means)

x = Surface_effect (negative = the mutation reduces surface expression); the
per-position value is the mean over that position's missense variants, matching
plots/pca/variant/pca_pc1_vs_surface.py.

y comes from the robust 3-parameter refit (curve_refitting/refit_3param_robust_
{drug}.csv), the same source as the variant_morphine_vs_fentanyl_* panels:
  EC50  fit_ok, missense, new_curve_type == "sigmoid"
  Emax  fit_ok, missense, new_curve_type in {sigmoid, flat, no_baseline},
        rescaled to the Activity axis (synonymous-sigmoid Emin -> 0, Emax -> 1)
The two columns therefore carry different variant counts on purpose: EC50 is
only defined for a fitted sigmoid, while a flat curve still has a real top.

Dashed references: synonymous-mean surface effect (vertical), synonymous-mean
EC50 / Activity = 1 (horizontal). Spearman rho and n are printed in each panel.

Usage:  python3 surface_vs_ec50_emax.py [Drug]      (default Morphine)
Output: surface_vs_{drug}_ec50_emax.pdf  (+ _stats.csv)
"""
import matplotlib
matplotlib.use("Agg")

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy.stats import gaussian_kde, pearsonr, spearmanr
from statsmodels.nonparametric.smoothers_lowess import lowess

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
DRUG = sys.argv[1] if len(sys.argv) > 1 else "Morphine"

plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 6,
    "mathtext.default": "regular",
    "axes.labelsize": 6, "axes.titlesize": 6,
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "text.color": "black", "axes.edgecolor": "black",
    "axes.labelcolor": "black", "axes.titlecolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "lines.linewidth": 0.5, "patch.linewidth": 0.5,
    "pdf.fonttype": 42,
})
MM = 1 / 25.4
EMAX_OK = {"sigmoid", "flat", "no_baseline"}
LOESS_COL = "#9b1c1d"
LOESS_FRAC = 0.3      # neighbourhood fraction; 0.3 tracks curvature without chasing noise
LOESS_IT = 3          # robustifying iterations, so outliers do not drag the fit

# ── Data ─────────────────────────────────────────────────────────────────────
dms = pd.read_csv(ROOT / "dms_scores" / "composite_dms_scores.csv")
fit = pd.read_csv(ROOT / "curve_refitting" / f"refit_3param_robust_{DRUG.lower()}.csv")

surf = dms[["hgvs", "Surface_effect"]].dropna().rename(
    columns={"Surface_effect": "surface"})
syn_surface = dms.loc[dms.type == "synonymous", "Surface_effect"].mean()

syn_fit = fit[(fit.type == "synonymous") & (fit.new_curve_type == "sigmoid")]
e0, e1 = syn_fit.fitted_emin.mean(), syn_fit.fitted_emax.mean()
syn_ec50 = syn_fit.fitted_ec50_logM.mean()

ec50 = (fit[fit.fit_ok & (fit.type == "missense") &
            (fit.new_curve_type == "sigmoid")]
        [["hgvs", "position", "fitted_ec50_logM"]]
        .rename(columns={"fitted_ec50_logM": "y"})
        .merge(surf, on="hgvs"))

emax = fit[fit.fit_ok & (fit.type == "missense") &
           fit.new_curve_type.isin(EMAX_OK)][["hgvs", "position", "fitted_emax"]].copy()
emax["y"] = (emax.fitted_emax - e0) / (e1 - e0)
emax = emax.drop(columns="fitted_emax").merge(surf, on="hgvs")

COLUMNS = [
    dict(key="ec50", df=ec50, ylab=f"EC$_{{50}}$ {DRUG} (log[M])", href=syn_ec50),
    dict(key="emax", df=emax, ylab=f"E$_{{max}}$ {DRUG} (activity)", href=1.0),
]

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2)
fig.set_size_inches(90 * MM, 85 * MM)

stats = []
for col, spec in enumerate(COLUMNS):
    var = spec["df"]
    pos = var.groupby("position")[["surface", "y"]].mean().reset_index()

    for row, (d, grain) in enumerate([(var, "variant"), (pos, "position")]):
        ax = axes[row, col]
        x, y = d.surface.to_numpy(float), d.y.to_numpy(float)

        ax.axvline(syn_surface, color="black", lw=0.5, ls="--", zorder=1)
        ax.axhline(spec["href"], color="black", lw=0.5, ls="--", zorder=1)

        if grain == "variant":
            # Thousands of points: colour by local density, dense drawn on top.
            z = gaussian_kde(np.vstack([x, y]))(np.vstack([x, y]))
            o = z.argsort()
            ax.scatter(x[o], y[o], c=z[o], s=1, cmap="magma",
                       edgecolors="none", alpha=0.85, zorder=2)
        else:
            ax.scatter(x, y, s=1.5, color="black", alpha=0.8,
                       edgecolors="none", zorder=2)

        # LOWESS trend. The white stroke keeps it readable where it crosses the
        # bright core of the density cloud as well as the sparse black tails.
        sm = lowess(y, x, frac=LOESS_FRAC, it=LOESS_IT, return_sorted=True)
        ln, = ax.plot(sm[:, 0], sm[:, 1], color=LOESS_COL, lw=0.75,
                      solid_capstyle="round", zorder=3)
        ln.set_path_effects([pe.withStroke(linewidth=1.75, foreground="white"),
                             pe.Normal()])

        rho, p_rho = spearmanr(x, y)
        r, _ = pearsonr(x, y)
        stats.append(dict(drug=DRUG, parameter=spec["key"], grain=grain, n=len(d),
                          spearman_rho=rho, spearman_p=p_rho,
                          pearson_r=r, pearson_r2=r ** 2))

        ax.text(0.04, 0.04, f"$\\rho$ = {rho:.2f}\nn = {len(d)}",
                transform=ax.transAxes, size=6, ha="left", va="bottom", zorder=4,
                bbox=dict(boxstyle="square,pad=0.15", facecolor="white",
                          edgecolor="none", alpha=0.8))
        ax.set_title(f"{grain.capitalize()} level", pad=2)
        ax.set_xlabel("Surface expression score")
        ax.set_ylabel(spec["ylab"], labelpad=1)
        ax.tick_params(axis="both", pad=1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

fig.subplots_adjust(left=0.10, bottom=0.075, right=0.985, top=0.965,
                    wspace=0.34, hspace=0.40)
out_pdf = OUT / f"surface_vs_{DRUG.lower()}_ec50_emax.pdf"
fig.savefig(out_pdf)
plt.close(fig)

stats_df = pd.DataFrame(stats)
stats_df.to_csv(OUT / f"surface_vs_{DRUG.lower()}_ec50_emax_stats.csv", index=False)
print(stats_df.to_string(index=False))
print(f"\nSaved: {out_pdf.name}")

#!/usr/bin/env python3
"""
Per-position dose-response grid with the Y-AXIS REVERSED, for Morphine, Fentanyl
and DAMGO -- new companion to per_position_curves_{drug} / refit_3param_per_
position_curves.py (originals left untouched).

Only difference from the originals: the y-axis is flipped (top<->bottom) so the
fitted curves read INCREASING instead of decreasing. The curves themselves are
the unchanged raw Hill-3 fits (same Emin/Emax values, same syn reference, same
95% CI, same colours/layout) -- nothing is rescaled or normalized; the axis is
simply reversed via set_ylim(y_hi, y_lo).

Shows all positions present in any drug's missense data (~399), one panel each.
Curve-type sets and dose ranges are kept identical to each drug's original.

Outputs (one per ligand):
  curve_refitting/per_position_curves_reversed_morphine.{pdf,png}
  curve_refitting/per_position_curves_reversed_fentanyl.{pdf,png}
  curve_refitting/per_position_curves_reversed_damgo.{pdf,png}
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT  = ROOT / "curve_refitting"

plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 4,
    "axes.labelsize": 4, "axes.titlesize": 4,
    "xtick.labelsize": 3, "ytick.labelsize": 3,
    "text.color": "black", "axes.edgecolor": "black",
    "axes.labelcolor": "black", "axes.titlecolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "axes.linewidth": 0.3,
    "xtick.major.width": 0.3, "ytick.major.width": 0.3,
    "xtick.major.size": 1.2, "ytick.major.size": 1.2,
    "lines.linewidth": 0.4, "patch.linewidth": 0.3,
})

def hill3(x, emin, emax, log_ec50):
    return emin + (emax - emin) / (1 + 10**(log_ec50 - x))

DRUG_INFO = {
    "Morphine": dict(csv=OUT / "refit_3param_robust_morphine.csv",
                     dose_lo=-12.0, dose_hi=-4.0,
                     keep={"sigmoid", "flat"}, var_color="#3779b9"),
    "Fentanyl": dict(csv=OUT / "refit_3param_robust_fentanyl.csv",
                     dose_lo=-13.0, dose_hi=-4.0,
                     keep={"sigmoid", "flat"}, var_color="#3779b9"),
    "DAMGO":    dict(csv=OUT / "refit_3param_robust_damgo.csv",
                     dose_lo=-12.5, dose_hi=-4.5,
                     keep={"sigmoid", "flat", "no_baseline", "no_plateau"},
                     var_color="#3779b9"),
}

NCOLS    = 20
N_X      = 100
PANEL_MM = 11

# Common position list = union of every drug's missense positions (~all 399)
all_pos = set()
for info in DRUG_INFO.values():
    d = pd.read_csv(info["csv"])
    all_pos |= set(d[d.type == "missense"].position.dropna().astype(int).unique())
POSITIONS = sorted(all_pos)
print(f"Common position grid: {len(POSITIONS)} positions "
      f"({POSITIONS[0]}..{POSITIONS[-1]})")

def plot_drug(drug):
    info = DRUG_INFO[drug]
    df   = pd.read_csv(info["csv"])

    mis = df[df.fit_ok & (df.type == "missense") & df.new_curve_type.isin(info["keep"])].copy()
    syn = df[df.fit_ok & (df.type == "synonymous") & df.new_curve_type.isin(info["keep"])].copy()

    xs = np.linspace(info["dose_lo"], info["dose_hi"], N_X)

    # Synonymous mean fit + 95% CI (raw Hill-3, exactly as the originals).
    y_syn = hill3(xs, syn.fitted_emin.mean(), syn.fitted_emax.mean(),
                  syn.fitted_ec50_logM.mean())
    syn_preds = np.array([hill3(xs, r.fitted_emin, r.fitted_emax, r.fitted_ec50_logM)
                          for _, r in syn.iterrows()])
    syn_p025 = np.percentile(syn_preds, 2.5,  axis=0)
    syn_p975 = np.percentile(syn_preds, 97.5, axis=0)

    # y-range from evaluated curve values (robust to param-pinned fits), + syn band.
    var_preds = np.array([hill3(xs, r.fitted_emin, r.fitted_emax, r.fitted_ec50_logM)
                          for _, r in mis.iterrows()])
    y_lo = min(float(np.percentile(var_preds, 2)), float(syn_p025.min()))
    y_hi = max(float(np.percentile(var_preds, 98)), float(syn_p975.max()))
    y_lo -= 0.1; y_hi += 0.1

    wt_lookup = (df[df.position.notna()].drop_duplicates("position")
                 .set_index("position").wildtype.to_dict())
    by_pos = {int(p): g for p, g in mis.groupby("position")}

    nrows = int(np.ceil(len(POSITIONS) / NCOLS))
    MM = 1/25.4
    fig, axes = plt.subplots(nrows, NCOLS,
                             figsize=(NCOLS*PANEL_MM*MM, (nrows*PANEL_MM + 8)*MM),
                             squeeze=False, sharex=True, sharey=True)

    for idx, pos in enumerate(POSITIONS):
        ax  = axes[idx // NCOLS, idx % NCOLS]
        wt  = wt_lookup.get(float(pos), "?")
        sub = by_pos.get(pos)
        n_var = 0 if sub is None else len(sub)

        ax.fill_between(xs, syn_p025, syn_p975, color="grey", alpha=0.30, lw=0, zorder=1)
        if sub is not None:
            for _, r in sub.iterrows():
                ax.plot(xs, hill3(xs, r.fitted_emin, r.fitted_emax, r.fitted_ec50_logM),
                        color=info["var_color"], lw=0.25, alpha=0.35, zorder=2)
        ax.plot(xs, y_syn, color="black", lw=0.7, zorder=4)
        ax.axhline(0, color="black", lw=0.2, alpha=0.4, zorder=1)

        ax.set_xlim(info["dose_lo"], info["dose_hi"])
        ax.set_ylim(y_hi, y_lo)                 # <-- REVERSED y-axis (curves go up)
        ax.set_title(f"{wt}{int(pos)}  n={n_var}", fontsize=3.5, pad=1)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        for s in ("left", "bottom"): ax.spines[s].set_linewidth(0.3)
        ax.tick_params(axis="both", which="both", length=0, labelsize=0)

    for j in range(len(POSITIONS), nrows * NCOLS):
        axes[j // NCOLS, j % NCOLS].axis("off")

    fig.suptitle(f"{drug}: per-position Hill-3 fits, y-axis reversed (curves increasing). "
                 f"Black = synonymous mean fit; grey band = syn 95% CI.",
                 fontsize=6, y=0.995)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.97, bottom=0.01,
                        hspace=0.6, wspace=0.2)

    out_pdf = OUT / f"per_position_curves_reversed_{drug.lower()}.pdf"
    fig.savefig(out_pdf, dpi=600)
    fig.savefig(out_pdf.with_suffix(".png"), dpi=600)
    print(f"{drug}: {len(mis)} variants, {len(syn)} syn, y-range[{y_lo:.2f},{y_hi:.2f}] "
          f"reversed -> {out_pdf.name}")
    plt.close()

for drug in ["Morphine", "Fentanyl", "DAMGO"]:
    plot_drug(drug)

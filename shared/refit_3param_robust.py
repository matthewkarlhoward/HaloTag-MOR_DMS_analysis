#!/usr/bin/env python3
"""
Robust 3-parameter Hill refit with outlier rejection and principled
classification.

Pipeline per variant:
  1. Initial Hill-3 fit on all available dose points (≥4 valid).
  2. Compute per-point residuals; drop any point with |residual| > OUTLIER_K × σ
     (where σ is the within-fit RMS residual). At most OUTLIER_MAX_DROP points
     dropped; must leave ≥4 for the refit.
  3. Refit on the cleaned data.
  4. Classify the cleaned fit:

       boundary    — EC50 within EC50_EDGE_LOG of the assay dose range, OR any
                     parameter pinned at fit bounds.
       no_plateau  — top-2 dose scores differ by > GAP_THRESHOLD
                     (right-side asymptote not reached at the highest tested dose).
       no_baseline — bottom-2 dose scores differ by > GAP_THRESHOLD
                     (left-side asymptote not reached at the lowest tested dose;
                     EC50 is below the assayed range).
       super_emax  — fitted Emax is deeper than the syn distribution minimum
                     (anchored to the empirical synonymous reference).
       flat        — dynamic range |Emax − Emin| < RANGE_FLAT.
       sigmoid     — everything else (clean fit, both asymptotes within range).

Outputs:
  curve_refitting/refit_3param_robust_morphine.csv
  curve_refitting/refit_3param_robust_fentanyl.csv
  curve_refitting/refit_3param_robust_summary.{pdf,png}
"""
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT  = ROOT / "curve_refitting"

plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 6,
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
})

# ── Hyperparameters ─────────────────────────────────────────────────────────
OUTLIER_K          = 2.0     # drop points with |residual| > K * RMS_residual
OUTLIER_MAX_DROP   = 2       # at most 2 points dropped per variant
MIN_POINTS         = 4       # need at least this many points for a fit
EC50_EDGE_LOG      = 0.5     # EC50 within 0.5 log of dose range = boundary
GAP_THRESHOLD      = 0.15    # |top1-top2| or |bot1-bot2| > this = no asymptote
RANGE_FLAT         = 0.20    # |Emax-Emin| < this = flat
PARAM_BOUND_PAD    = 0.1     # |param| within this of fit bound = pinned

DOSE_RANGE = {
    "Morphine": (-11.5, -4.5),   # log[M]
    "Fentanyl": (-12.0, -5.0),
}

CLASS_COLORS = {
    "sigmoid":     "#2c7c2c",
    "flat":        "#e08214",
    "no_plateau":  "#7a2e7a",
    "no_baseline": "#225a9c",
    "super_emax":  "#1f6f99",
    "boundary":    "#c0273a",
    "fit_fail":    "#888888",
}

# ── Hill-3 ───────────────────────────────────────────────────────────────────
def hill3(x, emin, emax, log_ec50):
    return emin + (emax - emin) / (1 + 10**(log_ec50 - x))

def fit_once(x, y, bounds):
    p0 = [
        float(np.mean(y[:max(1, len(y)//3)])),
        float(np.mean(y[-max(1, len(y)//3):])),
        float(np.median(x)),
    ]
    popt, _ = curve_fit(hill3, x, y, p0=p0, bounds=bounds, maxfev=5000)
    return popt

def robust_fit(y_obs, doses, ec50_bounds_pad=2.0):
    y_obs = np.asarray(y_obs, dtype=float)
    valid = ~np.isnan(y_obs)
    if valid.sum() < MIN_POINTS:
        return None
    x = doses[valid].astype(float)
    y = y_obs[valid].astype(float)
    bounds = ([-2.0, -2.0, x.min() - ec50_bounds_pad],
              [ 2.0,  2.0, x.max() + ec50_bounds_pad])

    # Initial fit
    try:
        popt = fit_once(x, y, bounds)
    except (RuntimeError, ValueError):
        return None
    ypred = hill3(x, *popt)
    resid = y - ypred
    rms   = np.sqrt(np.mean(resid ** 2)) if len(resid) else 0.0

    # Identify outliers
    n_dropped = 0
    if rms > 0:
        # threshold relative to RMS, capped to ≥0.1 to avoid flagging tight fits
        sigma = max(rms, 0.05)
        bad   = np.abs(resid) > OUTLIER_K * sigma
        # Limit how many we drop
        if bad.sum() > 0 and bad.sum() <= OUTLIER_MAX_DROP and (valid.sum() - bad.sum()) >= MIN_POINTS:
            n_dropped = int(bad.sum())
            x = x[~bad]; y = y[~bad]
            # Refit on cleaned data
            try:
                popt = fit_once(x, y, bounds)
            except (RuntimeError, ValueError):
                pass

    # Final stats
    ypred = hill3(x, *popt)
    ss_res = float(np.sum((y - ypred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    emin, emax, ec50 = popt

    # Plateau / baseline checks on the CLEANED data, sorted by dose
    order = np.argsort(x)
    x_sorted = x[order]; y_sorted = y[order]
    if len(y_sorted) >= 2:
        plateau_gap  = float(abs(y_sorted[-1] - y_sorted[-2]))  # right asymptote
        baseline_gap = float(abs(y_sorted[1]  - y_sorted[0]))   # left asymptote
    else:
        plateau_gap = baseline_gap = np.nan

    # Boundary checks
    hit_param_bound = (abs(popt[0] - bounds[0][0]) < PARAM_BOUND_PAD
                       or abs(popt[0] - bounds[1][0]) < PARAM_BOUND_PAD
                       or abs(popt[1] - bounds[0][1]) < PARAM_BOUND_PAD
                       or abs(popt[1] - bounds[1][1]) < PARAM_BOUND_PAD
                       or abs(popt[2] - bounds[0][2]) < PARAM_BOUND_PAD
                       or abs(popt[2] - bounds[1][2]) < PARAM_BOUND_PAD)

    return dict(fitted_emin=emin, fitted_emax=emax,
                fitted_ec50_logM=ec50, r2=r2,
                n_points_used=int(len(x)),
                n_outliers_dropped=n_dropped,
                plateau_gap=plateau_gap,
                baseline_gap=baseline_gap,
                hit_param_bound=hit_param_bound,
                fit_ok=True)

# ── Per-drug dose columns ────────────────────────────────────────────────────
dms = pd.read_csv(ROOT / "dms_scores/composite_dms_scores.csv")

def get_dose_cols(prefix, divisor=1.0):
    pattern = re.compile(prefix + r"(-\d+)M_rescaled_score$")
    pairs = []
    for c in dms.columns:
        m = pattern.search(c)
        if m:
            pairs.append((c, int(m.group(1)) / divisor))
    pairs.sort(key=lambda x: x[1])
    return pairs

DRUGS = {
    "Morphine": get_dose_cols("MOR_Morphine_", divisor=10.0),
    "Fentanyl": get_dose_cols("MOR_Fentanyl_", divisor=1.0),
}

# ── Refit + classify ─────────────────────────────────────────────────────────
def classify(df, drug):
    lo, hi = DOSE_RANGE[drug]
    syn = df[(df.fit_ok) & (df.type == "synonymous")]
    syn_emax_min = float(syn.fitted_emax.min()) if len(syn) else -1.0

    out = df.copy()
    out["dynamic_range"] = (out.fitted_emax - out.fitted_emin).abs()

    fit_failed = ~out.fit_ok
    boundary   = (out.fit_ok &
                  ((out.fitted_ec50_logM < lo + EC50_EDGE_LOG) |
                   (out.fitted_ec50_logM > hi - EC50_EDGE_LOG) |
                   out.hit_param_bound))
    flat        = out.fit_ok & ~boundary & (out.dynamic_range < RANGE_FLAT)
    no_plateau  = (out.fit_ok & ~boundary & ~flat
                   & (out.plateau_gap > GAP_THRESHOLD))
    no_baseline = (out.fit_ok & ~boundary & ~flat & ~no_plateau
                   & (out.baseline_gap > GAP_THRESHOLD))
    super_emax  = (out.fit_ok & ~boundary & ~flat & ~no_plateau & ~no_baseline
                   & (out.fitted_emax < syn_emax_min - 0.05))
    sigmoid     = (out.fit_ok & ~boundary & ~flat
                   & ~no_plateau & ~no_baseline & ~super_emax)

    out["new_curve_type"] = "fit_fail"
    out.loc[boundary,    "new_curve_type"] = "boundary"
    out.loc[flat,        "new_curve_type"] = "flat"
    out.loc[no_plateau,  "new_curve_type"] = "no_plateau"
    out.loc[no_baseline, "new_curve_type"] = "no_baseline"
    out.loc[super_emax,  "new_curve_type"] = "super_emax"
    out.loc[sigmoid,     "new_curve_type"] = "sigmoid"
    out["syn_emax_min_reference"] = syn_emax_min
    return out

def refit_drug(drug, pairs):
    cols  = [c for c, _ in pairs]
    doses = np.array([d for _, d in pairs], dtype=float)
    rows = []
    for _, r in dms.iterrows():
        y = r[cols].values
        fit = robust_fit(y, doses)
        if fit is None:
            fit = dict(fitted_emin=np.nan, fitted_emax=np.nan,
                       fitted_ec50_logM=np.nan, r2=np.nan,
                       n_points_used=0, n_outliers_dropped=0,
                       plateau_gap=np.nan, baseline_gap=np.nan,
                       hit_param_bound=False, fit_ok=False)
        row = dict(hgvs=r.hgvs, type=r.type, position=r.position,
                   wildtype=r.wildtype, mutation=r.mutation,
                   original_curve_type=r.get(f"{drug}_curve_type", "n/a"))
        row.update(fit)
        rows.append(row)
    df = pd.DataFrame(rows)
    df = classify(df, drug)
    return df

for drug in ["Morphine","Fentanyl"]:
    print(f"\nRefitting {drug} with outlier rejection...")
    df = refit_drug(drug, DRUGS[drug])
    out_csv = OUT / f"refit_3param_robust_{drug.lower()}.csv"
    df.to_csv(out_csv, index=False)
    print(f"  Saved: {out_csv}")
    # Summary
    print(f"  Variants total: {len(df)}, fit_ok: {df.fit_ok.sum()}")
    print(f"  Variants with outliers dropped: {(df.n_outliers_dropped > 0).sum()}")
    print(f"  syn_emax_min reference: {df.syn_emax_min_reference.iloc[0]:.3f}")
    print("  Classification:")
    for ct, n in df.new_curve_type.value_counts().items():
        print(f"    {ct:>12}: {n}")
    print("  Classification (missense only):")
    sub = df[df.type == "missense"]
    for ct, n in sub.new_curve_type.value_counts().items():
        print(f"    {ct:>12}: {n}")

# ── Summary plot ─────────────────────────────────────────────────────────────
MM = 1/25.4
fig, axes = plt.subplots(2, 4, figsize=(220*MM, 110*MM))

ORDER = ["sigmoid","flat","no_plateau","no_baseline","super_emax","boundary"]
for col_i, drug in enumerate(["Morphine", "Fentanyl"]):
    df = pd.read_csv(OUT / f"refit_3param_robust_{drug.lower()}.csv")
    sub = df[(df.fit_ok) & (df.type == "missense")]

    # Panel 1: Emax distribution by class
    ax = axes[col_i, 0]
    bins = np.linspace(-2, 1, 60)
    for ct in ORDER:
        m = sub[sub.new_curve_type == ct]
        if not len(m): continue
        ax.hist(m.fitted_emax, bins=bins, color=CLASS_COLORS[ct], alpha=0.65,
                edgecolor="none", label=f"{ct} (n={len(m)})")
    ax.set_xlabel(r"fitted $E_\mathrm{max}$")
    ax.set_ylabel("variants")
    ax.set_title(f"{drug} (missense)", color="black", pad=2)
    ax.set_xlim(-2.1, 1.1)
    ax.legend(frameon=False, loc="upper right", fontsize=4,
              handlelength=0.6, handletextpad=0.3)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    # Panel 2: Plateau gap distribution by class
    ax = axes[col_i, 1]
    bins = np.linspace(0, 1, 50)
    for ct in ORDER:
        m = sub[sub.new_curve_type == ct]
        if not len(m): continue
        ax.hist(m.plateau_gap.dropna(), bins=bins, color=CLASS_COLORS[ct],
                alpha=0.65, edgecolor="none", label=f"{ct} (n={len(m)})")
    ax.axvline(GAP_THRESHOLD, color="black", lw=0.4, ls="--", alpha=0.5)
    ax.set_xlabel("plateau gap (|top1 − top2|)")
    ax.set_ylabel("variants")
    ax.set_title(f"{drug} plateau", color="black", pad=2)
    ax.legend(frameon=False, loc="upper right", fontsize=4,
              handlelength=0.6, handletextpad=0.3)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    # Panel 3: Baseline gap distribution by class
    ax = axes[col_i, 2]
    for ct in ORDER:
        m = sub[sub.new_curve_type == ct]
        if not len(m): continue
        ax.hist(m.baseline_gap.dropna(), bins=bins, color=CLASS_COLORS[ct],
                alpha=0.65, edgecolor="none", label=f"{ct} (n={len(m)})")
    ax.axvline(GAP_THRESHOLD, color="black", lw=0.4, ls="--", alpha=0.5)
    ax.set_xlabel("baseline gap (|bot1 − bot2|)")
    ax.set_ylabel("variants")
    ax.set_title(f"{drug} baseline", color="black", pad=2)
    ax.legend(frameon=False, loc="upper right", fontsize=4,
              handlelength=0.6, handletextpad=0.3)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    # Panel 4: R² distribution by class
    ax = axes[col_i, 3]
    for ct in ORDER:
        m = sub[sub.new_curve_type == ct]
        if not len(m): continue
        ax.hist(m.r2.dropna(), bins=bins, color=CLASS_COLORS[ct],
                alpha=0.65, edgecolor="none", label=f"{ct} (n={len(m)})")
    ax.set_xlabel(r"refit R$^2$")
    ax.set_ylabel("variants")
    ax.set_title(f"{drug} fit quality", color="black", pad=2)
    ax.set_xlim(0, 1.05)
    ax.legend(frameon=False, loc="upper left", fontsize=4,
              handlelength=0.6, handletextpad=0.3)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

fig.tight_layout()
out_pdf = OUT / "refit_3param_robust_summary.pdf"
out_png = OUT / "refit_3param_robust_summary.png"
fig.savefig(out_pdf, dpi=600)
fig.savefig(out_png, dpi=600)
print(f"\nSaved: {out_pdf}")
print(f"Saved: {out_png}")

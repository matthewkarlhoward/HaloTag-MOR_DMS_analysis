#!/usr/bin/env python3
"""
Robust 3-parameter Hill refit for DAMGO dose-response data.

Mirrors refit_3param_robust.py (Morphine / Fentanyl) adapted for DAMGO's
6-point dose schedule. The 0 µM condition is treated as log[M] = -12 (per
project convention), serving as a low-dose baseline anchor 3 logs below the
next-lowest titration point.

Dose mapping (composite column → log[M]):
    DAMGO_mor_dms_1 →  -5   (10 µM)
    DAMGO_mor_dms_2 →  -6   (1 µM)
    DAMGO_mor_dms_3 →  -7   (0.1 µM)
    DAMGO_mor_dms_4 →  -8   (0.01 µM)
    DAMGO_mor_dms_5 →  -9   (0.001 µM)
    DAMGO_mor_dms_6 → -12   (0 µM, anchor)

Differences from the Mor/Fent pipeline
--------------------------------------
* No `boundary` class. The 3-log gap between the -12 anchor and -9 makes the
  EC50-edge rule used for Mor/Fent uninformative here. Variants whose EC50
  would have been flagged as boundary stay in their otherwise-assigned class
  (typically sigmoid). The `hit_param_bound` column is still written to the
  output so a downstream filter can drop pinned fits if needed.
* Six dose points instead of eight (one is the anchor); MIN_POINTS=4 still
  applies, so any variant with ≤3 measured doses fails out.

Pipeline per variant
--------------------
  1. Initial Hill-3 fit on all available dose points (≥4 valid).
  2. Compute per-point residuals; drop any point with |residual| > OUTLIER_K × σ
     (σ = within-fit RMS residual, floored at 0.05). At most OUTLIER_MAX_DROP
     points dropped; must leave ≥4 for the refit.
  3. Refit on the cleaned data.
  4. Classify the cleaned fit:

       no_plateau  — top-2 dose scores differ by > GAP_THRESHOLD
                     (right-side asymptote not reached at the highest tested dose).
       no_baseline — bottom-2 dose scores differ by > GAP_THRESHOLD
                     (left-side asymptote not reached; the -12 anchor doesn't
                     join cleanly with the next-lowest measured dose).
       super_emax  — fitted Emax deeper than the syn distribution minimum
                     (anchored to the empirical synonymous reference).
       flat        — dynamic range |Emax − Emin| < RANGE_FLAT.
       sigmoid     — everything else (clean fit, both asymptotes within range).

Outputs:
  curve_refitting/refit_3param_robust_damgo.csv
  curve_refitting/refit_3param_robust_damgo_summary.{pdf,png}
"""
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

# ── Hyperparameters (identical to robust.py for Mor/Fent) ────────────────────
OUTLIER_K          = 2.0     # drop points with |residual| > K * RMS_residual
OUTLIER_MAX_DROP   = 2       # at most 2 points dropped per variant
MIN_POINTS         = 4       # need at least this many points for a fit
GAP_THRESHOLD      = 0.15    # |top1-top2| or |bot1-bot2| > this = no asymptote
RANGE_FLAT         = 0.20    # |Emax-Emin| < this = flat
PARAM_BOUND_PAD    = 0.1     # |param| within this of fit bound = pinned (recorded only)

# Dose mapping: DAMGO_mor_dms_N → log[M]
DAMGO_DOSE_MAP = {
    1:  -5.0,   # 10 µM
    2:  -6.0,   # 1 µM
    3:  -7.0,   # 0.1 µM
    4:  -8.0,   # 0.01 µM
    5:  -9.0,   # 0.001 µM
    6: -12.0,   # 0 µM anchor
}

CLASS_COLORS = {
    "sigmoid":     "#2c7c2c",
    "flat":        "#e08214",
    "no_plateau":  "#7a2e7a",
    "no_baseline": "#225a9c",
    "super_emax":  "#1f6f99",
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

    try:
        popt = fit_once(x, y, bounds)
    except (RuntimeError, ValueError):
        return None
    ypred = hill3(x, *popt)
    resid = y - ypred
    rms   = np.sqrt(np.mean(resid ** 2)) if len(resid) else 0.0

    n_dropped = 0
    if rms > 0:
        sigma = max(rms, 0.05)
        bad   = np.abs(resid) > OUTLIER_K * sigma
        if bad.sum() > 0 and bad.sum() <= OUTLIER_MAX_DROP and (valid.sum() - bad.sum()) >= MIN_POINTS:
            n_dropped = int(bad.sum())
            x = x[~bad]; y = y[~bad]
            try:
                popt = fit_once(x, y, bounds)
            except (RuntimeError, ValueError):
                pass

    ypred = hill3(x, *popt)
    ss_res = float(np.sum((y - ypred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    emin, emax, ec50 = popt

    order = np.argsort(x)
    x_sorted = x[order]; y_sorted = y[order]
    if len(y_sorted) >= 2:
        plateau_gap  = float(abs(y_sorted[-1] - y_sorted[-2]))   # right asymptote
        baseline_gap = float(abs(y_sorted[1]  - y_sorted[0]))    # left asymptote
    else:
        plateau_gap = baseline_gap = np.nan

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

# ── Load composite + select DAMGO dose columns ───────────────────────────────
dms = pd.read_csv(ROOT / "dms_scores/composite_dms_scores.csv")

def get_damgo_dose_cols():
    pairs = []
    for n, dose in DAMGO_DOSE_MAP.items():
        col = f"DAMGO_mor_dms_{n}_rescaled_score"
        if col in dms.columns:
            pairs.append((col, dose))
        else:
            print(f"  WARNING: column {col} not found in composite_dms_scores.csv")
    pairs.sort(key=lambda x: x[1])
    return pairs

DAMGO_PAIRS = get_damgo_dose_cols()
print(f"DAMGO: {len(DAMGO_PAIRS)} doses at log[M] = "
      + ", ".join(f"{d:g}" for _, d in DAMGO_PAIRS))

# ── Refit + classify ─────────────────────────────────────────────────────────
def classify(df):
    """5-class classifier (no `boundary` — DAMGO dose density is too low for
    the EC50-edge rule to be informative). fit_fail covers variants that
    didn't converge or had <4 valid points."""
    syn = df[(df.fit_ok) & (df.type == "synonymous")]
    syn_emax_min = float(syn.fitted_emax.min()) if len(syn) else -1.0

    out = df.copy()
    out["dynamic_range"] = (out.fitted_emax - out.fitted_emin).abs()

    flat        = out.fit_ok & (out.dynamic_range < RANGE_FLAT)
    no_plateau  = (out.fit_ok & ~flat
                   & (out.plateau_gap > GAP_THRESHOLD))
    no_baseline = (out.fit_ok & ~flat & ~no_plateau
                   & (out.baseline_gap > GAP_THRESHOLD))
    super_emax  = (out.fit_ok & ~flat & ~no_plateau & ~no_baseline
                   & (out.fitted_emax < syn_emax_min - 0.05))
    sigmoid     = (out.fit_ok & ~flat & ~no_plateau
                   & ~no_baseline & ~super_emax)

    out["new_curve_type"] = "fit_fail"
    out.loc[flat,        "new_curve_type"] = "flat"
    out.loc[no_plateau,  "new_curve_type"] = "no_plateau"
    out.loc[no_baseline, "new_curve_type"] = "no_baseline"
    out.loc[super_emax,  "new_curve_type"] = "super_emax"
    out.loc[sigmoid,     "new_curve_type"] = "sigmoid"
    out["syn_emax_min_reference"] = syn_emax_min
    return out

def refit_damgo():
    cols  = [c for c, _ in DAMGO_PAIRS]
    doses = np.array([d for _, d in DAMGO_PAIRS], dtype=float)
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
                   original_curve_type=r.get("DAMGO_curve_type", "n/a"))
        row.update(fit)
        rows.append(row)
    df = pd.DataFrame(rows)
    df = classify(df)
    return df

print("\nRefitting DAMGO with outlier rejection...")
df = refit_damgo()
out_csv = OUT / "refit_3param_robust_damgo.csv"
df.to_csv(out_csv, index=False)
print(f"  Saved: {out_csv}")
print(f"  Variants total:                  {len(df)}")
print(f"  fit_ok:                          {df.fit_ok.sum()}")
print(f"  Variants with outliers dropped:  {(df.n_outliers_dropped > 0).sum()}")
print(f"  Variants with param pinned:      {df.hit_param_bound.sum()}")
print(f"  syn_emax_min reference:          {df.syn_emax_min_reference.iloc[0]:.3f}")
print("  Classification (all variants):")
for ct, n in df.new_curve_type.value_counts().items():
    print(f"    {ct:>12}: {n}")
print("  Classification (missense only):")
sub_mis = df[df.type == "missense"]
for ct, n in sub_mis.new_curve_type.value_counts().items():
    print(f"    {ct:>12}: {n}")

# ── Summary plot ─────────────────────────────────────────────────────────────
MM = 1/25.4
fig, axes = plt.subplots(1, 4, figsize=(220*MM, 60*MM))

ORDER = ["sigmoid","flat","no_plateau","no_baseline","super_emax"]
sub = df[(df.fit_ok) & (df.type == "missense")]

# Panel 1: Emax distribution by class
ax = axes[0]
bins = np.linspace(-2, 1, 60)
for ct in ORDER:
    m = sub[sub.new_curve_type == ct]
    if not len(m): continue
    ax.hist(m.fitted_emax, bins=bins, color=CLASS_COLORS[ct], alpha=0.65,
            edgecolor="none", label=f"{ct} (n={len(m)})")
ax.set_xlabel(r"fitted $E_\mathrm{max}$")
ax.set_ylabel("variants")
ax.set_title("DAMGO (missense)", color="black", pad=2)
ax.set_xlim(-2.1, 1.1)
ax.legend(frameon=False, loc="upper right", fontsize=4,
          handlelength=0.6, handletextpad=0.3)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

# Panel 2: Plateau gap distribution by class
ax = axes[1]
bins_gap = np.linspace(0, 1, 50)
for ct in ORDER:
    m = sub[sub.new_curve_type == ct]
    if not len(m): continue
    ax.hist(m.plateau_gap.dropna(), bins=bins_gap, color=CLASS_COLORS[ct],
            alpha=0.65, edgecolor="none", label=f"{ct} (n={len(m)})")
ax.axvline(GAP_THRESHOLD, color="black", lw=0.4, ls="--", alpha=0.5)
ax.set_xlabel("plateau gap (|top1 − top2|)")
ax.set_ylabel("variants")
ax.set_title("DAMGO plateau", color="black", pad=2)
ax.legend(frameon=False, loc="upper right", fontsize=4,
          handlelength=0.6, handletextpad=0.3)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

# Panel 3: Baseline gap distribution by class
ax = axes[2]
for ct in ORDER:
    m = sub[sub.new_curve_type == ct]
    if not len(m): continue
    ax.hist(m.baseline_gap.dropna(), bins=bins_gap, color=CLASS_COLORS[ct],
            alpha=0.65, edgecolor="none", label=f"{ct} (n={len(m)})")
ax.axvline(GAP_THRESHOLD, color="black", lw=0.4, ls="--", alpha=0.5)
ax.set_xlabel("baseline gap (|bot1 − bot2|)")
ax.set_ylabel("variants")
ax.set_title("DAMGO baseline", color="black", pad=2)
ax.legend(frameon=False, loc="upper right", fontsize=4,
          handlelength=0.6, handletextpad=0.3)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

# Panel 4: R² distribution by class
ax = axes[3]
bins_r2 = np.linspace(0, 1, 50)
for ct in ORDER:
    m = sub[sub.new_curve_type == ct]
    if not len(m): continue
    ax.hist(m.r2.dropna(), bins=bins_r2, color=CLASS_COLORS[ct],
            alpha=0.65, edgecolor="none", label=f"{ct} (n={len(m)})")
ax.set_xlabel(r"refit R$^2$")
ax.set_ylabel("variants")
ax.set_title("DAMGO fit quality", color="black", pad=2)
ax.set_xlim(0, 1.05)
ax.legend(frameon=False, loc="upper left", fontsize=4,
          handlelength=0.6, handletextpad=0.3)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

fig.tight_layout()
out_pdf = OUT / "refit_3param_robust_damgo_summary.pdf"
out_png = OUT / "refit_3param_robust_damgo_summary.png"
fig.savefig(out_pdf, dpi=600)
fig.savefig(out_png, dpi=600)
print(f"\nSaved: {out_pdf}")
print(f"Saved: {out_png}")

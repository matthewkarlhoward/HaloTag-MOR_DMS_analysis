#!/usr/bin/env python3
"""
Normalise the K100 / V175 / R278 raw BRET curves the way the WT / A119L /
double-mutant curves were normalised.

The rule was reverse-engineered from the 20260427 plate, whose raw BRET is in
figures/singles_raw_points.csv (source "doubles_raw") and whose normalised form
is figures/doubles1_points.csv:

    for EACH curve (ligand x variant) separately,
        fit  Y = Bottom + (Top - Bottom) / (1 + 10^(LogEC50 - X))   (Hill = 1)
        to that curve's own raw replicates, then divide every raw point by the
        fitted Top.

So each curve's own no-drug plateau becomes 1.0 and the agonist response reads
downward from it (signal-down assay). Running --validate reproduces all 32
normalised 20260427 curves from the raw plate to 5 decimals (31/32 exactly; the
flat Nalbuphine_100N non-responder differs by 2.6%, where the plateau is not
identifiable and Prism's optimiser landed elsewhere).

This is the per-curve scaling used for the dose-response panels. It is NOT the
WT-DAMGO anchoring used for the efficacy bars/heatmaps, which is applied later
on top of these normalised spans.

Inputs : figures/singles_raw_points.csv   (extract_singles_prism.py)
Outputs: figures/singles_norm_points.csv  long points + norm column
         figures/singles_norm_params.csv  per-curve Top / Bottom / LogEC50 / fit
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import f as _fdist

REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "figures" / "singles_raw_points.csv"
OUT_POINTS = REPO / "figures" / "singles_norm_points.csv"
OUT_PARAMS = REPO / "figures" / "singles_norm_params.csv"

KEY = ["source", "run", "ligand", "variant"]


def sigmoid(x, bottom, top, logec50):
    return bottom + (top - bottom) / (1.0 + 10.0 ** (logec50 - x))


def fit_curve(x, y):
    """3-parameter fit + F-test against a flat line.

    Returns (baseline, plateau, logec50, r2, pval): `baseline` is the plateau at
    the DILUTE end (no drug, Prism's Bottom) and `plateau` the one at the
    saturating end (Prism's Top). This is a signal-down assay, so baseline is
    the larger of the two for a real agonist. NaNs if the fit will not converge.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    lo, hi = float(np.min(x)), float(np.max(x))
    try:
        p, _ = curve_fit(sigmoid, x, y,
                         p0=[float(np.max(y)), float(np.min(y)), -8.0],
                         bounds=([0.0, 0.0, lo], [np.inf, np.inf, hi]),
                         maxfev=40000)
    except Exception:
        return (np.nan,) * 5
    ss_sig = float(np.sum((y - sigmoid(x, *p)) ** 2))
    ss_flat = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_sig / ss_flat if ss_flat > 0 else np.nan
    df1, df2 = 2, len(y) - 3
    if ss_sig > 0 and df2 > 0 and ss_flat > ss_sig:
        F = ((ss_flat - ss_sig) / df1) / (ss_sig / df2)
        pval = float(1.0 - _fdist.cdf(F, df1, df2))
    else:
        pval = 1.0
    return float(p[0]), float(p[1]), float(p[2]), r2, pval


def normalise(df, reps=None, dose_range=None):
    """Add a `norm` column = bret / (that curve's fitted Top).

    reps        : keep only these replicate indices (None = all)
    dose_range  : (lo, hi) logM window used for the fit; the vehicle well
                  (logM == 0) is always excluded. None = every real dose.
    """
    d = df[df.logM != 0].copy()
    if reps is not None:
        d = d[d.rep.isin(reps)]
    if dose_range is not None:
        lo, hi = dose_range
        d = d[(d.logM >= lo) & (d.logM <= hi)]

    params = []
    for key, g in d.groupby(KEY, sort=False):
        base, plateau, lec50, r2, pval = fit_curve(g.logM, g.bret)
        params.append(dict(zip(KEY, key)) | dict(
            n=len(g), n_rep=g.rep.nunique(), n_dose=g.logM.nunique(),
            bret_baseline=base, bret_plateau=plateau, logec50=lec50,
            r2=r2, pval=pval, divisor=base, responsive=bool(pval < 0.05),
            pinned=bool(np.isfinite(lec50) and
                        min(abs(lec50 - g.logM.min()),
                            abs(lec50 - g.logM.max())) < 1e-6),
            emax_norm=plateau / base if base and np.isfinite(base) else np.nan,
            window=1 - plateau / base if base and np.isfinite(base) else np.nan))
    par = pd.DataFrame(params)

    out = d.merge(par[KEY + ["divisor"]], on=KEY, how="left")
    out["norm"] = out.bret / out.divisor
    return out, par


def validate():
    """Reproduce figures/doubles1_points.csv from the raw 20260427 plate."""
    raw = pd.read_csv(RAW)
    raw = raw[raw.source == "doubles_raw"]
    ref = pd.read_csv(REPO / "figures" / "doubles1_points.csv")
    # doubles1 used the first four replicates and doses -13 .. -4
    got, par = normalise(raw, reps=[0, 1, 2, 3], dose_range=(-13, -4))
    # doubles1 spells a curve "<ligand>_<variant>"; here they are two columns
    got = got.assign(variant=got.ligand + "_" + got.variant)
    m = ref.merge(got, on=["ligand", "variant", "logM", "rep"],
                  suffixes=("_ref", ""))
    if len(m) != len(ref):
        print(f"!! matched {len(m)} of {len(ref)} reference rows")
    err = (m.norm - m.response).abs() / m.response
    per = m.assign(err=err).groupby("variant").err.max().sort_values()
    print(f"max relative error overall : {err.max():.2e}")
    print(f"curves within 1e-4         : {(per < 1e-4).sum()} / {len(per)}")
    print("worst curves:")
    print(per.tail(3).to_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true",
                    help="check the rule against figures/doubles1_points.csv")
    args = ap.parse_args()
    if args.validate:
        validate()
        return

    raw = pd.read_csv(RAW)
    out, par = normalise(raw)
    out.to_csv(OUT_POINTS, index=False)
    par.to_csv(OUT_PARAMS, index=False)
    print(f"-> {OUT_POINTS.relative_to(REPO)}  ({len(out)} rows)")
    print(f"-> {OUT_PARAMS.relative_to(REPO)}  ({len(par)} curves)")
    show = par[par.variant.str.contains("100|175|278|WT")]
    cols = ["source", "ligand", "variant", "bret_baseline", "emax_norm",
            "window", "logec50", "r2", "pval", "responsive", "pinned"]
    print(show[cols].round(4).to_string(index=False))


if __name__ == "__main__":
    main()

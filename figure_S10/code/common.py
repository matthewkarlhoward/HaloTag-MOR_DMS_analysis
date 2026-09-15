"""Shared paths, zone definitions and helpers for the modeling project."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

PROJ = Path(__file__).resolve().parents[1]
REPO = PROJ.parent
PROC = PROJ / "data" / "processed"          # -> mor_efficacy/data/processed
NORM = PROJ / "norm"
ESM = PROJ / "esm"
DIFF = PROJ / "difference"
EVAL = PROJ / "eval"
CONFIGS = PROJ / "configs"
for _d in (NORM, ESM, DIFF, EVAL, CONFIGS):
    _d.mkdir(parents=True, exist_ok=True)

# --- structured receptor only ---------------------------------------------
# TM1 through helix 8. The N-terminus (1-64) and C-terminal tail (356-400) are
# disordered in every mu-OR structure, carry no BW number, and have no reliable
# structural covariates. They are filtered BEFORE analysis, not after: they
# change the permutation/FDR calibration, the conditioning quantiles, the
# expression LOWESS fit and the zone denominators. 65-355 is the same window
# used for the ligand-profile clustering elsewhere in this repo.
POS_MIN, POS_MAX = 65, 355


def core(df: "pd.DataFrame", col: str = "position") -> "pd.DataFrame":
    return df[df[col].between(POS_MIN, POS_MAX)].copy()


def core_mask(positions) -> "np.ndarray":
    p = np.asarray(positions)
    return (p >= POS_MIN) & (p <= POS_MAX)


# Position-level noise ceiling on the efficacy axis, measured by replicate
# split-half in the companion project (mor_efficacy/RESULTS.md). Every efficacy
# correlation reported here is divided by sqrt of this.
CEILING_POS = 0.7606
DISATT = float(np.sqrt(CEILING_POS))

# --- zones (handoff section 4: stratify EVERY correlation by zone) ---------
# gate = the conserved activation microswitches. pocket = orthosteric shell.
# periphery = transducer-facing, defined by distance to the Gi alpha5 footprint.
GATE_ELEMENTS = ["DRY", "NPxxY", "PIF", "CWxP", "Na_site", "TM6_lock"]
PERIPHERY_A = 8.0
ZONE_ORDER = ["pocket", "gate", "periphery", "other"]


def assign_zone(row) -> str:
    if row.element in GATE_ELEMENTS:
        return "gate"
    if row.element == "orthosteric":
        return "pocket"
    if pd.notna(row.gprot_dist_min) and row.gprot_dist_min <= PERIPHERY_A:
        return "periphery"
    return "other"


def load_positions() -> pd.DataFrame:
    p = pd.read_parquet(PROC / "positions.parquet")
    p["zone"] = p.apply(assign_zone, axis=1)
    return p


def disattenuate(rho: float, axis: str) -> float:
    """Efficacy correlations get the measured ceiling; others are returned as-is."""
    return rho / DISATT if axis in ("efficacy", "efficacy_resid") else rho


def spearman_ci(x, y, n_boot=1000, seed=0):
    from scipy.stats import spearmanr
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 8:
        return np.nan, np.nan, np.nan, len(x)
    rho = spearmanr(x, y).statistic
    rng = np.random.default_rng(seed)
    b = [spearmanr(x[i], y[i]).statistic
         for i in (rng.integers(0, len(x), len(x)) for _ in range(n_boot))]
    lo, hi = np.nanpercentile(b, [2.5, 97.5])
    return rho, lo, hi, len(x)


def perm_p(x, y, n_perm=1000, seed=0):
    from scipy.stats import spearmanr
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 8:
        return np.nan
    obs = abs(spearmanr(x, y).statistic)
    rng = np.random.default_rng(seed)
    null = [abs(spearmanr(x, rng.permutation(y)).statistic) for _ in range(n_perm)]
    return float(np.mean(np.asarray(null) >= obs))

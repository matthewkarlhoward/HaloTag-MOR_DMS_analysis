#!/usr/bin/env python3
"""
Shared computation for the EC50/Emax overlap analysis (morphine).

Answers two different questions that the "separable" claim conflates:

  (1) Do EC50-altering and Emax-altering positions form DISTINCT SETS?
      No. At the paper's own LoF definition (>1 SD of the synonymous
      distribution) the two sets overlap ~2x more than chance, and 64% of
      EC50-LoF positions are also Emax-LoF.

  (2) Is there a POTENCY-SPECIFIC component once the shared efficacy axis is
      removed?  Yes, and it is localised to the orthosteric pocket.

Position effects use the same filters as figure 3B/C (make_chimerax_param_maps.py):
  EC50  -> sigmoid + no_baseline   (a flat curve has no identifiable EC50)
  Emax  -> sigmoid + flat + no_baseline  (a dead variant is a real Emax effect)
  n >= 5 missense variants per position.

Emax is on the Activity scale: synonymous mean Emin -> 0, mean Emax -> 1.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[2]
REFIT = ROOT / "curve_refitting/refit_3param_robust_morphine.csv"
DIST  = ROOT / "structures/processed/ligand_distances/per_pdb_legacy/8ef6_distances.csv"

EC50_KEEP = {"sigmoid", "no_baseline"}
EMAX_KEEP = {"sigmoid", "flat", "no_baseline"}
MIN_N     = 5

MOTIF = {166: "DRY", 167: "DRY", 168: "DRY",
         334: "NPxxY", 335: "NPxxY", 336: "NPxxY", 337: "NPxxY", 338: "NPxxY",
         157: "PIF", 246: "PIF", 291: "PIF",
         116: "Na+", 149: "Na+",
         293: "CWxP", 294: "CWxP", 295: "CWxP", 296: "CWxP"}

COL = {"EC50": "#1f4fd8", "Emax": "#d42b2b", "both": "#7b3fa0", "neither": "#d9d9d9"}


def _split_half(df, col, seed=0, nrep=200):
    """Spearman-Brown corrected split-half reliability of the position means."""
    rng = np.random.default_rng(seed)
    groups = [v[col].to_numpy() for _, v in df.groupby("position") if len(v) >= 6]
    out = []
    for _ in range(nrep):
        a, b = [], []
        for v in groups:
            idx = rng.permutation(len(v)); h = len(v) // 2
            a.append(v[idx[:h]].mean()); b.append(v[idx[h:2 * h]].mean())
        r = pearsonr(a, b)[0]
        out.append(2 * r / (1 + r))
    return float(np.median(out))


def load():
    """-> (table, meta).  One row per position with both parameters measured."""
    d = pd.read_csv(REFIT)
    syn = d[(d.type == "synonymous") & (d.new_curve_type.isin(EC50_KEEP))]
    e0, e1 = syn.fitted_emin.mean(), syn.fitted_emax.mean()
    ec_syn = syn.fitted_ec50_logM.mean()

    d = d.copy()
    d["act"] = (d.fitted_emax - e0) / (e1 - e0)

    # 1 SD of the synonymous distribution -> the paper's LoF threshold
    syn_all = d[d.type == "synonymous"]
    sd_ec = syn_all[syn_all.new_curve_type.isin(EC50_KEEP)].fitted_ec50_logM.std()
    sd_em = syn_all[syn_all.new_curve_type.isin(EMAX_KEEP)].act.std()

    mis = d[d.type == "missense"]
    ec = mis[mis.new_curve_type.isin(EC50_KEEP)].dropna(subset=["fitted_ec50_logM"])
    em = mis[mis.new_curve_type.isin(EMAX_KEEP)].dropna(subset=["act"])

    E = ec.groupby("position").fitted_ec50_logM.agg(["mean", "count"])
    E = E[E["count"] >= MIN_N]["mean"]
    M = em.groupby("position").act.agg(["mean", "count"])
    M = M[M["count"] >= MIN_N]["mean"]

    pos = sorted(set(E.index) & set(M.index))
    t = pd.DataFrame({"ec50": E.loc[pos], "act": M.loc[pos]})
    t["dec50"] = t.ec50 - ec_syn          # rightward shift, log units
    t["eloss"] = 1 - t.act                # efficacy loss, activity units
    t["wt"] = d.groupby("position").wildtype.first().reindex(pos)
    t["motif"] = pd.Series(MOTIF).reindex(pos).fillna("")

    dist = pd.read_csv(DIST).set_index("residue_number_human")
    t["dist"] = dist.distance_sidechain_angstrom.reindex(pos)

    # ---- classes at the paper's LoF definition -----------------------------
    t["lof_ec50"] = t.dec50 > sd_ec
    t["lof_emax"] = t.eloss > sd_em
    t["cls"] = np.where(t.lof_ec50 & t.lof_emax, "both",
               np.where(t.lof_ec50, "EC50",
               np.where(t.lof_emax, "Emax", "neither")))

    # ---- potency-specific residual (noise-corrected Deming regression) -----
    rel_ec = _split_half(ec, "fitted_ec50_logM")
    rel_em = _split_half(em, "act")
    lam = (t.dec50.var() * (1 - rel_ec)) / (t.eloss.var() * (1 - rel_em))
    sxx, syy = t.eloss.var(), t.dec50.var()
    sxy = np.cov(t.eloss, t.dec50)[0, 1]
    slope = ((syy - lam * sxx) + np.sqrt((syy - lam * sxx) ** 2 + 4 * lam * sxy ** 2)) / (2 * sxy)
    icept = t.dec50.mean() - slope * t.eloss.mean()
    t["resid"] = t.dec50 - (icept + slope * t.eloss)

    meta = dict(ec_syn=ec_syn, sd_ec=sd_ec, sd_em=sd_em,
                rel_ec=rel_ec, rel_em=rel_em, lam=lam,
                slope=slope, icept=icept,
                r_obs=pearsonr(t.dec50, t.eloss)[0],
                r_dis=pearsonr(t.dec50, t.eloss)[0] / np.sqrt(rel_ec * rel_em),
                n=len(t))
    return t, meta


def label(t, p):
    r = t.loc[p]
    return f"{r.wt}{int(p)}"


if __name__ == "__main__":
    t, m = load()
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in m.items()})
    print(t.cls.value_counts().to_string())


def bootstrap_stability(t, meta, nrep=300, seed=0, cache=True):
    """How often does a position keep its class when the variants at that
    position are resampled?  Returns a frame of class frequencies."""
    out = Path(__file__).resolve().parent / "class_stability.csv"
    if cache and out.exists():
        return pd.read_csv(out, index_col=0)
    d = pd.read_csv(REFIT)
    syn = d[(d.type == "synonymous") & (d.new_curve_type.isin(EC50_KEEP))]
    e0, e1 = syn.fitted_emin.mean(), syn.fitted_emax.mean()
    ec_syn = syn.fitted_ec50_logM.mean()
    d = d.copy(); d["act"] = (d.fitted_emax - e0) / (e1 - e0)
    mis = d[d.type == "missense"]
    ec = {p: v.fitted_ec50_logM.dropna().to_numpy()
          for p, v in mis[mis.new_curve_type.isin(EC50_KEEP)].groupby("position")}
    em = {p: v.act.dropna().to_numpy()
          for p, v in mis[mis.new_curve_type.isin(EMAX_KEEP)].groupby("position")}
    pos = [p for p in t.index if len(ec.get(p, [])) >= MIN_N and len(em.get(p, [])) >= MIN_N]
    rng = np.random.default_rng(seed)
    keep = pd.DataFrame(0.0, index=pos, columns=["EC50", "Emax", "both", "neither"])
    for _ in range(nrep):
        de = np.array([rng.choice(ec[p], len(ec[p]), replace=True).mean() for p in pos]) - ec_syn
        dm = 1 - np.array([rng.choice(em[p], len(em[p]), replace=True).mean() for p in pos])
        a, b = de > meta["sd_ec"], dm > meta["sd_em"]
        cl = np.where(a & b, "both", np.where(a, "EC50", np.where(b, "Emax", "neither")))
        for p, c in zip(pos, cl):
            keep.loc[p, c] += 1
    keep /= nrep
    if cache:
        keep.to_csv(out)
    return keep

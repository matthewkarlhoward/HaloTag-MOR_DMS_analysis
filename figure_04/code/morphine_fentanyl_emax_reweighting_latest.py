#!/usr/bin/env python3
"""
Morphine vs Fentanyl Emax reweighting — recomputed on the LATEST robust 3-param
refits (curve_refitting/refit_3param_robust_{drug}.csv), replacing the stale
composite_dms_scores.csv (April 7) used by the original
morphine_fentanyl_emax_lof_comparison.py / _reweighting_bootstrap.py.

Reports, so the "shared network, different weights" claim can be stated with
noise control:

  POSITION LEVEL (per-position mean ΔEmax, sigmoid, n>=5 both ligands)
    - overall Spearman ρ (all shared positions)          -> paper "overall" number
    - Spearman ρ within LoF union / intersection         -> paper "deleterious subset"
    - within-ligand bootstrap NOISE CEILING vs cross-ligand ρ
        (the drop is meaningful only if cross-ligand ρ < within-ligand ρ)
    - # positions significantly reweighted beyond noise (BH q<0.05, LoF union)

  VARIANT LEVEL (matches the paper's current "variants" wording)
    - overall Spearman ρ (all anchored missense)
    - Spearman ρ in the >1 SD LoF-union subset

Convention: fitted_emax scale has syn-sigmoid mean ~ -0.5; ΔEmax = emax - syn_mean,
POSITIVE ΔEmax = loss of function.
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CR   = ROOT / "curve_refitting"
MIN_VARIANTS = 5
LOF_K        = 1.0
N_BOOT       = 10000
RNG          = np.random.default_rng(20260727)

def load(drug):
    return pd.read_csv(CR / f"refit_3param_robust_{drug.lower()}.csv")

mor_raw = load("morphine")
fen_raw = load("fentanyl")

# ───────────────────────── POSITION-LEVEL (lof_comparison logic) ──────────────
def per_drug_positions(df):
    syn = df[(df.type == "synonymous") & (df.new_curve_type == "sigmoid")]
    syn_mean, syn_sd = syn.fitted_emax.mean(), syn.fitted_emax.std(ddof=1)
    mis = df[(df.type == "missense") & (df.new_curve_type == "sigmoid")]
    pos = (mis.groupby("position")
              .agg(emax=("fitted_emax", "mean"),
                   n=("fitted_emax", "count"),
                   wt=("wildtype", "first"))
              .dropna())
    pos = pos[pos.n >= MIN_VARIANTS]
    pos["delta_emax"] = pos.emax - syn_mean
    return pos.reset_index(), syn_mean, syn_sd

mor_pos, mor_syn, mor_sd = per_drug_positions(mor_raw)
fen_pos, fen_syn, fen_sd = per_drug_positions(fen_raw)
th_mor, th_fen = LOF_K * mor_sd, LOF_K * fen_sd

shared = (mor_pos[["position", "wt", "delta_emax"]]
          .rename(columns={"delta_emax": "d_mor"})
          .merge(fen_pos[["position", "delta_emax"]]
                 .rename(columns={"delta_emax": "d_fen"}), on="position"))
shared["lof_mor"] = shared.d_mor > th_mor
shared["lof_fen"] = shared.d_fen > th_fen

rho_all, p_all = spearmanr(shared.d_mor, shared.d_fen)
uni = shared[shared.lof_mor | shared.lof_fen]
inter = shared[shared.lof_mor & shared.lof_fen]
rho_uni, p_uni = spearmanr(uni.d_mor, uni.d_fen)
rho_int, p_int = spearmanr(inter.d_mor, inter.d_fen)

print("="*70)
print("POSITION LEVEL  (latest robust refits)")
print("="*70)
print(f"syn emax SD: Mor {mor_sd:.4f} (thr {th_mor:.4f}) | Fen {fen_sd:.4f} (thr {th_fen:.4f})")
print(f"shared positions (sigmoid, n>=5 both): {len(shared)}")
print(f"  LoF union: {len(uni)}   LoF intersection: {len(inter)}")
print(f"  overall  Spearman ρ = {rho_all:.3f}  (p={p_all:.1e}, n={len(shared)})")
print(f"  LoF-union ρ         = {rho_uni:.3f}  (p={p_uni:.1e}, n={len(uni)})")
print(f"  LoF-intersect ρ     = {rho_int:.3f}  (p={p_int:.1e}, n={len(inter)})")

# ───────────────────────── BOOTSTRAP NOISE CEILING ────────────────────────────
# Resample missense variants within each shared position; two independent draws
# per ligand give the within-ligand noise ceiling (max ρ if weights identical).
positions = shared.position.tolist()
def var_index(df):
    m = df[(df.type == "missense") & (df.new_curve_type == "sigmoid")]
    return {p: m.loc[m.position == p, "fitted_emax"].values for p in positions}
mor_v, fen_v = var_index(mor_raw), var_index(fen_raw)

def boot_matrix(vmap):
    """(N_BOOT x n_pos) bootstrap means, vectorized per position."""
    out = np.empty((N_BOOT, len(positions)))
    for i, p in enumerate(positions):
        a = vmap[p]
        idx = RNG.integers(0, len(a), size=(N_BOOT, len(a)))
        out[:, i] = a[idx].mean(axis=1)
    return out

mor_b1, mor_b2 = boot_matrix(mor_v), boot_matrix(mor_v)
fen_b1, fen_b2 = boot_matrix(fen_v), boot_matrix(fen_v)

def rho_rows(A, B):
    return np.array([spearmanr(A[i], B[i]).statistic for i in range(A.shape[0])])

rho_within_mor = rho_rows(mor_b1, mor_b2)
rho_within_fen = rho_rows(fen_b1, fen_b2)
rho_cross      = rho_rows(mor_b1, fen_b1)

def ci(x): return np.quantile(x, 0.025), np.quantile(x, 0.975)
print("\n" + "="*70)
print(f"NOISE CEILING vs CROSS-LIGAND ρ  ({N_BOOT} bootstraps, all shared positions)")
print("="*70)
print(f"  within-Mor ρ : mean {rho_within_mor.mean():.3f}  95%CI [{ci(rho_within_mor)[0]:.3f}, {ci(rho_within_mor)[1]:.3f}]")
print(f"  within-Fen ρ : mean {rho_within_fen.mean():.3f}  95%CI [{ci(rho_within_fen)[0]:.3f}, {ci(rho_within_fen)[1]:.3f}]")
print(f"  cross    ρ   : mean {rho_cross.mean():.3f}  95%CI [{ci(rho_cross)[0]:.3f}, {ci(rho_cross)[1]:.3f}]")
gap = min(rho_within_mor.mean(), rho_within_fen.mean()) - rho_cross.mean()
print(f"  => cross-ligand ρ sits {gap:.3f} BELOW the lower noise ceiling")
print(f"     ({(rho_cross < np.quantile(rho_within_mor,0.05)).mean()*100:.1f}% of cross iters below within-Mor 5th pctile)")

# ───────────────────────── PER-POSITION REWEIGHTING TEST ──────────────────────
# %syn-lost units; diff = Mor - Fen; bootstrap p that 0 is excluded; BH on LoF union.
mor_pct = 100 * (1 - mor_b1 / mor_syn)
fen_pct = 100 * (1 - fen_b1 / fen_syn)
diff_boot = mor_pct - fen_pct
mor_obs = np.array([100*(1 - mor_v[p].mean()/mor_syn) for p in positions])
fen_obs = np.array([100*(1 - fen_v[p].mean()/fen_syn) for p in positions])
diff_obs = mor_obs - fen_obs

def p_two(sam, obs):
    if obs > 0: return max(2*np.mean(sam <= 0), 1/N_BOOT)
    if obs < 0: return max(2*np.mean(sam >= 0), 1/N_BOOT)
    return 1.0
pvals = np.array([p_two(diff_boot[:, i], diff_obs[i]) for i in range(len(positions))])

mor_pct_sd = 100 * mor_sd / abs(mor_syn)
fen_pct_sd = 100 * fen_sd / abs(fen_syn)
lof_union_mask = (mor_obs > LOF_K*mor_pct_sd) | (fen_obs > LOF_K*fen_pct_sd)

def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p)
    q = p[o]*n/(np.arange(n)+1); q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = q; return np.clip(out, 0, 1)
q = np.ones(len(positions)); q[lof_union_mask] = bh(pvals[lof_union_mask])

res = pd.DataFrame({"position": positions, "wt": shared.wt.values,
                    "d_mor%": mor_obs, "d_fen%": fen_obs,
                    "diff_MorMinusFen%": diff_obs, "p": pvals, "q": q,
                    "lof_union": lof_union_mask})
res["sig"] = res.q < 0.05
n_sig = res.sig.sum()
print("\n" + "="*70)
print("PER-POSITION REWEIGHTING (BH q<0.05 within LoF union)")
print("="*70)
print(f"  significantly reweighted positions: {n_sig}  of {lof_union_mask.sum()} LoF-union tested")
_sig = res[res.sig].copy()
top = _sig.reindex(_sig["diff_MorMinusFen%"].abs().sort_values(ascending=False).index)
print("  top reweighted (|Mor-Fen %syn|):")
for _, r in top.head(12).iterrows():
    who = "MOR more LoF" if r["diff_MorMinusFen%"] > 0 else "FEN more LoF"
    print(f"    {r.wt}{int(r.position):<4} diff={r['diff_MorMinusFen%']:+6.1f}%  q={r.q:.3g}  ({who})")

# ───────────────────────── VARIANT LEVEL (paper's wording) ────────────────────
EMAX_OK = {"sigmoid", "flat", "no_baseline"}
def act_miss(df, pre, syn_mean, syn_e0):
    sub = df[df.fit_ok & (df.type=="missense") & df.new_curve_type.isin(EMAX_OK)].copy()
    sub[f"{pre}_emax"] = (sub.fitted_emax - syn_e0) / (syn_mean - syn_e0)  # 1=WT,0=dead
    sub[f"{pre}_delta"] = sub.fitted_emax - syn_mean  # +=LoF
    return sub[["hgvs", f"{pre}_emax", f"{pre}_delta"]]

mor_e0 = mor_raw[(mor_raw.type=="synonymous")&(mor_raw.new_curve_type=="sigmoid")].fitted_emin.mean()
fen_e0 = fen_raw[(fen_raw.type=="synonymous")&(fen_raw.new_curve_type=="sigmoid")].fitted_emin.mean()
v = act_miss(mor_raw,"Mor",mor_syn,mor_e0).merge(act_miss(fen_raw,"Fen",fen_syn,fen_e0), on="hgvs")
rho_v_all, p_v_all = spearmanr(v.Mor_emax, v.Fen_emax)
sub = v[(v.Mor_delta > th_mor) | (v.Fen_delta > th_fen)]
rho_v_sub, p_v_sub = spearmanr(sub.Mor_emax, sub.Fen_emax)
print("\n" + "="*70)
print("VARIANT LEVEL (anchored missense, Activity-scaled Emax)")
print("="*70)
print(f"  overall Spearman ρ = {rho_v_all:.3f}  (p={p_v_all:.1e}, n={len(v)})")
print(f"  >1 SD LoF-union ρ  = {rho_v_sub:.3f}  (p={p_v_sub:.1e}, n={len(sub)})")

res.to_csv(ROOT/"dms_scores/mor_fent_reweighting_per_position_latest.csv", index=False)
print(f"\nSaved per-position table: dms_scores/mor_fent_reweighting_per_position_latest.csv")

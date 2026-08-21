#!/usr/bin/env python3
"""
Morphine vs Fentanyl EC50 reweighting on the LATEST robust refits — the potency
analogue of morphine_fentanyl_emax_reweighting_latest.py.

EC50 is only defined for sigmoid fits, so we use sigmoid, n>=5 missense/position.
ΔEC50 = per-position mean logM EC50 − synonymous mean logM EC50.
Positive ΔEC50 = right-shift = potency loss (LOF), matching the figure's
"EC50 LOF bias" convention.

Reports: overall/subset Spearman ρ, within-ligand noise ceiling vs cross-ligand
ρ, and per-position positions significantly reweighted (BH q<0.05, LoF union),
so the text can name specific EC50 differentiators.
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CR   = ROOT / "curve_refitting"
COL  = "fitted_ec50_logM"
MIN_VARIANTS = 5
LOF_K        = 1.0
N_BOOT       = 10000
RNG          = np.random.default_rng(20260727)

def load(d): return pd.read_csv(CR / f"refit_3param_robust_{d}.csv")
mor_raw, fen_raw = load("morphine"), load("fentanyl")

def per_drug_positions(df):
    syn = df[(df.type=="synonymous") & (df.new_curve_type=="sigmoid")]
    syn_mean, syn_sd = syn[COL].mean(), syn[COL].std(ddof=1)
    mis = df[(df.type=="missense") & (df.new_curve_type=="sigmoid")].dropna(subset=[COL])
    pos = (mis.groupby("position")
              .agg(v=(COL,"mean"), n=(COL,"count"), wt=("wildtype","first")).dropna())
    pos = pos[pos.n >= MIN_VARIANTS]
    pos["delta"] = pos.v - syn_mean
    return pos.reset_index(), syn_mean, syn_sd

mor_pos, mor_syn, mor_sd = per_drug_positions(mor_raw)
fen_pos, fen_syn, fen_sd = per_drug_positions(fen_raw)
th_mor, th_fen = LOF_K*mor_sd, LOF_K*fen_sd

shared = (mor_pos[["position","wt","delta"]].rename(columns={"delta":"d_mor"})
          .merge(fen_pos[["position","delta"]].rename(columns={"delta":"d_fen"}), on="position"))
shared["lof_mor"] = shared.d_mor > th_mor
shared["lof_fen"] = shared.d_fen > th_fen
rho_all,_ = spearmanr(shared.d_mor, shared.d_fen)
uni  = shared[shared.lof_mor | shared.lof_fen]
rho_uni,_ = spearmanr(uni.d_mor, uni.d_fen)

print("="*68); print("POSITION-LEVEL EC50 (latest robust refits)"); print("="*68)
print(f"syn EC50 SD (logM): Mor {mor_sd:.3f} (thr {th_mor:.3f}) | Fen {fen_sd:.3f} (thr {th_fen:.3f})")
print(f"shared positions (sigmoid n>=5 both): {len(shared)}   LoF union: {len(uni)}")
print(f"  overall  Spearman ρ = {rho_all:.3f}  (n={len(shared)})")
print(f"  LoF-union ρ         = {rho_uni:.3f}  (n={len(uni)})")

# bootstrap noise ceiling
positions = shared.position.tolist()
def vindex(df):
    m = df[(df.type=="missense") & (df.new_curve_type=="sigmoid")].dropna(subset=[COL])
    return {p: m.loc[m.position==p, COL].values for p in positions}
mor_v, fen_v = vindex(mor_raw), vindex(fen_raw)
def bootM(vm):
    out = np.empty((N_BOOT, len(positions)))
    for i,p in enumerate(positions):
        a=vm[p]; idx=RNG.integers(0,len(a),size=(N_BOOT,len(a))); out[:,i]=a[idx].mean(1)
    return out
mor_b1,mor_b2,fen_b1,fen_b2 = bootM(mor_v),bootM(mor_v),bootM(fen_v),bootM(fen_v)
def rr(A,B): return np.array([spearmanr(A[i],B[i]).statistic for i in range(A.shape[0])])
wM,wF,cx = rr(mor_b1,mor_b2), rr(fen_b1,fen_b2), rr(mor_b1,fen_b1)
ci=lambda x:(np.quantile(x,.025),np.quantile(x,.975))
print("\n"+"="*68); print(f"NOISE CEILING vs CROSS-LIGAND ρ ({N_BOOT} boot)"); print("="*68)
print(f"  within-Mor {wM.mean():.3f} {ci(wM)} | within-Fen {wF.mean():.3f} {ci(wF)} | cross {cx.mean():.3f} {ci(cx)}")

# per-position reweighting (raw logM diff)
diffb = (mor_b1-mor_syn)-(fen_b1-fen_syn)
mor_o = np.array([mor_v[p].mean()-mor_syn for p in positions])
fen_o = np.array([fen_v[p].mean()-fen_syn for p in positions])
diffo = mor_o - fen_o
def p2(s,o):
    if o>0: return max(2*np.mean(s<=0),1/N_BOOT)
    if o<0: return max(2*np.mean(s>=0),1/N_BOOT)
    return 1.0
pv = np.array([p2(diffb[:,i],diffo[i]) for i in range(len(positions))])
lof_union = (mor_o>th_mor)|(fen_o>th_fen)
def bh(p):
    p=np.asarray(p);n=len(p);o=np.argsort(p);q=p[o]*n/(np.arange(n)+1)
    q=np.minimum.accumulate(q[::-1])[::-1];out=np.empty(n);out[o]=q;return np.clip(out,0,1)
q=np.ones(len(positions)); q[lof_union]=bh(pv[lof_union])
res=pd.DataFrame({"position":positions,"wt":shared.wt.values,"d_mor":mor_o,"d_fen":fen_o,
                  "diff_MorMinusFen_logM":diffo,"q":q,"lof_union":lof_union})
res["sig"]=res.q<0.05
res["dir"]=np.where(res.diff_MorMinusFen_logM>0,"MOR more LOF","FEN more LOF")
print("\n"+"="*68); print("PER-POSITION EC50 REWEIGHTING (BH q<0.05, LoF union)"); print("="*68)
print(f"  significant positions: {res.sig.sum()} of {lof_union.sum()} tested")
sig=res[res.sig].reindex(res[res.sig]["diff_MorMinusFen_logM"].abs().sort_values(ascending=False).index)
for _,r in sig.iterrows():
    print(f"    {r.wt}{int(r.position):<4} Δlog={r['diff_MorMinusFen_logM']:+5.2f}  q={r.q:.3g}  ({r.dir})")

# cross-check the figure labels
print("\n  figure-label check:")
for pos,nm in {149:"D149",153:"M153",126:"Q126",328:"Y328",299:"H299"}.items():
    r=res[res.position==pos]
    if len(r)==0: print(f"    {nm:<6} EXCLUDED (no sigmoid n>=5 both)")
    else:
        r=r.iloc[0]; print(f"    {nm:<6} Δlog={r['diff_MorMinusFen_logM']:+5.2f} q={r.q:.3g} sig={r.sig} ({r.dir})")

# variant level (matches panel b scatter)
EMAX_OK={"sigmoid"}
def vm(df,pre):
    s=df[df.fit_ok&(df.type=="missense")&df.new_curve_type.isin(EMAX_OK)].dropna(subset=[COL]).copy()
    return s[["hgvs",COL]].rename(columns={COL:f"{pre}"})
v=vm(mor_raw,"Mor").merge(vm(fen_raw,"Fen"),on="hgvs")
rho_v,_=spearmanr(v.Mor,v.Fen)
print(f"\n  VARIANT-level overall Spearman ρ = {rho_v:.3f}  (n={len(v)})")
res.to_csv(ROOT/"dms_scores/mor_fent_ec50_reweighting_per_position_latest.csv",index=False)
print("\nSaved: dms_scores/mor_fent_ec50_reweighting_per_position_latest.csv")

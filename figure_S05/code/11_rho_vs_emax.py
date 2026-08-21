#!/usr/bin/env python3
"""
Scatter: efficacy (Emax) vs per-receptor coupling (rho, gate-scaled), morphine, colored by
surface expression. Shows that variants with the SAME low Emax split into two mechanisms:
  - rho ~ 0 + low expression  -> expression-limited (trafficking)
  - rho < 0 + normal expression -> genuine coupling loss
This is what rho adds beyond raw Emax.
"""
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; BIN=ROOT/"dms_data/damgo_counts/counts 6"
REPFLUOR=np.array([3.9,4.9,5.22,5.5])
per=[]
for r in ["R1","R2","R3","R4"]:
    cols=[pd.read_csv(BIN/f"MOR_surface_{r}_B{b}.csv")[["hgvs","count"]].pipe(
          lambda x:(x.set_index("hgvs")["count"]/x["count"].sum()).rename(f"b{b}")) for b in range(4)]
    M=pd.concat(cols,axis=1).fillna(0);P=M.div(M.sum(axis=1),axis=0)
    per.append(pd.Series((P.values*REPFLUOR).sum(1),index=M.index,name=r))
mf=pd.concat(per,axis=1).mean(axis=1).rename("mean_log10F").reset_index()
op=pd.read_csv(ROOT/"operational_model/morphine_operational_fits.csv").merge(mf,on="hgvs",how="left")
op["log10_S"]=op.mean_log10F-op[op.type=="synonymous"].mean_log10F.mean()
tw=np.log10(op[(op.type=="synonymous")&op.tau.notna()&~op.saturating.astype(bool)].tau.median())
op["rho"]=(np.log10(op.tau)-tw)-op.log10_S
wtE=op[op.type=="synonymous"].emax_obs_data.median()
op["Emax"]=op.emax_obs_data/wtE                       # efficacy, WT=1
c=op[(op.curve_type=="sigmoid")&(op.rmse<0.12)&~op.saturating.astype(bool)&op.type.eq("missense")].dropna(subset=["rho","Emax","log10_S"])

plt.rcParams.update({"font.family":"Helvetica","font.size":6,"mathtext.default":"regular",
    "axes.labelsize":6,"axes.titlesize":6,"xtick.labelsize":6,"ytick.labelsize":6,
    "text.color":"black","axes.edgecolor":"black","axes.labelcolor":"black","axes.titlecolor":"black",
    "xtick.color":"black","ytick.color":"black","axes.linewidth":0.5,
    "xtick.major.width":0.5,"ytick.major.width":0.5,"lines.linewidth":0.5})
MM=1/25.4
fig,ax=plt.subplots(figsize=(70*MM,60*MM))
ax.axhline(0,color="grey",lw=.5,ls=":"); ax.axvline(1,color="grey",lw=.5,ls=":")
sc=ax.scatter(c.Emax,c.rho,c=c.log10_S,cmap="PuOr",vmin=-0.4,vmax=0.4,s=2,lw=0,alpha=.7)
ax.set_ylim(-1.5,1.5)
ax.set_xlabel("E$_{max}$ (efficacy, WT=1)",labelpad=1)
ax.set_ylabel("log$_{10}$ ρ  (per-receptor coupling)",labelpad=1)
ax.tick_params(pad=1); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
cb=fig.colorbar(sc,ax=ax,fraction=.05,pad=.02); cb.set_label("log$_{10}$ S (expression)",fontsize=5); cb.ax.tick_params(labelsize=4)
# guide text
ax.text(0.30,0.03,"low Emax, ρ≈0\n= expression-limited",fontsize=4.5,color="#7f3b08",va="center")
ax.text(0.28,-0.55,"low Emax, ρ<0\n= coupling loss",fontsize=4.5,color="#2d004b",va="center")
fig.tight_layout()
out=ROOT/"operational_model/plots/fig12_rho_vs_emax"
fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png",dpi=400)
print(f"n={len(c)}  saved {out}")

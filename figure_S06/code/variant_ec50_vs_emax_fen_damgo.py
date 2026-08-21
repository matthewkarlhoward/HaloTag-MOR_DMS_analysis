#!/usr/bin/env python3
"""
Variant-level EC50 (x, log[M]) vs Emax (y, Activity) scatter for FENTANYL and DAMGO,
matching morphine_variant_ec50_vs_emax.py. One point = one sigmoid missense variant;
Emax on the Activity scale (syn Emin->0, Emax->1); density-colored, dense-on-top;
grey dotted lines = synonymous mean EC50 (vertical) and Emax=1 (horizontal).
Syn reference (Emin/Emax/EC50) uses the {sigmoid,no_baseline} synonymous set so the
scale is well anchored for DAMGO (super-agonist; syn are mostly no_baseline).
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from pathlib import Path
HERE=Path(__file__).resolve().parent
REF={"sigmoid","no_baseline"}
# per-drug: curve types to plot, drop pinned artifacts, and the ASSAYED EC50 range
# (points with EC50 outside the tested concentrations are extrapolated -> removed)
CFG={"fentanyl":dict(keep={"sigmoid"},              drop_pinned=False, assay=(-12.0,-5.0)),
     "damgo":   dict(keep={"sigmoid","no_baseline"},drop_pinned=True,  assay=(-12.0,-5.0))}
plt.rcParams.update({"font.family":"Helvetica","font.size":6,"mathtext.default":"regular","axes.labelsize":6,
    "axes.titlesize":6,"xtick.labelsize":6,"ytick.labelsize":6,"legend.fontsize":6,
    "text.color":"black","axes.edgecolor":"black","axes.labelcolor":"black",
    "axes.titlecolor":"black","xtick.color":"black","ytick.color":"black",
    "axes.linewidth":0.5,"xtick.major.width":0.5,"ytick.major.width":0.5,
    "xtick.major.size":2,"ytick.major.size":2,"lines.linewidth":0.5,"patch.linewidth":0.5})
MM=1/25.4

DRUGX={"fentanyl":"fentanyl","damgo":"DAMGO"}
for drug in ["fentanyl","damgo"]:
    cfg=CFG[drug]
    df=pd.read_csv(HERE/f"refit_3param_robust_{drug}.csv"); df=df[df.fit_ok]
    syn=df[(df.type=="synonymous")&(df.new_curve_type.isin(REF))]
    syn_emin=float(syn.fitted_emin.mean()); syn_emax=float(syn.fitted_emax.mean())
    syn_ec50=float(df[(df.type=="synonymous")&(df.new_curve_type=="sigmoid")].fitted_ec50_logM.mean())
    mis=df[(df.type=="missense")&(df.new_curve_type.isin(cfg["keep"]))].copy()
    mis["emax_act"]=(mis.fitted_emax-syn_emin)/(syn_emax-syn_emin)
    n_all=len(mis)
    if cfg["drop_pinned"]: mis=mis[~mis.hit_param_bound.astype(bool)]
    lo,hi=cfg["assay"]
    mis=mis[(mis.fitted_ec50_logM>=lo)&(mis.fitted_ec50_logM<=hi)&(mis.emax_act>=0.0)]
    ec50=mis.fitted_ec50_logM.to_numpy(); emax=mis.emax_act.to_numpy()
    z=gaussian_kde(np.vstack([ec50,emax]))(np.vstack([ec50,emax])); o=z.argsort()
    ec50,emax,z=ec50[o],emax[o],z[o]
    fig,ax=plt.subplots(figsize=(40*MM,40*MM))
    ax.axvline(syn_ec50,color="grey",lw=0.5,ls=":",zorder=1)
    ax.axhline(1.0,color="grey",lw=0.5,ls=":",zorder=1)
    ax.scatter(ec50,emax,c=z,s=1,cmap="magma",edgecolors="none",alpha=0.85,zorder=2)
    ax.set_xlabel(f"EC$_{{50}}$ log[{DRUGX[drug]}]",labelpad=1)
    ax.set_ylabel("E$_{max}$ (Activity)",labelpad=1)
    ax.set_xlim(*cfg["assay"])
    ax.tick_params(axis="both",pad=1)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    for s in ("left","bottom"): ax.spines[s].set_linewidth(0.5)
    fig.subplots_adjust(left=0.205,bottom=0.175,right=0.97,top=0.98)
    out=HERE/f"{drug}_variant_ec50_vs_emax.pdf"
    fig.savefig(out,dpi=600); fig.savefig(out.with_suffix(".png"),dpi=600); plt.close()
    print(f"{DRUGX[drug]}: n={len(ec50)} plotted ({n_all-len(mis)} removed) | curve types {cfg['keep']} "
          f"| syn_ec50={syn_ec50:.2f} -> {out.name}")

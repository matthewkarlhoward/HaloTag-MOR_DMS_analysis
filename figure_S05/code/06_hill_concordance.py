#!/usr/bin/env python3
"""
Supplemental concordance figure: the operational-model fit recovers the same
potency and observed efficacy as the robust three-parameter Hill fit (morphine).
Compares EQUIVALENT/observed parameters (not tau, which is the deconvolved,
reserve-corrected efficacy and is deliberately a nonlinear transform of Emax).

Panels (operational y vs robust-Hill x, with 1:1 line):
  A  potency        : operational logEC50   vs  Hill fitted_ec50_logM
  B  observed Emax  : operational emax_obs   vs  Hill activation (emin-emax)
  C  baseline (Emin): operational B          vs  Hill fitted_emin
"""
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, spearmanr, pearsonr
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
op=pd.read_csv(ROOT/"operational_model/morphine_operational_fits.csv")
rob=pd.read_csv(ROOT/"curve_refitting/refit_3param_robust_morphine.csv")
m=op.merge(rob[["hgvs","fitted_ec50_logM","fitted_emax","fitted_emin","new_curve_type"]],
           on="hgvs",how="inner")
m=m[(m.new_curve_type=="sigmoid")&(m.curve_type=="sigmoid")&(m.rmse<0.08)
    &(~m.saturating.astype(bool))].copy()   # well-fit curves (rmse<0.08)
m["hill_activation"]=m.fitted_emin-m.fitted_emax   # +=active, matches op emax_obs sign

plt.rcParams.update({"font.family":"Helvetica","font.size":6,"axes.labelsize":6,
    "axes.titlesize":6,"xtick.labelsize":6,"ytick.labelsize":6,
    "text.color":"black","axes.edgecolor":"black","axes.labelcolor":"black",
    "axes.titlecolor":"black","xtick.color":"black","ytick.color":"black",
    "axes.linewidth":0.5,"xtick.major.width":0.5,"ytick.major.width":0.5,
    "xtick.major.size":2,"ytick.major.size":2})
MM=1/25.4

def panel(ax,x,y,xlab,ylab,title,lims):
    x=np.asarray(x,float);y=np.asarray(y,float)
    ok=np.isfinite(x)&np.isfinite(y);x=x[ok];y=y[ok]
    ax.plot(lims,lims,color="grey",ls=":",lw=0.5,zorder=1)
    z=gaussian_kde(np.vstack([x,y]))(np.vstack([x,y]));o=z.argsort()
    ax.scatter(x[o],y[o],c=z[o],s=1.5,cmap="magma",lw=0,alpha=0.85,zorder=2)
    sp=spearmanr(x,y).statistic; pe=pearsonr(x,y).statistic
    ax.set_xlim(lims);ax.set_ylim(lims);ax.set_box_aspect(1)
    ax.set_xlabel(xlab,labelpad=1);ax.set_ylabel(ylab,labelpad=1)
    ax.set_title(f"{title}\nρ={sp:.2f}, r={pe:.2f}, n={len(x)}",pad=2)
    ax.tick_params(pad=1)
    for s in ("top","right"): ax.spines[s].set_visible(False)

fig,axes=plt.subplots(1,3,figsize=(150*MM,52*MM))
panel(axes[0],m.fitted_ec50_logM,m.logEC50,
      "Hill logEC$_{50}$ (M)","operational logEC$_{50}$ (M)","Potency",(-11,-4))
panel(axes[1],m.hill_activation,m.emax_obs_fit,
      "Hill E$_{max}$ (activation)","operational E$_{max}$ (observed)","Efficacy (observed)",(0,0.75))
panel(axes[2],m.fitted_emin,m.B,
      "Hill E$_{min}$","operational baseline","Basal",(-0.5,0.5))
fig.tight_layout(w_pad=1.2)
out=ROOT/"operational_model/plots/fig9_hill_concordance"
fig.savefig(f"{out}.pdf");fig.savefig(f"{out}.png",dpi=300)
print("saved",out)
print(f"n(clean sigmoid both) = {len(m)}")
for lab,a,b in [("EC50","fitted_ec50_logM","logEC50"),
                ("Emax(activation)","hill_activation","emax_obs_fit"),
                ("Emin","fitted_emin","B")]:
    x=m[a].values;y=m[b].values;ok=np.isfinite(x)&np.isfinite(y)
    print(f"  {lab:18s} Spearman {spearmanr(x[ok],y[ok]).statistic:+.3f}  Pearson {pearsonr(x[ok],y[ok]).statistic:+.3f}")

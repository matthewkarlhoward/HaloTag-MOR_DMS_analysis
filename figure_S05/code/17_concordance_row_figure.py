#!/usr/bin/env python3
"""
Final single-row figure (a-d), 165 mm wide x <=50 mm tall:
  a  Potency         operational logEC50 vs Hill logEC50        (from 06_hill_concordance)
  b  Efficacy (obs)  operational Emax     vs Hill activation     (from 06_hill_concordance)
  c  Basal           operational baseline vs Hill Emin           (from 06_hill_concordance)
  d  rho vs Emax      per-receptor coupling vs efficacy, by expr (from 11_rho_vs_emax)
Panel letters: bold Helvetica size 8. All other text: regular Helvetica size 6.
Exact width: Agg backend, fixed figsize, constrained_layout, NO bbox_inches='tight'.
"""
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, spearmanr, pearsonr
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
MM = 1/25.4

# ---------- data for a,b,c : operational vs robust-Hill concordance ----------
op  = pd.read_csv(ROOT/"operational_model/morphine_operational_fits.csv")
rob = pd.read_csv(ROOT/"curve_refitting/refit_3param_robust_morphine.csv")
m = op.merge(rob[["hgvs","fitted_ec50_logM","fitted_emax","fitted_emin","new_curve_type"]],
             on="hgvs", how="inner")
m = m[(m.new_curve_type=="sigmoid")&(m.curve_type=="sigmoid")&(m.rmse<0.08)
      &(~m.saturating.astype(bool))].copy()
m["hill_activation"] = m.fitted_emin - m.fitted_emax

# ---------- data for d : rho (gate-scaled) vs Emax, colored by expression ----
BIN = ROOT/"dms_data/damgo_counts/counts 6"; REPFLUOR = np.array([3.9,4.9,5.22,5.5])
per=[]
for r in ["R1","R2","R3","R4"]:
    cols=[pd.read_csv(BIN/f"MOR_surface_{r}_B{b}.csv")[["hgvs","count"]].pipe(
          lambda x:(x.set_index("hgvs")["count"]/x["count"].sum()).rename(f"b{b}")) for b in range(4)]
    M=pd.concat(cols,axis=1).fillna(0); P=M.div(M.sum(axis=1),axis=0)
    per.append(pd.Series((P.values*REPFLUOR).sum(1),index=M.index,name=r))
mf=pd.concat(per,axis=1).mean(axis=1).rename("mean_log10F").reset_index()
d=pd.read_csv(ROOT/"operational_model/morphine_operational_fits.csv").merge(mf,on="hgvs",how="left")
d["log10_S"]=d.mean_log10F-d[d.type=="synonymous"].mean_log10F.mean()
tw=np.log10(d[(d.type=="synonymous")&d.tau.notna()&~d.saturating.astype(bool)].tau.median())
d["rho"]=(np.log10(d.tau)-tw)-d.log10_S
wtE=d[d.type=="synonymous"].emax_obs_data.median(); d["Emax"]=d.emax_obs_data/wtE
d=d[(d.curve_type=="sigmoid")&(d.rmse<0.12)&~d.saturating.astype(bool)&d.type.eq("missense")]\
   .dropna(subset=["rho","Emax","log10_S"])

# ---------- style ----------
plt.rcParams.update({"font.family":"Helvetica","mathtext.default":"regular",
    "font.size":6,"axes.labelsize":6,"axes.titlesize":6,"xtick.labelsize":6,"ytick.labelsize":6,
    "text.color":"black","axes.edgecolor":"black","axes.labelcolor":"black","axes.titlecolor":"black",
    "xtick.color":"black","ytick.color":"black","axes.linewidth":0.5,
    "xtick.major.width":0.5,"ytick.major.width":0.5,"xtick.major.size":2,"ytick.major.size":2})

# exact manual placement (mm -> figure fraction); square panels, colorbar hugging d
W,H = 165.0, 50.0
side = 29.0                      # square panel side (mm)
y0   = 8.0                       # bottom margin (xlabel+ticks)
xL   = 11.0                      # left margin (a's ylabel)
gap  = 9.0                       # inter-panel gap (ylabels of b,c,d)
xs   = [xL + i*(side+gap) for i in range(4)]
fig  = plt.figure(figsize=(W*MM, H*MM))
def rect(x): return [x/W, y0/H, side/W, side/H]
axes = [fig.add_axes(rect(x)) for x in xs]
cax  = fig.add_axes([(xs[3]+side+1.5)/W, y0/H, 1.6/W, side/H])   # right of panel d

def concord(ax,x,y,xlab,ylab,title,lims):
    x=np.asarray(x,float);y=np.asarray(y,float);ok=np.isfinite(x)&np.isfinite(y);x,y=x[ok],y[ok]
    ax.plot(lims,lims,color="grey",ls=":",lw=0.5,zorder=1)
    z=gaussian_kde(np.vstack([x,y]))(np.vstack([x,y]));o=z.argsort()
    ax.scatter(x[o],y[o],c=z[o],s=1.2,cmap="magma",lw=0,alpha=0.85,zorder=2)
    sp=spearmanr(x,y).statistic;pe=pearsonr(x,y).statistic
    ax.set_xlim(lims);ax.set_ylim(lims)
    ax.set_xlabel(xlab,labelpad=1);ax.set_ylabel(ylab,labelpad=1)
    ax.set_title(f"{title}\nρ={sp:.2f}, r={pe:.2f}, n={len(x)}",pad=2)
    ax.tick_params(pad=1)
    for s in ("top","right"): ax.spines[s].set_visible(False)

concord(axes[0],m.fitted_ec50_logM,m.logEC50,
        "Hill logEC$_{50}$ (M)","operational logEC$_{50}$ (M)","Potency",(-11,-4))
concord(axes[1],m.hill_activation,m.emax_obs_fit,
        "Hill E$_{max}$ (activation)","operational E$_{max}$ (observed)","Efficacy (observed)",(0,0.75))
concord(axes[2],m.fitted_emin,m.B,
        "Hill E$_{min}$","operational baseline","Basal",(-0.5,0.5))

# panel d
ax=axes[3]
ax.axhline(0,color="grey",lw=.5,ls=":"); ax.axvline(1,color="grey",lw=.5,ls=":")
sc=ax.scatter(d.Emax,d.rho,c=d.log10_S,cmap="PuOr",vmin=-0.4,vmax=0.4,s=1.2,lw=0,alpha=.75,zorder=2)
ax.set_xlim(0.25,1.5);ax.set_ylim(-1.5,1.5)
ax.set_xlabel("E$_{max}$ (efficacy, WT=1)",labelpad=1)
ax.set_ylabel("log$_{10}$ ρ (per-receptor coupling)",labelpad=1)
ax.tick_params(pad=1)
for s in ("top","right"): ax.spines[s].set_visible(False)
ax.set_title(" ",pad=2)   # keep top aligned with concordance panels (which have 2-line titles)
ax.text(0.30,1.22,"low E$_{max}$, ρ≈0\n= expression-limited",fontsize=4.5,color="#7f3b08",va="center")
ax.text(0.30,-1.24,"low E$_{max}$, ρ<0\n= coupling loss",fontsize=4.5,color="#2d004b",va="center")
cb=fig.colorbar(sc,cax=cax); cb.set_label("log$_{10}$ S (expression)",fontsize=5)
cb.ax.tick_params(labelsize=4,width=0.5,length=2); cb.outline.set_linewidth(0.5)

# ---------- panel letters: bold Helvetica size 8 ----------
for ax,letter in zip(axes,"abcd"):
    ax.annotate(letter, xy=(0,1), xycoords="axes fraction", xytext=(-24,12),
                textcoords="offset points", fontsize=8, fontweight="bold",
                family="Helvetica", va="bottom", ha="left")

out=ROOT/"operational_model/plots/fig_operational_concordance_row"
fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png",dpi=400)
w,h=fig.get_size_inches()
print(f"saved {out}  ({w/MM:.1f} x {h/MM:.1f} mm)  n_abc={len(m)}  n_d={len(d)}")

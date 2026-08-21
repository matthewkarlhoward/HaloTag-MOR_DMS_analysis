#!/usr/bin/env python3
"""
Supplement: three-way Venn of positions that ALTER EC50 (and Emax) for each
titrated agonist (morphine, fentanyl, DAMGO). "Alter" = per-position mean effect
> 1 SD of that ligand's synonymous distribution (right-shifted EC50 / reduced Emax),
over a common set of testable positions.
  EC50  -> sigmoid n>=5 (EC50 needs a defined sigmoid; DAMGO-limited, ~75 common)
  Emax  -> emax-anchored n>=5 (sigmoid/flat/no_baseline; ~399 common)
"""
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn.layout.venn3 import DefaultLayoutAlgorithm
from pathlib import Path
CR=Path(__file__).resolve().parents[1]/"curve_refitting"
LIG=["morphine","fentanyl","damgo"]; LAB=["Morphine","Fentanyl","DAMGO"]

def altering(prop):
    # EC50: sigmoid + no_baseline (no_baseline keeps a usable inflection; recovers
    # DAMGO super-potent variants whose low-dose baseline is off-scale).
    # Emax: sigmoid + flat + no_baseline (activation is defined for all).
    col="fitted_ec50_logM" if prop=="EC50" else "fitted_emax"
    mis_types=["sigmoid","no_baseline"] if prop=="EC50" else ["sigmoid","flat","no_baseline"]
    syn_types=["sigmoid","no_baseline"]
    perlig={}
    for lg in LIG:
        df=pd.read_csv(CR/f"refit_3param_robust_{lg}.csv")
        syn=df[(df.type=="synonymous")&(df.new_curve_type.isin(syn_types))]
        sm,sd=syn[col].mean(),syn[col].std(ddof=1)
        mis=df[(df.type=="missense")&(df.new_curve_type.isin(mis_types))].dropna(subset=[col])
        g=mis.groupby("position").agg(v=(col,"mean"),n=(col,"count")).dropna()
        g=g[g.n>=5]; g["delta"]=g.v-sm
        perlig[lg]=(g,sd)
    common=set(perlig["morphine"][0].index)
    for lg in LIG[1:]: common&=set(perlig[lg][0].index)
    sets=[]
    for lg in LIG:
        g,sd=perlig[lg]
        sets.append({p for p in common if g.loc[p,"delta"]>1.0*sd})
    return sets,len(common)

plt.rcParams.update({
    "font.family":"Helvetica","font.size":6,"mathtext.default":"regular",
    "text.color":"black","axes.edgecolor":"black","axes.labelcolor":"black",
    "axes.titlecolor":"black","xtick.color":"black","ytick.color":"black",
    "axes.linewidth":0.5,"lines.linewidth":0.5,"patch.linewidth":0.5})
MM=1/25.4
DISP={"EC50":r"EC$_{50}$","Emax":r"E$_\mathrm{max}$"}
fig,axes=plt.subplots(1,2,figsize=(60*MM,40*MM))   # 60 mm wide x 40 mm tall
IDS=['100','010','001','110','101','011','111']
EQ=DefaultLayoutAlgorithm(fixed_subset_sizes=(1,1,1,1,1,1,1))   # equal-size circles
for ax,prop in zip(axes,["EC50","Emax"]):
    sets,ntot=altering(prop)
    v=venn3(sets,set_labels=LAB,ax=ax,layout_algorithm=EQ)
    for pid in IDS:                                          # strip region fills
        p=v.get_patch_by_id(pid)
        if p is not None: p.set_alpha(0); p.set_facecolor("none"); p.set_edgecolor("none")
    venn3_circles((1,1,1,1,1,1,1),ax=ax,layout_algorithm=EQ,color="black",linewidth=0.5)
    for t in (v.set_labels or []):
        if t: t.set_fontsize(6); t.set_color("black")
    for t in (v.subset_labels or []):
        if t: t.set_fontsize(6); t.set_color("black")
    ax.set_title(f"Positions altering {DISP[prop]}", fontsize=6, color="black")
    print(f"{prop}: {ntot} common testable | "+", ".join(f"{LAB[i]} {len(sets[i])}" for i in range(3)))
fig.tight_layout(pad=0.6)
out=CR/"venn_altering_positions"
fig.savefig(f"{out}.pdf");fig.savefig(f"{out}.png",dpi=600)   # fixed figure size (165 mm wide)
print("saved",out)

#!/usr/bin/env python3
"""
ChimeraX maps of per-position EC50 / Emax effect for each titrated agonist, each on
its OWN structure, in the user's captured views (side + top-down), with fentanyl/DAMGO
receptors superposed onto 8EF6 so the same cameras reproduce identical viewpoints.
  Morphine -> 8EF6 chain R (MOI)   Fentanyl -> 8EF5 chain M (7V7)   DAMGO -> 8EFQ chain R (peptide P)
Views: Emax -> side ; EC50 -> side + top-down pocket.
Coloring: white(low) -> Emax red / EC50 blue (saturated at strong LoF); grey=untestable.
Outputs -> curve_refitting/chimerax_maps/
"""
import numpy as np, pandas as pd, gemmi
from pathlib import Path
ROOT=Path("/Users/mkh/GitHub/mor_dms_analysis")
CR=ROOT/"curve_refitting"; OUT=CR/"chimerax_maps"; OUT.mkdir(exist_ok=True)
CIFD=ROOT/"structures/raw/experimental"
LIG={"morphine":"Morphine","fentanyl":"Fentanyl","damgo":"DAMGO"}
REF={"sigmoid","no_baseline"}; EC50_KEEP={"sigmoid","no_baseline"}; EMAX_KEEP={"sigmoid","flat","no_baseline"}
EC50_MAX=1.5; EMAX_MAX=1.0
PAL={"ec50":"white:white:blue","emax":"white:white:red"}   # white low half -> saturated blue/red

# --- user-captured 8EF6 views (chain R) -------------------------------------
CAM_SIDE="-0.66813,0.056198,-0.74192,38.605,0.20656,0.97196,-0.11239,102.95,0.7148,-0.22835,-0.661,35.991"
M6_SIDE=np.array([[-0.85867,-0.19853,0.47252,221.61],[0.20054,0.71828,0.66622,-76.677],[-0.47166,0.66683,-0.57695,181.08]])
CAM_TOP="-0.011909,0.9134,-0.40688,87.959,-0.25662,0.39051,0.88411,208,0.96644,0.11494,0.22975,156.05"
# top view was captured on DAMGO (8EFQ chain R); reference M6_TOP is backed out below
M_DAMGO_TOP=np.array([[-0.86546,-0.31701,-0.38793,366.93],
                      [-0.39492,-0.04475,0.91762,0.74857],
                      [-0.30826,0.94737,-0.08646,59.39]])

STRUCT={
 "morphine":dict(cif="8ef6.cif", rchain="R", lig=":MOI", ligname="MOI"),
 "fentanyl":dict(cif="8ef5.cif", rchain="M", lig=":7V7", ligname="7V7"),
 "damgo":   dict(cif="8efq.cif", rchain="R", lig="/P",   ligname=None),
}

def ca_dict(cif,chain):
    m=gemmi.read_structure(str(CIFD/cif))[0]
    return {r.seqid.num:np.array([a.pos.x,a.pos.y,a.pos.z]) for r in m[chain] for a in r if a.name=="CA"}
REFCA=ca_dict("8ef6.cif","R")

def align_A(cif,chain):
    if cif=="8ef6.cif" and chain=="R": return np.eye(3),np.zeros(3)
    S=ca_dict(cif,chain); common=sorted(set(S)&set(REFCA))
    P=np.array([S[i] for i in common]); Q=np.array([REFCA[i] for i in common])
    pc,qc=P.mean(0),Q.mean(0); H=(P-pc).T@(Q-qc); U,_,Vt=np.linalg.svd(H)
    d=np.sign(np.linalg.det(Vt.T@U.T)); Ra=Vt.T@np.diag([1,1,d])@U.T; ta=qc-Ra@pc
    return Ra,ta
def compose(M6,Ra,ta):
    Rb,tb=M6[:,:3],M6[:,3]; Rm,tm=Rb@Ra, Rb@ta+tb
    return ",".join(f"{x:.5f}" for i in range(3) for x in [Rm[i,0],Rm[i,1],Rm[i,2],tm[i]])

# back out reference (8EF6) top matrix from the DAMGO-captured one: M6 = Rd·Ra^T, td - M6·ta
_Rad,_tad=align_A("8efq.cif","R"); _Rd,_td=M_DAMGO_TOP[:,:3],M_DAMGO_TOP[:,3]
_R6t=_Rd@_Rad.T; M6_TOP=np.hstack([_R6t,(_td-_R6t@_tad)[:,None]])

def per_pos(lc):
    df=pd.read_csv(CR/f"refit_3param_robust_{lc}.csv")
    syn=df[(df.type=="synonymous")&(df.new_curve_type.isin(REF))]
    ec_syn=syn.fitted_ec50_logM.mean(); e0,e1=syn.fitted_emin.mean(),syn.fitted_emax.mean()
    mis=df[df.type=="missense"].copy()
    ec=mis[mis.new_curve_type.isin(EC50_KEEP)].dropna(subset=["fitted_ec50_logM"])
    ecp=ec.groupby("position").fitted_ec50_logM.agg(["mean","count"]); ecp=ecp[ecp["count"]>=5]; ecp=ecp["mean"]-ec_syn
    em=mis[mis.new_curve_type.isin(EMAX_KEEP)].copy(); em["loss"]=1-(em.fitted_emax-e0)/(e1-e0)
    emp=em.groupby("position").loss.agg(["mean","count"]); emp=emp[emp["count"]>=5]; emp=emp["mean"]
    return ecp,emp

def write_defattr(path,attr,series,chain):
    with open(path,"w") as f:
        f.write(f"attribute: {attr}\nrecipient: residues\n")
        for pos,val in series.items(): f.write(f"\t/{chain}:{int(pos)}\t{val:.4f}\n")

def strip_cmd(s):
    rc=s["rchain"]
    return (f"delete ~/{rc} & ~/P\ndelete solvent" if s["ligname"] is None
            else f"delete ~/{rc}\ndelete :CLR\ndelete solvent")

def cxc(path,s,defattr,attr,pal,vmax,title,cam,mmat,view):
    rc=s["rchain"]; peptide=(s["ligname"] is None)
    ligstyle="stick" if (view=="top" or peptide) else "sphere"
    path.write_text(f"""# {title} ({view} view)
close session
open {CIFD/s['cif']}
{strip_cmd(s)}
hide atoms
show /{rc} cartoons
color /{rc} #cccccc
show {s['lig']} atoms
style {s['lig']} {ligstyle}
color {s['lig']} yellow
color {s['lig']} byhetero
open {defattr}
color byattribute r:{attr} /{rc} palette {pal} range 0,{vmax} target c noValueColor #cccccc
set bgColor white
lighting soft
graphics silhouettes true
view matrix camera {cam}
view matrix models #1,{mmat}
2dlabel create t text "{title} ({view})" xpos .03 ypos .95 size 12 color black
""")

DATA={}   # lc -> dict(ec,em,side_mat,top_mat,rc,cif,lig,ligname)
for lc,Lig in LIG.items():
    s=STRUCT[lc]; ec,em=per_pos(lc); rc=s["rchain"]
    Ra,ta=align_A(s["cif"],rc); side_mat=compose(M6_SIDE,Ra,ta); top_mat=compose(M6_TOP,Ra,ta)
    DATA[lc]=dict(ec=ec,em=em,side_mat=side_mat,top_mat=top_mat,**s)
    d=OUT/f"attr_{lc}_emax.defattr"; write_defattr(d,"emaxlof",em,rc)
    cxc(OUT/f"{lc}_emax.cxc",s,d,"emaxlof",PAL["emax"],EMAX_MAX,f"{Lig} Emax LoF",CAM_SIDE,side_mat,"side")
    d=OUT/f"attr_{lc}_ec50.defattr"; write_defattr(d,"ec50lof",ec,rc)
    cxc(OUT/f"{lc}_ec50.cxc",    s,d,"ec50lof",PAL["ec50"],EC50_MAX,f"{Lig} EC50 LoF",CAM_SIDE,side_mat,"side")
    cxc(OUT/f"{lc}_ec50_top.cxc",s,d,"ec50lof",PAL["ec50"],EC50_MAX,f"{Lig} EC50 LoF",CAM_TOP,top_mat,"top")
    print(f"{Lig:9s}: chain {rc}  [{s['cif']}]  ec/em positions {len(ec)}/{len(em)}")

master=OUT/"render_all.cxc"; lines=["# Render all maps to PNG"]
for lc in LIG:
    for name in [f"{lc}_emax",f"{lc}_ec50",f"{lc}_ec50_top"]:
        lines+=[f"open {OUT}/{name}.cxc", f"save {OUT}/{name}.png width 1200 height 1500 supersample 3"]
master.write_text("\n".join(lines)+"\n"); print("master:",master.name)

# ---------- INTERACTIVE WORKSPACE: all 3 structures superposed, both attributes ----------
MODELS=[("morphine",1),("fentanyl",2),("damgo",3)]   # fixed open order -> model ids
for param,attr in [("emax","emaxlof"),("ec50","ec50lof")]:
    with open(OUT/f"attr_all_{param}.defattr","w") as f:
        f.write(f"attribute: {attr}\nrecipient: residues\n")
        for lc,mid in MODELS:
            series=DATA[lc]["em"] if param=="emax" else DATA[lc]["ec"]
            ch=DATA[lc]["rchain"]
            for pos,val in series.items(): f.write(f"\t#{mid}/{ch}:{int(pos)}\t{val:.4f}\n")

def mats_line(key):   # "#1,<12> #2,<12> #3,<12>"
    return " ".join(f"#{mid},{DATA[lc][key]}" for lc,mid in MODELS)

ws=OUT/"workspace.cxc"
lig_lines=[]; color_lines=[]
for lc,mid in MODELS:
    d=DATA[lc]; rc=d["rchain"]; lig=d["lig"]; pep=(d["ligname"] is None)
    lig_lines+=[f"show #{mid}{lig if lig.startswith(':') else lig} atoms",
                f"style #{mid}{lig} {'stick' if pep else 'sphere'}",
                f"color #{mid}{lig} yellow", f"color #{mid}{lig} byhetero"]
    color_lines.append(f"color byattribute r:emaxlof #{mid}/{rc} palette {PAL['emax']} range 0,{EMAX_MAX} target c noValueColor #cccccc")
ws.write_text(f"""# Interactive workspace: morphine(#1,8EF6/R)  fentanyl(#2,8EF5/M)  damgo(#3,8EFQ/R)
# All receptors superposed; both attributes loaded. Switch views: `view side` / `view top`.
close session
open {CIFD/'8ef6.cif'}
open {CIFD/'8ef5.cif'}
open {CIFD/'8efq.cif'}
delete #1 & ~/R
delete #2 & ~/M
delete #3 & ~/R & ~/P
delete :CLR
delete solvent
hide atoms
show #1/R cartoons
show #2/M cartoons
show #3/R cartoons
color #1/R #cccccc
color #2/M #cccccc
color #3/R #cccccc
{chr(10).join(lig_lines)}
open {OUT/'attr_all_emax.defattr'}
open {OUT/'attr_all_ec50.defattr'}
{chr(10).join(color_lines)}
set bgColor white
lighting soft
graphics silhouettes true
# --- named views (camera + all three model matrices) ---
view matrix camera {CAM_TOP}
view matrix models {mats_line('top_mat')}
view name top
view matrix camera {CAM_SIDE}
view matrix models {mats_line('side_mat')}
view name side
""")
print("workspace:",ws.name)

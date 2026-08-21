#!/usr/bin/env python3
"""
Build the cAMP + TruPath parameter table (LogEC50, Emax %DAMGO, each +/- SEM, n).

CENTRAL values = the exact estimates used in the figures (so table <-> figures agree):
  cAMP   LogEC50/Emax : pooled 3-param Hill=1 fit (Emax = span/DAMGO_span*100)
  TruPath LogEC50      : pooled 3-param Hill=1 fit
  TruPath Emax         : Prism 'Data 4' peak-based value
SEM = across independent replicate curves (per-replicate Hill=1 fits); TruPath Emax
      SEM = Data 4. n = number of replicate curves.
"""
import re, csv, warnings, numpy as np, openpyxl
from scipy.optimize import curve_fit
warnings.filterwarnings("ignore")
def model(x,top,bot,le): return bot+(top-bot)/(1+10.0**(le-x))
def fit(pts):
    x=np.array([p[0] for p in pts]); y=np.array([p[1] for p in pts])
    p,_=curve_fit(model,x,y,p0=[np.min(y),np.max(y),-8.0],maxfev=40000); return p  # top,bot,le

# ---------------- cAMP ----------------
CAMP_BLOCKS={2:"DAMGO",11:"Fentanyl",20:"Methadone",29:"Buprenorphine",38:"Morphine",
   47:"Naloxone",56:"PZM21",65:"Uninduced cells",74:"TRV130",83:"MP",92:"Carfentanil",101:"Nalbuphine",110:"Butorphanol"}
ws=openpyxl.load_workbook("transcriptional_drc.xlsx",data_only=True)["Sheet1"]
starts=sorted(CAMP_BLOCKS)
camp_pool={}; camp_rep={}; camp_nrep={}
for i,c in enumerate(starts):
    nxt=starts[i+1] if i+1<len(starts) else 119; nm=CAMP_BLOCKS[c]
    if nm=="Uninduced cells": continue
    allpts=[]; cols=[]; ncol=0
    for cc in range(c,nxt):
        pts=[]
        for r in range(2,ws.max_row+1):
            xa=ws.cell(r,1).value; v=ws.cell(r,cc).value
            if xa in (None,"") or v in (None,""): continue
            if isinstance(v,str) and "*" in v: continue
            try: pts.append((float(xa),float(str(v).replace("*","").strip())))
            except: continue
        allpts+=pts
        if len(pts)>=5:
            ncol+=1                                   # data replicate curve
            try: cols.append(tuple(fit(pts)))
            except: pass
    camp_pool[nm]=fit(allpts); camp_rep[nm]=cols; camp_nrep[nm]=ncol
camp_pool_dspan=camp_pool["DAMGO"][1]-camp_pool["DAMGO"][0]
camp_rep_dspan=np.mean([b-t for (t,b,l) in camp_rep["DAMGO"]])

# ---------------- TruPath ----------------
PZFX="20260117_mor_wt_trupath_compilation.pzfx"
txt=open(PZFX,encoding="utf-8",errors="replace").read()
strip=lambda s:re.sub(r"<[^>]+>","",s).strip()
tabs=re.split(r"(?=<Table )",txt)
def get(t): return next(tb for tb in tabs if re.search(r"<Title>(.*?)</Title>",tb,re.S) and strip(re.search(r"<Title>(.*?)</Title>",tb,re.S).group(1))==t)
def subs(cb):
    return [[(strip(d) if d!="" else "") for d in re.findall(r"<d[^>]*/>|<d[^>]*>(.*?)</d>",s,re.S)]
            for s in re.findall(r"<Subcolumn>(.*?)</Subcolumn>",cb,re.S)]
rt=get("20260117 µOR WT Gi1 TRUPATH")
Xr=subs(re.search(r"<XColumn.*?</XColumn>",rt,re.S).group(0))[0]
EXCLUDE={("MP",-4.0)}
d4=get("Data 4"); D4={}
for yc in re.findall(r"<YColumn.*?</YColumn>",d4,re.S):
    nm=strip(re.search(r"<Title>(.*?)</Title>",yc,re.S).group(1))
    vals=[strip(x) for x in re.findall(r"<d[^>]*>(.*?)</d>",yc,re.S) if strip(x)!=""]
    if len(vals)>=2: D4[nm]=(float(vals[0]),float(vals[1]))
tp_pool={}; tp_rep={}; tp_nrep={}
for yc in re.findall(r"<YColumn.*?</YColumn>",rt,re.S):
    nm=strip(re.search(r"<Title>(.*?)</Title>",yc,re.S).group(1))
    if nm not in D4: continue
    allpts=[]; les=[]; nsub=0
    for sub in subs(yc):
        pts=[]
        for k,v in enumerate(sub):
            if k<len(Xr) and v!="" and Xr[k]!="" and "*" not in v:
                try: x=float(Xr[k]); y=float(v)
                except: continue
                if (nm,x) in EXCLUDE: continue
                pts.append((x,y))
        allpts+=pts
        if len(pts)>=5:
            nsub+=1                                   # data replicate curve
            try:
                p=fit(pts)
                if -14<p[2]<-3: les.append(p[2])
            except: pass
    tp_pool[nm]=fit(allpts); tp_rep[nm]=les; tp_nrep[nm]=nsub

def sem(a): a=np.array(a); return a.std(ddof=1)/np.sqrt(len(a)) if len(a)>1 else float('nan')

LIGS=["DAMGO","Morphine","Carfentanil","Fentanyl","C6Guano","Methadone","PZM21","TRV130",
      "MP","Butorphanol","Buprenorphine","Nalbuphine","Naltrexone","Naloxone"]
rows=[]
for nm in LIGS:
    row=dict(ligand=nm)
    # cAMP
    if nm in camp_pool:
        t,b,le=camp_pool[nm]
        row["camp_logec50"]=round(le,2)
        row["camp_emax"]=round((b-t)/camp_pool_dspan*100,1)
        row["camp_logec50_sem"]=round(sem([l for (tt,bb,l) in camp_rep[nm]]),2)
        row["camp_emax_sem"]=round(sem([(bb-tt)/camp_rep_dspan*100 for (tt,bb,l) in camp_rep[nm]]),1)
        row["camp_n"]=camp_nrep[nm]
    else:
        for k in ["camp_logec50","camp_emax","camp_logec50_sem","camp_emax_sem"]: row[k]=""
        row["camp_n"]=""
    # TruPath
    t,b,le=tp_pool[nm]
    row["tp_logec50"]=round(le,2)
    row["tp_logec50_sem"]=round(sem(tp_rep[nm]),2)
    row["tp_emax"]=round(D4[nm][0],1)
    row["tp_emax_sem"]=round(D4[nm][1],1)
    row["tp_n"]=tp_nrep[nm]
    rows.append(row)

cols=["ligand","camp_logec50","camp_logec50_sem","camp_emax","camp_emax_sem","camp_n",
      "tp_logec50","tp_logec50_sem","tp_emax","tp_emax_sem","tp_n"]
with open("param_table.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=cols); w.writeheader(); w.writerows(rows)

def f(v,d): return ("%%.%df"%d)%v if v!="" else "--"
print(f"{'Ligand':13s} | {'cAMP LogEC50':>14} {'cAMP Emax%':>13} {'n':>3} | {'TP LogEC50':>13} {'TP Emax%':>13} {'n':>3}")
print("-"*95)
for r in rows:
    ce=f"{f(r['camp_logec50'],2)}±{f(r['camp_logec50_sem'],2)}" if r['camp_logec50']!="" else "--"
    cm=f"{f(r['camp_emax'],1)}±{f(r['camp_emax_sem'],1)}" if r['camp_emax']!="" else "--"
    te=f"{f(r['tp_logec50'],2)}±{f(r['tp_logec50_sem'],2)}"
    tm=f"{f(r['tp_emax'],1)}±{f(r['tp_emax_sem'],1)}"
    print(f"{r['ligand']:13s} | {ce:>14} {cm:>13} {str(r['camp_n']):>3} | {te:>13} {tm:>13} {r['tp_n']:>3}")
print("\nSaved param_table.csv")

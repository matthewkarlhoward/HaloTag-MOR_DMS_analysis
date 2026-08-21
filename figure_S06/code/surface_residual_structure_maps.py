#!/usr/bin/env python3
"""
Expression-corrected per-position LoF, painted on the morphine-bound MOR (8EF6).

Both structure panels in the existing figure show raw LoF, which surface
expression alone can explain a large part of (position-level rho: EC50 -0.52,
Emax +0.80 against Surface_effect). Here each metric is regressed on the
per-position surface score with LOWESS and the RESIDUAL is mapped instead, so
the colour is the LoF that expression does NOT account for:

  resid > 0   more loss than a receptor of that abundance should show
  resid <= 0  fully explained by expression -> clipped to white

Metric definitions and filters match curve_refitting/make_chimerax_param_maps.py
so the residual maps sit alongside the raw ones:
  EC50 LoF = mean position logEC50 - synonymous logEC50   (sigmoid, no_baseline)
  Emax LoF = 1 - mean position activity                   (sigmoid, flat, no_baseline)
  surface  = mean missense Surface_effect at that position
  positions need >= MIN_N scored variants for that metric

Structure/camera are the captured 8EF6 chain-R side view, ligand MOI as yellow
spheres, EC50 white->blue and Emax white->red, as in the raw maps.

Outputs -> plots/scatter/surface/chimerax/
  attr_morphine_{ec50,emax}_resid.defattr
  morphine_{ec50,emax}_resid.cxc          (+ render_residual_maps.cxc)
  surface_residuals_by_position.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

ROOT = Path(__file__).resolve().parents[3]
CR = ROOT / "curve_refitting"
CIFD = ROOT / "structures" / "raw" / "experimental"
OUT = Path(__file__).resolve().parent / "chimerax"
OUT.mkdir(parents=True, exist_ok=True)

DRUG = "morphine"
CIF, RCHAIN, LIG = "8ef6.cif", "R", ":MOI"
MIN_N = 5
REF = {"sigmoid", "no_baseline"}
EC50_KEEP = {"sigmoid", "no_baseline"}
EMAX_KEEP = {"sigmoid", "flat", "no_baseline"}
LOESS_FRAC, LOESS_IT = 0.3, 3

# Captured 8EF6 chain-R side view (same as the raw maps).
CAM_SIDE = ("-0.66813,0.056198,-0.74192,38.605,0.20656,0.97196,-0.11239,102.95,"
            "0.7148,-0.22835,-0.661,35.991")
MOD_SIDE = ("-0.85867,-0.19853,0.47252,221.61000,"
            "0.20054,0.71828,0.66622,-76.67700,"
            "-0.47166,0.66683,-0.57695,181.08000")

# ── Per-position LoF metrics ─────────────────────────────────────────────────
fit = pd.read_csv(CR / f"refit_3param_robust_{DRUG}.csv")
syn = fit[(fit.type == "synonymous") & (fit.new_curve_type.isin(REF))]
ec_syn = syn.fitted_ec50_logM.mean()
e0, e1 = syn.fitted_emin.mean(), syn.fitted_emax.mean()
mis = fit[fit.type == "missense"]

ec = (mis[mis.new_curve_type.isin(EC50_KEEP)].dropna(subset=["fitted_ec50_logM"])
      .groupby("position").fitted_ec50_logM.agg(["mean", "count"]))
ec50_lof = (ec[ec["count"] >= MIN_N]["mean"] - ec_syn).rename("ec50_lof")

em = mis[mis.new_curve_type.isin(EMAX_KEEP)].copy()
em["loss"] = 1 - (em.fitted_emax - e0) / (e1 - e0)
emg = em.groupby("position").loss.agg(["mean", "count"])
emax_lof = emg[emg["count"] >= MIN_N]["mean"].rename("emax_lof")

# ── Per-position surface expression ──────────────────────────────────────────
dms = pd.read_csv(ROOT / "dms_scores" / "composite_dms_scores.csv")
surface = (dms[dms.type == "missense"].groupby("position")["Surface_effect"]
           .mean().rename("surface"))

df = pd.concat([surface, ec50_lof, emax_lof], axis=1)

# ── LOWESS residuals ─────────────────────────────────────────────────────────
def add_residuals(frame, col):
    """LOWESS-predicted LoF from surface score, and the residual around it."""
    sub = frame[["surface", col]].dropna()
    sm = lowess(sub[col], sub.surface, frac=LOESS_FRAC, it=LOESS_IT,
                return_sorted=True)
    pred = np.interp(sub.surface, sm[:, 0], sm[:, 1])
    frame.loc[sub.index, f"{col}_pred"] = pred
    frame.loc[sub.index, f"{col}_resid"] = sub[col] - pred
    # Linear residual kept for comparison; only the LOWESS one is painted.
    b, a = np.polyfit(sub.surface, sub[col], 1)
    frame.loc[sub.index, f"{col}_resid_lin"] = sub[col] - (a + b * sub.surface)


for col in ("ec50_lof", "emax_lof"):
    add_residuals(df, col)

df.index.name = "position"
df.to_csv(Path(__file__).resolve().parent / "surface_residuals_by_position.csv")

# ── ChimeraX assets ──────────────────────────────────────────────────────────
def write_defattr(path, attr, series):
    with open(path, "w") as f:
        f.write(f"attribute: {attr}\nrecipient: residues\n")
        for pos, val in series.dropna().items():
            f.write(f"\t/{RCHAIN}:{int(pos)}\t{val:.4f}\n")


def write_cxc(path, defattr, attr, palette, vmax, title):
    path.write_text(f"""# {title} (side view)
close session
open {CIFD / CIF}
delete ~/{RCHAIN}
delete :CLR
delete solvent
hide atoms
show /{RCHAIN} cartoons
color /{RCHAIN} #cccccc
show {LIG} atoms
style {LIG} sphere
color {LIG} yellow
color {LIG} byhetero
open {defattr}
color byattribute r:{attr} /{RCHAIN} palette {palette} range 0,{vmax} target c noValueColor #cccccc
set bgColor white
lighting soft
graphics silhouettes true
windowsize 1600 1200
view matrix camera {CAM_SIDE}
view matrix models #1,{MOD_SIDE}
""")


PANELS = [
    dict(key="ec50", attr="ec50resid", palette="white:white:blue",
         title="Morphine EC50 LoF, expression-corrected"),
    dict(key="emax", attr="emaxresid", palette="white:white:red",
         title="Morphine Emax LoF, expression-corrected"),
]

summary = []
for p in PANELS:
    resid = df[f"{p['key']}_lof_resid"].dropna()
    # Only excess loss is painted; gains relative to the trend go to white.
    painted = resid.clip(lower=0)
    vmax = float(np.round(np.nanpercentile(resid, 97.5), 2))
    d = OUT / f"attr_{DRUG}_{p['key']}_resid.defattr"
    write_defattr(d, p["attr"], painted)
    write_cxc(OUT / f"{DRUG}_{p['key']}_resid.cxc", d, p["attr"],
              p["palette"], vmax, p["title"])
    summary.append(dict(panel=p["key"], n_positions=len(resid),
                        resid_min=resid.min(), resid_max=resid.max(),
                        vmax_97p5=vmax, n_above_vmax=int((resid > vmax).sum())))

master = OUT / "render_residual_maps.cxc"
master.write_text("\n".join(
    line for p in PANELS for line in
    (f"open {OUT}/{DRUG}_{p['key']}_resid.cxc",
     f"save {OUT}/{DRUG}_{p['key']}_resid.png width 1200 height 1500 "
     f"supersample 3 transparentBackground true")
) + "\n")

print(pd.DataFrame(summary).to_string(index=False))
print("\ncorrelation of residual with surface score (should be ~0):")
for p in PANELS:
    sub = df[["surface", f"{p['key']}_lof_resid"]].dropna()
    print(f"  {p['key']:5s} r = {np.corrcoef(sub.surface, sub.iloc[:, 1])[0, 1]:+.3f}")
print(f"\nWrote {OUT}")

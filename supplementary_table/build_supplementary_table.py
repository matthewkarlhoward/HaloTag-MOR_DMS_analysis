#!/usr/bin/env python3
"""
Build the consolidated supplementary data table for the muOR deep mutational
pharmacology paper.

Reads only from ../data/ (this repository), writes:
  MOR_DMS_supplementary_data.xlsx   one workbook, one sheet per grain
  csv/<sheet>.csv                   the same sheets as flat CSVs

Grains are kept separate on purpose: variants, positions, ligands and sorted
samples are different units of observation and do not belong in one flat file.

Run:  python3 build_supplementary_table.py
"""
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
CSV_OUT = HERE / "csv"
CSV_OUT.mkdir(parents=True, exist_ok=True)
XLSX = HERE / "MOR_DMS_supplementary_data.xlsx"

# ---------------------------------------------------------------------------
# Conventions
# ---------------------------------------------------------------------------
# The 15 ligands screened in the paper. MCAM is present in the analysis master
# table but was NOT part of this study and is excluded everywhere below.
LIGANDS = ["Naltrexone", "Naloxone", "Nalbuphine", "Buprenorphine", "Butorphanol",
           "MP", "TRV130", "PZM21", "Methadone", "C6guano", "Morphine",
           "SR17018", "Fentanyl", "Carfentanil", "DAMGO"]
EXCLUDE_LIGANDS = ["MCAM"]

EFFICACY_CLASS = {
    "Naltrexone": "antagonist", "Naloxone": "antagonist",
    "Nalbuphine": "weak", "Buprenorphine": "weak", "Butorphanol": "weak", "MP": "weak",
    "TRV130": "intermediate", "PZM21": "intermediate", "Methadone": "intermediate",
    "C6guano": "strong", "Morphine": "strong", "SR17018": "strong",
    "Fentanyl": "strong", "Carfentanil": "strong", "DAMGO": "strong",
}
CHEMOTYPE = {
    "Naltrexone": "morphinan", "Naloxone": "morphinan", "Nalbuphine": "morphinan",
    "Buprenorphine": "morphinan", "Butorphanol": "morphinan", "Morphine": "morphinan",
    "MP": "indole alkaloid", "TRV130": "synthetic", "PZM21": "synthetic",
    "SR17018": "synthetic", "Methadone": "diphenylheptane",
    "C6guano": "fentanyl", "Fentanyl": "fentanyl", "Carfentanil": "fentanyl",
    "DAMGO": "peptide",
}
# name normalisation applied across the analysis repo
RENAME = {"Oliceridine": "TRV130", "Lofentanil": "Carfentanil",
          "Mitragynine_Pseudoindoxyl": "MP", "C6-guano": "C6guano",
          "C6Guano": "C6guano"}

DRC_LIGANDS = ["Morphine", "Fentanyl", "DAMGO"]

# (composite column, ligand, log10 [ligand] M).  DAMGO_mor_dms_6 is the
# ligand-free (forskolin-only) arm of the DAMGO screen, NOT a DAMGO dose.
DOSE_COLS = (
    [(f"MOR_Morphine_{t}M_rescaled_score", "Morphine", c) for t, c in
     [("-115", -11.5), ("-105", -10.5), ("-95", -9.5), ("-85", -8.5),
      ("-75", -7.5), ("-65", -6.5), ("-55", -5.5), ("-45", -4.5)]]
    + [(f"MOR_Fentanyl_{t}M_rescaled_score", "Fentanyl", c) for t, c in
       [("-12", -12.0), ("-11", -11.0), ("-10", -10.0), ("-9", -9.0),
        ("-8", -8.0), ("-7", -7.0), ("-6", -6.0), ("-5", -5.0)]]
    + [(f"DAMGO_mor_dms_{i}_rescaled_score", "DAMGO", c) for i, c in
       [(1, -5.0), (2, -6.0), (3, -7.0), (4, -8.0), (5, -9.0)]]
)
NO_LIGAND_COL = "DAMGO_mor_dms_6_rescaled_score"


def log(msg):
    print(msg, file=sys.stderr)


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
log("loading master table ...")
dms = pd.read_csv(DATA / "dms_scores/composite_dms_scores.csv", low_memory=False)

hill = {d: pd.read_csv(DATA / f"curve_fits/refit_3param_robust_{d.lower()}.csv")
        for d in DRC_LIGANDS}
oper = {d: pd.read_csv(DATA / f"operational_model/{d.lower()}_operational_fits.csv")
        for d in DRC_LIGANDS}

pca_var = pd.read_csv(DATA / "pca/pca_variant_scores.csv")
pca_pos = pd.read_csv(DATA / "pca/pca_position_scores.csv")
pca_load = pd.read_csv(DATA / "pca/pca_loadings_table.csv").rename(
    columns={"Unnamed: 0": "ligand"})
pc_table = pd.read_csv(DATA / "pca/per_position_pc_scores.csv")

lofgof = pd.read_csv(DATA / "dms_scores/lof_gof_scores.csv")
param = pd.read_csv(DATA / "pharmacology/param_table.csv")

# ---------------------------------------------------------------------------
# Sheet 01 — variant scores
# ---------------------------------------------------------------------------
log("sheet 01: variant scores ...")
ident = ["hgvs", "type", "position", "wildtype", "mutation",
         "GPCRdb", "SSE", "pos_mouse"]
v = dms[ident].copy()

# surface expression
v = v.join(dms[["Surface_effect", "Surface_effect_se", "Surface_lfsr"]].rename(
    columns={"Surface_effect": "surface_effect",
             "Surface_effect_se": "surface_effect_se",
             "Surface_lfsr": "surface_lfsr"}))

# 15-ligand saturating screen + the ligand-free (forskolin) arm
for lig in LIGANDS + ["FSK"]:
    name = "NoLigand" if lig == "FSK" else lig
    for src, dst in [("_effect", "_effect"), ("_effect_se", "_effect_se"),
                     ("_lfsr", "_lfsr")]:
        col = f"{lig}{src}"
        if col in dms.columns:
            v[f"{name}{dst}"] = dms[col]

# final Hill parameters (the published EC50 / Emax)
for d in DRC_LIGANDS:
    h = hill[d][["hgvs", "fitted_emin", "fitted_ec50_logM", "fitted_emax",
                 "r2", "n_points_used", "n_outliers_dropped", "new_curve_type",
                 "fit_ok"]].rename(columns={
        "fitted_emin": f"{d}_hill_emin",
        "fitted_ec50_logM": f"{d}_hill_logEC50_M",
        "fitted_emax": f"{d}_hill_emax",
        "r2": f"{d}_hill_r2",
        "n_points_used": f"{d}_hill_n_points",
        "n_outliers_dropped": f"{d}_hill_n_outliers_dropped",
        "new_curve_type": f"{d}_hill_curve_type",
        "fit_ok": f"{d}_hill_fit_ok"})
    v = v.merge(h, on="hgvs", how="left")

# Operational model parameters. Only morphine carries the expression
# decomposition (log10 rho); fentanyl and DAMGO were fit for concordance only.
OP_COLS = {"logtau": "op_logtau", "logKA": "op_logKA", "logEC50": "op_logEC50_M",
           "emax_obs_fit": "op_emax_obs", "log10_rho": "op_log10_rho",
           "rmse": "op_rmse", "saturating": "op_saturating",
           "ka_unidentified": "op_KA_unidentified"}
for d in DRC_LIGANDS:
    have = [c for c in OP_COLS if c in oper[d].columns]
    o = oper[d][["hgvs"] + have].rename(
        columns={c: f"{d}_{OP_COLS[c]}" for c in have})
    v = v.merge(o, on="hgvs", how="left")

# variant-level PCA scores
v = v.merge(pca_var[["hgvs", "PC1", "PC2", "PC3"]].rename(
    columns={"PC1": "PCA_variant_PC1", "PC2": "PCA_variant_PC2",
             "PC3": "PCA_variant_PC3"}), on="hgvs", how="left")

# external annotation
v = v.join(dms[["alphamissense_score", "alphamissense_class", "MTR"]])
v = v.sort_values(["position", "mutation"]).reset_index(drop=True)

# ---------------------------------------------------------------------------
# Sheet 02 — per-dose scores (concentration-response screens)
# ---------------------------------------------------------------------------
log("sheet 02: dose-response scores ...")
d2 = dms[["hgvs", "type", "position", "wildtype", "mutation"]].copy()
for col, lig, conc in DOSE_COLS:
    d2[f"{lig}_score_logM{conc:+.1f}"] = dms[col]
    d2[f"{lig}_sd_logM{conc:+.1f}"] = dms[col.replace("_score", "_sd")]
d2["NoLigand_DAMGOscreen_score"] = dms[NO_LIGAND_COL]
d2["NoLigand_DAMGOscreen_sd"] = dms[NO_LIGAND_COL.replace("_score", "_sd")]
d2 = d2.sort_values(["position", "mutation"]).reset_index(drop=True)

# ---------------------------------------------------------------------------
# Sheet 03 — position summary
# ---------------------------------------------------------------------------
log("sheet 03: position summary ...")
mis = dms[dms["type"] == "missense"]
pos = (dms.groupby("position")
          .agg(wildtype=("wildtype", "first"), GPCRdb=("GPCRdb", "first"),
               SSE=("SSE", "first"), pos_mouse=("pos_mouse", "first"))
          .reset_index())
pos["n_missense"] = mis.groupby("position").size().reindex(pos.position).values

for lig in LIGANDS:
    pos[f"{lig}_mean_missense_effect"] = (
        mis.groupby("position")[f"{lig}_effect"].mean().reindex(pos.position).values)
pos["NoLigand_mean_missense_effect"] = (
    mis.groupby("position")["FSK_effect"].mean().reindex(pos.position).values)
pos["surface_mean_missense_effect"] = (
    mis.groupby("position")["Surface_effect"].mean().reindex(pos.position).values)

for d in DRC_LIGANDS:
    h = hill[d]
    hm = h[h["type"] == "missense"]
    pos[f"{d}_mean_logEC50_M"] = (
        hm.groupby("position")["fitted_ec50_logM"].mean().reindex(pos.position).values)
    pos[f"{d}_mean_emax"] = (
        hm.groupby("position")["fitted_emax"].mean().reindex(pos.position).values)

pos = pos.merge(pca_pos[["position", "PC1", "PC2", "PC3"]].rename(
    columns={"PC1": "PCA_position_PC1", "PC2": "PCA_position_PC2",
             "PC3": "PCA_position_PC3"}), on="position", how="left")
pos = pos.merge(pc_table[["position", "variant_PC1", "variant_PC2", "variant_PC3"]].rename(
    columns={"variant_PC1": "PCA_variantmean_PC1",
             "variant_PC2": "PCA_variantmean_PC2",
             "variant_PC3": "PCA_variantmean_PC3"}), on="position", how="left")

# morphine-vs-fentanyl LOF bias (figures 4c, 4e)
for f, keep, pref in [
        ("pharmacology/mor_fent_ec50_reweighting_per_position_latest.csv",
         ["diff_MorMinusFen_logM", "q", "sig"], "EC50bias"),
        ("pharmacology/mor_fent_reweighting_per_position_latest.csv",
         ["diff_MorMinusFen%", "q", "sig"], "Emaxbias")]:
    p = DATA / f
    if p.exists():
        t = pd.read_csv(p)
        t = t[["position"] + [c for c in keep if c in t.columns]]
        t.columns = ["position"] + [f"{pref}_{c}" for c in t.columns[1:]]
        pos = pos.merge(t, on="position", how="left")

# G-protein contact class (active / intermediate / both)
gp_dir = DATA / "structures/processed/gprotein_distances"
act = set()
for f in sorted(gp_dir.glob("*_gprotein_distances.csv")):
    if "9PXW" in f.name:
        continue
    t = pd.read_csv(f)
    act |= set(t.loc[(t.shell == "first_shell") &
                     t.residue_number_human.notna(), "residue_number_human"].astype(int))
inter_f = DATA / "structures/processed/intermediate_gprotein_distances/intermediate_contacts_combined.csv"
t = pd.read_csv(inter_f)
inter = set(t.loc[(t.shell == "first_shell") & (t.partner_type == "galpha") &
                  t.residue_number_human.notna(), "residue_number_human"].astype(int))


def gp_class(p):
    a, i = p in act, p in inter
    return "Both" if a and i else "Active" if a else "Intermediate" if i else ""


pos["gprotein_contact_class"] = pos.position.map(gp_class)
pos = pos.sort_values("position").reset_index(drop=True)

# ---------------------------------------------------------------------------
# Sheet 04 — per-position LOF / GOF fractions (long)
# ---------------------------------------------------------------------------
log("sheet 04: LOF/GOF ...")
lg = lofgof.rename(columns={"drug": "ligand"}).copy()
lg["ligand"] = lg["ligand"].replace(RENAME)
lg = lg[~lg["ligand"].isin(EXCLUDE_LIGANDS)]
lg["efficacy_class"] = lg["ligand"].map(EFFICACY_CLASS)
lg = lg[["position", "ligand", "efficacy_class", "n_variants",
         "pct_lof", "pct_gof", "lof_thresh", "gof_thresh"]]
lg = lg.sort_values(["position", "ligand"]).reset_index(drop=True)

# ---------------------------------------------------------------------------
# Sheet 05 — ligand summary
# ---------------------------------------------------------------------------
log("sheet 05: ligand summary ...")
lig_tbl = pd.DataFrame({"ligand": LIGANDS})
lig_tbl["chemotype"] = lig_tbl.ligand.map(CHEMOTYPE)
lig_tbl["efficacy_class"] = lig_tbl.ligand.map(EFFICACY_CLASS)
lig_tbl["dms_screen"] = np.where(lig_tbl.ligand.isin(DRC_LIGANDS),
                                 "concentration-response + 10 uM panel", "10 uM panel")

p = param.copy()
p["ligand"] = p["ligand"].replace(RENAME)
lig_tbl = lig_tbl.merge(p, on="ligand", how="left")

# PCA ligand loadings, position-level fit
pl = pca_load.copy()
pl["ligand"] = pl["ligand"].str.replace("_effect$", "", regex=True).replace(RENAME)
lig_tbl = lig_tbl.merge(
    pl[["ligand", "PC1", "PC2", "PC3"]].rename(
        columns={"PC1": "PCA_position_loading_PC1", "PC2": "PCA_position_loading_PC2",
                 "PC3": "PCA_position_loading_PC3"}), on="ligand", how="left")

# PCA ligand loadings, variant-level fit (recomputed exactly as in figure 5d-f)
ec = [c for c in pca_var.columns if c.endswith("_effect")]
X = pca_var[ec].to_numpy(float)
keep = ~np.isnan(X).any(axis=1)
Xz = X[keep]
mu, sd = Xz.mean(0), Xz.std(0)
sd[sd == 0] = 1
Xz = (Xz - mu) / sd
U, S, Vt = np.linalg.svd(Xz - Xz.mean(0), full_matrices=False)
ld = pd.DataFrame((Vt.T * S)[:, :3], columns=["PC1", "PC2", "PC3"],
                  index=[c[:-7] for c in ec])
scores = (U * S)[:, :3]
# sign-align to the position PCA, the same way the figure script does
vp = pca_var.loc[keep, ["position"]].copy()
for j, col in enumerate(["PC1", "PC2", "PC3"]):
    vp[col] = scores[:, j]
    mm = pca_pos[["position", col]].merge(
        vp.groupby("position")[col].mean().rename("v").reset_index(), on="position").dropna()
    if np.corrcoef(mm[col], mm["v"])[0, 1] < 0:
        ld[col] *= -1
ld = ld.reset_index().rename(columns={"index": "ligand"})
ld["ligand"] = ld["ligand"].replace(RENAME)
lig_tbl = lig_tbl.merge(ld.rename(columns={
    "PC1": "PCA_variant_loading_PC1", "PC2": "PCA_variant_loading_PC2",
    "PC3": "PCA_variant_loading_PC3"}), on="ligand", how="left")

order = {c: i for i, c in enumerate(["antagonist", "weak", "intermediate", "strong"])}
lig_tbl = lig_tbl.sort_values(
    ["efficacy_class", "tp_emax"],
    key=lambda s: s.map(order) if s.name == "efficacy_class" else s).reset_index(drop=True)

# ---------------------------------------------------------------------------
# Sheet 06 — ligand contacts (long)
# ---------------------------------------------------------------------------
log("sheet 06: ligand contacts ...")
con = pd.read_csv(DATA / "structures/processed/ligand_distances/experimental_and_chai_combined.csv")
bup = DATA / "structures/processed/ligand_distances/per_pdb_legacy/buprenorphine_distances.csv"
if bup.exists():
    b = pd.read_csv(bup)
    b["drug_name"] = "Buprenorphine"
    con = pd.concat([con, b], ignore_index=True)
con["ligand"] = con["drug_name"].replace(RENAME)
con = con[con["ligand"].isin(LIGANDS)]
con = con.rename(columns={"residue_number_human": "position",
                          "pdb_id": "structure",
                          "distance_angstrom": "min_distance_angstrom"})
con = (con[["position", "ligand", "structure", "min_distance_angstrom",
            "distance_sidechain_angstrom", "shell"]]
       .dropna(subset=["position"]))
con["position"] = con["position"].astype(int)
con = (con.sort_values(["ligand", "position", "min_distance_angstrom"])
          .drop_duplicates(["ligand", "position"])
          .sort_values(["position", "ligand"]).reset_index(drop=True))

# ---------------------------------------------------------------------------
# Sheet 07 — validation pharmacology: fitted parameters (singles + doubles)
# ---------------------------------------------------------------------------
log("sheet 07: validation pharmacology (fitted) ...")
sys.path.insert(0, str(HERE))
try:
    from drc_fit import fit_drc
except ImportError:
    sys.path.insert(0, str(HERE.parent / "figure_08" / "code"))
    from drc_fit import fit_drc

# --- doubles: fit here, as plotted in figures 8d/f/g and S11a ---------------
pts = pd.read_csv(DATA / "pharmacology/doubles/doubles_merged_points.csv")


def per_rep_window_sem(g, keys, xcol, ycol):
    per = []
    for _, r in g.groupby(keys):
        m = r.groupby(xcol)[ycol].mean().reset_index()
        f = fit_drc(m[xcol].values, m[ycol].values)
        if f["ok"] and np.isfinite(f["span"]):
            per.append(-f["span"])
    return float(np.std(per, ddof=1) / np.sqrt(len(per))) if len(per) > 1 else np.nan


rows = []
for (lig, var), g in pts.groupby(["ligand", "variant"]):
    m = g.groupby("logM")["response"].mean().reset_index()
    r = fit_drc(m.logM.values, m.response.values)
    rows.append(dict(
        pipeline="doubles", figure="8d, 8f, 8g, S11a",
        source="doubles1+doubles2", run="20260421+20260427",
        ligand=lig, variant=var.replace(f"{lig}_", ""),
        is_double=bool(g.is_double.iloc[0]), partner=g.partner.iloc[0],
        n_replicates=int(g.groupby(["experiment", "rep"]).ngroups),
        n_concentrations=int(g.logM.nunique()),
        logEC50_M=r["params"][2] if r["params"] is not None else np.nan,
        activation_window=-r["span_used"],
        activation_window_sem=per_rep_window_sem(g, ["experiment", "rep"], "logM", "response"),
        responsive=r["responsive"], f_test_p=r["pval"], r2=r["r2"],
        normalization="two runs co-scaled on shared WT/A119L arms; not renormalised to DAMGO"))
dbl = pd.DataFrame(rows)
wt_damgo = dbl.loc[(dbl.ligand == "DAMGO") & (dbl.variant == "WT"), "activation_window"].iloc[0]
dbl["pct_activation_vs_WT_DAMGO"] = dbl.activation_window / wt_damgo * 100
dbl["pct_activation_sem"] = dbl.activation_window_sem / wt_damgo * 100

# --- singles: parameters as fitted by the singles pipeline ------------------
sp = pd.read_csv(DATA / "pharmacology/singles/singles_norm_params.csv")
spts = pd.read_csv(DATA / "pharmacology/singles/singles_norm_points.csv")
FIGMAP = {"ICL": "7d, 7e", "R278_run1": "7f", "doubles_raw": "7d-f (same raw run as doubles1)"}
sing = pd.DataFrame(dict(
    pipeline="singles", figure=sp.source.map(FIGMAP).fillna("7d-f"),
    source=sp.source, run=sp.run.astype(str),
    ligand=sp.ligand, variant=sp.variant,
    is_double=sp.variant.str.contains("_"), partner=np.nan,
    n_replicates=sp.n_rep, n_concentrations=sp.n_dose,
    logEC50_M=sp.logec50, activation_window=sp.window,
    activation_window_sem=np.nan,
    responsive=sp.responsive, f_test_p=sp.pval, r2=sp.r2,
    normalization="each curve divided by its own fitted no-drug plateau"))
# percent is taken within each source, against that source's own WT + DAMGO
ref = (sp[(sp.variant == "WT") & (sp.ligand == "DAMGO")]
       .set_index("source")["window"].to_dict())
sing["pct_activation_vs_WT_DAMGO"] = [
    w / ref[s] * 100 if s in ref else np.nan for w, s in zip(sing.activation_window, sing.source)]
sing["pct_activation_sem"] = np.nan

COLS = ["pipeline", "figure", "source", "run", "ligand", "variant", "is_double", "partner",
        "n_replicates", "n_concentrations", "logEC50_M", "activation_window",
        "activation_window_sem", "pct_activation_vs_WT_DAMGO", "pct_activation_sem",
        "responsive", "f_test_p", "r2", "normalization"]
val = pd.concat([sing[COLS], dbl[COLS]], ignore_index=True)
val["assay"] = "TRUPATH Gi1 BRET"
val = val.sort_values(["pipeline", "source", "ligand", "is_double", "variant"]).reset_index(drop=True)

# ---------------------------------------------------------------------------
# Sheet 09 — validation pharmacology: individual replicate points
# ---------------------------------------------------------------------------
log("sheet 09: validation points ...")
sraw = pd.read_csv(DATA / "pharmacology/singles/singles_raw_points.csv")
s_pts = spts.merge(sraw[["source", "run", "ligand", "variant", "logM", "rep", "bret"]],
                   on=["source", "run", "ligand", "variant", "logM", "rep", "bret"], how="left")
s_long = pd.DataFrame(dict(
    pipeline="singles", source=s_pts.source, run=s_pts.run.astype(str),
    ligand=s_pts.ligand, variant=s_pts.variant,
    is_double=s_pts.variant.str.contains("_"),
    log10_conc_M=s_pts.logM, replicate=s_pts.rep,
    raw_bret=s_pts.bret, normalized=s_pts.norm, divisor=s_pts.divisor))
d_long = pd.DataFrame(dict(
    pipeline="doubles", source=pts.set, run=pts.experiment.astype(str),
    ligand=pts.ligand, variant=[v.replace(f"{l}_", "") for v, l in zip(pts.variant, pts.ligand)],
    is_double=pts.is_double,
    log10_conc_M=pts.logM, replicate=pts.rep,
    raw_bret=np.nan, normalized=pts.response, divisor=np.nan))
vpts = pd.concat([s_long, d_long], ignore_index=True)
vpts = vpts.sort_values(["pipeline", "source", "ligand", "variant",
                         "log10_conc_M", "replicate"]).reset_index(drop=True)

# Restrict to what the paper actually plots: figure 7d-f (DAMGO / PZM21 /
# nalbuphine) and figures 8d,f,g + S11a (those three plus naloxone). The singles
# pipeline also measured five further ligands and two further I280 substitutions
# that appear in no figure; those are dropped here.
PUB_LIGANDS = {"DAMGO", "PZM21", "Nalbuphine", "Naloxone"}
PUB_DROP_VARIANTS = {"I280D", "I280P"}
val = val[val.ligand.isin(PUB_LIGANDS) & ~val.variant.isin(PUB_DROP_VARIANTS)].reset_index(drop=True)
vpts = vpts[vpts.ligand.isin(PUB_LIGANDS) & ~vpts.variant.isin(PUB_DROP_VARIANTS)].reset_index(drop=True)

# ---------------------------------------------------------------------------
# Sheet 08 — screen sample metadata
# ---------------------------------------------------------------------------
log("sheet 08: screen samples ...")
samp = pd.read_csv(DATA / "library_qc/per_bin_coverage_summary.csv")
samp = samp.rename(columns={"sample": "ligand_or_condition",
                            "conc": "log10_conc_M",
                            "bin": "sort_bin"})
samp = samp[["sample_id", "assay", "date", "replicate", "replicate_display",
             "ligand_or_condition", "log10_conc_M", "sort_bin",
             "mean_coverage", "median_coverage", "n_variants",
             "variant_counts_filename"]]

# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------
SHEETS = [
    ("01_variant_scores", v),
    ("02_dose_response_scores", d2),
    ("03_position_summary", pos),
    ("04_position_lof_gof", lg),
    ("05_ligand_summary", lig_tbl),
    ("06_ligand_contacts", con),
    ("07_validation_pharmacology", val),
    ("08_screen_samples", samp),
    ("09_validation_trupath_points", vpts),
]

for name, df in SHEETS:
    df.to_csv(CSV_OUT / f"{name}.csv", index=False)
    log(f"  {name:32s} {df.shape[0]:>7d} rows x {df.shape[1]:>3d} cols")

readme = pd.read_csv(HERE / "sheet_index.csv") if (HERE / "sheet_index.csv").exists() else None
with pd.ExcelWriter(XLSX, engine="openpyxl") as xl:
    if readme is not None:
        readme.to_excel(xl, sheet_name="00_README", index=False)
    for name, df in SHEETS:
        df.to_excel(xl, sheet_name=name, index=False, freeze_panes=(1, 0))
log(f"wrote {XLSX}")

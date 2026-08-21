#!/usr/bin/env python3
"""
Build a per-variant CSV of the morphine-vs-fentanyl Emax bias.

Per-position metric:
  bias = ΔEmax(morphine) − ΔEmax(fentanyl)
    where ΔEmax(drug) = mean missense Emax at position − syn mean Emax
                      (sigmoid-only, n ≥ 5 missense variants per position)

  Positive bias → morphine Emax more damaged than fentanyl
  Negative bias → fentanyl Emax more damaged than morphine

Each variant is mapped to its position's bias value (so every variant at a
given position carries the same bias).

Output: dms_scores/mor_fent_emax_bias_per_variant.csv
  Columns:
    hgvs                      — variant identifier
    position                  — residue position
    bias_emax_mor_minus_fent  — per-position bias (NaN if either ligand filtered)
    delta_emax_morphine       — per-position ΔEmax, Morphine
    delta_emax_fentanyl       — per-position ΔEmax, Fentanyl
    ca_min_to_ligand          — min Cα distance to {Mor, Fent} ligand heavy atoms
    outside_pocket            — True if ca_min_to_ligand > 7 Å for both ligands
"""
from pathlib import Path
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

ROOT    = Path(__file__).resolve().parents[1]
CIF_DIR = ROOT / "structures" / "raw" / "experimental"
SRC     = ROOT / "dms_scores" / "composite_dms_scores.csv"
OUT     = ROOT / "dms_scores" / "mor_fent_emax_bias_per_variant.csv"

MIN_VARIANTS = 5
DIST_THRESH  = 7.0  # Å

STRUCT = {
    "Morphine": dict(pdb="8ef6", ligand_code="MOI", receptor_chain="R"),
    "Fentanyl": dict(pdb="8ef5", ligand_code="7V7", receptor_chain="R"),
}

# ── Cα distance to nearest ligand heavy atom ─────────────────────────────────
def calpha_to_smallmol(pdb_id, ligand_code, receptor_chain):
    parser = MMCIFParser(QUIET=True)
    s = parser.get_structure(pdb_id, CIF_DIR / f"{pdb_id}.cif")
    lig = [a for model in s for chain in model for res in chain
           if res.get_resname() == ligand_code
           for a in res.get_atoms() if a.element != "H"]
    rows = []
    for model in s:
        for chain in model:
            if chain.id != receptor_chain: continue
            for res in chain:
                if not is_aa(res, standard=True) or "CA" not in res: continue
                ca = res["CA"]
                d  = min((ca - la) for la in lig)
                rows.append(dict(position=res.id[1], ca_dist=float(d)))
    return pd.DataFrame(rows)

dm  = calpha_to_smallmol(STRUCT["Morphine"]["pdb"],
                         STRUCT["Morphine"]["ligand_code"],
                         STRUCT["Morphine"]["receptor_chain"]
                         ).rename(columns={"ca_dist":"ca_mor"})
dfn = calpha_to_smallmol(STRUCT["Fentanyl"]["pdb"],
                         STRUCT["Fentanyl"]["ligand_code"],
                         STRUCT["Fentanyl"]["receptor_chain"]
                         ).rename(columns={"ca_dist":"ca_fen"})
dist = dm.merge(dfn, on="position", how="inner")
dist["ca_min_to_ligand"] = dist[["ca_mor","ca_fen"]].min(axis=1)
dist["outside_pocket"]   = dist["ca_min_to_ligand"] > DIST_THRESH

# ── Per-position ΔEmax for each drug ─────────────────────────────────────────
dms = pd.read_csv(SRC)
syn = dms[dms.type == "synonymous"]
mis = dms[dms.type == "missense"]

def per_position_demax(drug):
    emax_col = f"{drug}_emax"
    ct_col   = f"{drug}_curve_type"
    syn_emax = syn.loc[syn[ct_col]=="sigmoid", emax_col].mean()
    mis_sig  = mis[mis[ct_col]=="sigmoid"]
    g = mis_sig.groupby("position").agg(
        emax = (emax_col, "mean"),
        n    = (emax_col, "count"),
    ).dropna()
    g = g[g.n >= MIN_VARIANTS]
    return (g["emax"] - syn_emax).rename(f"delta_emax_{drug.lower()}")

mor = per_position_demax("Morphine")
fen = per_position_demax("Fentanyl")

per_pos = pd.concat([mor, fen], axis=1)
per_pos["bias_emax_mor_minus_fent"] = (
    per_pos[mor.name] - per_pos[fen.name]
)
per_pos = per_pos.reset_index()
per_pos = per_pos.merge(dist[["position","ca_min_to_ligand","outside_pocket"]],
                        on="position", how="left")

# ── Map onto every variant ───────────────────────────────────────────────────
out = dms[["hgvs","position"]].merge(per_pos, on="position", how="left")
out = out[["hgvs","position",
           "bias_emax_mor_minus_fent",
           "delta_emax_morphine",
           "delta_emax_fentanyl",
           "ca_min_to_ligand",
           "outside_pocket"]]

out.to_csv(OUT, index=False)
print(f"Saved: {OUT}")
print(f"  variant rows: {len(out):,}")
print(f"  with non-NaN bias: {out.bias_emax_mor_minus_fent.notna().sum():,}")
n_pos_outside = (per_pos.outside_pocket & per_pos.bias_emax_mor_minus_fent.notna()).sum()
print(f"  positions outside pocket (Cα > {DIST_THRESH:.0f} Å) with bias: {n_pos_outside}")
print()
print("Distribution of bias across positions (per-position, not per-variant):")
print(per_pos.bias_emax_mor_minus_fent.describe())

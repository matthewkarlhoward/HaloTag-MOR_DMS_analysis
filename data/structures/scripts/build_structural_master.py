"""
build_structural_master.py

Merges all per-position structural summary data into a single master dataframe
keyed on residue_number_human (one row per receptor position).
"""

from pathlib import Path
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
PROC = REPO_ROOT / "structures" / "processed"
OUT = PROC / "structural_master.csv"


# ---------------------------------------------------------------------------
# Circular statistics helpers (numpy-only, no scipy required)
# ---------------------------------------------------------------------------
def _circmean(x):
    """Circular mean of angles in degrees, range [-180, 180]."""
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    if len(x) == 0:
        return np.nan
    rad = np.radians(x)
    return float(np.degrees(np.arctan2(np.nanmean(np.sin(rad)),
                                       np.nanmean(np.cos(rad)))))


def _circstd(x):
    """Circular std of angles in degrees (Yamartino-style)."""
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    if len(x) <= 1:
        return np.nan
    rad = np.radians(x)
    S = np.nanmean(np.sin(rad))
    C = np.nanmean(np.cos(rad))
    R = np.sqrt(S**2 + C**2)
    R = min(R, 1.0)          # numerical guard
    return float(np.degrees(np.sqrt(-2.0 * np.log(R))))


# ---------------------------------------------------------------------------
# A. CA displacements
# ---------------------------------------------------------------------------
def load_ca_displacements():
    fpath = PROC / "ca_displacements" / "grand_average_ca_displacements.csv"
    df = pd.read_csv(fpath)
    df = df.rename(columns={
        "n_active_structures":  "ca_n_structures",
        "mean_ca_displacement": "ca_displacement_mean",
        "std_ca_displacement":  "ca_displacement_std",
        "min_ca_displacement":  "ca_displacement_min",
        "max_ca_displacement":  "ca_displacement_max",
    })
    # keep residue_name for backbone construction
    return df[["residue_number_human", "residue_name",
               "ca_n_structures", "ca_displacement_mean",
               "ca_displacement_std", "ca_displacement_min",
               "ca_displacement_max"]]


# ---------------------------------------------------------------------------
# B. Active vs inactive interactions
# ---------------------------------------------------------------------------
def load_interactions():
    fpath = PROC / "interactions" / "active_vs_inactive_interactions.csv"
    df = pd.read_csv(fpath)
    # keep all columns as-is; residue_number_human is the key
    return df


# ---------------------------------------------------------------------------
# C. Pi interactions (optional – file may not exist yet)
# ---------------------------------------------------------------------------
def load_pi_interactions():
    fpath = PROC / "pi_interactions" / "active_vs_inactive_pi_interactions.csv"
    try:
        df = pd.read_csv(fpath)
        return df
    except FileNotFoundError:
        print(f"WARNING: pi interactions file not found yet: {fpath}")
        return None


# ---------------------------------------------------------------------------
# D. Rotamers – pivot chi1 / chi2 to wide format
# ---------------------------------------------------------------------------
def load_rotamers():
    fpath = PROC / "rotamers" / "active_vs_inactive_chi_comparison.csv"
    df = pd.read_csv(fpath)

    rows = []
    for rn, grp in df.groupby("residue_number_human"):
        row = {"residue_number_human": rn}
        for chi in ["chi1", "chi2"]:
            sub = grp[grp["chi_angle"] == chi]
            if sub.empty:
                for col in ["active_mean", "inactive_mean", "delta",
                            "rotamer_shift", "active_dominant",
                            "inactive_dominant"]:
                    row[f"{chi}_{col}"] = np.nan
            else:
                s = sub.iloc[0]
                row[f"{chi}_active_mean"]       = s["active_circular_mean"]
                row[f"{chi}_inactive_mean"]     = s["inactive_circular_mean"]
                row[f"{chi}_delta"]             = s["abs_mean_difference"]
                row[f"{chi}_rotamer_shift"]     = s["rotamer_shift"]
                row[f"{chi}_active_dominant"]   = s["active_dominant_rotamer"]
                row[f"{chi}_inactive_dominant"] = s["inactive_dominant_rotamer"]
        rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# E. Backbone angles – aggregate to circular mean/std per position
# ---------------------------------------------------------------------------
def load_backbone_angles():
    fpath = PROC / "backbone_angles" / "all_structures_combined.csv"
    df = pd.read_csv(fpath)

    records = []
    for rn, grp in df.groupby("residue_number_human"):
        rec = {
            "residue_number_human": rn,
            "residue_name": grp["residue_name"].dropna().iloc[0]
                            if grp["residue_name"].notna().any() else np.nan,
        }
        for angle in ["phi", "psi", "omega"]:
            vals = grp[angle].dropna().values
            rec[f"{angle}_mean"] = _circmean(vals) if len(vals) > 0 else np.nan
            rec[f"{angle}_std"]  = _circstd(vals)  if len(vals) > 1 else np.nan
        records.append(rec)

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# F. G-protein distances
# ---------------------------------------------------------------------------
def load_gprotein_distances():
    fpath = PROC / "gprotein_distances" / "all_structures_combined.csv"
    df = pd.read_csv(fpath)

    records = []
    for rn, grp in df.groupby("residue_number_human"):
        dist_vals = grp["distance_sidechain_angstrom"].dropna()
        rec = {
            "residue_number_human": rn,
            "gprotein_dist_mean":   dist_vals.mean() if len(dist_vals) > 0 else np.nan,
            "gprotein_dist_min":    dist_vals.min()  if len(dist_vals) > 0 else np.nan,
            "gprotein_shell_any":   int((grp["shell"].str.lower() != "none").any()),
        }
        records.append(rec)

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# G. Ligand contacts
# ---------------------------------------------------------------------------
def load_ligand_contacts():
    fpath = PROC / "ligand_contacts" / "all_structures_contacts_5.0A.csv"
    df = pd.read_csv(fpath)

    records = []
    for rn, grp in df.groupby("residue_number_human"):
        dist_vals = grp["distance_sidechain_angstrom"].dropna()
        drugs     = sorted(grp["ligand_name"].dropna().unique().tolist())
        gpcrdb    = grp["GPCRdb"].dropna().iloc[0] if grp["GPCRdb"].notna().any() else np.nan
        sse       = grp["SSE"].dropna().iloc[0]    if grp["SSE"].notna().any()    else np.nan

        records.append({
            "residue_number_human":        rn,
            "ligand_contact_n_structures": grp["pdb_id"].nunique(),
            "ligand_contact_min_dist":     dist_vals.min() if len(dist_vals) > 0 else np.nan,
            "ligand_contact_drugs":        ";".join(drugs),
            "GPCRdb":                      gpcrdb,
            "SSE":                         sse,
        })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# H. Ligand distances – tested drugs
# ---------------------------------------------------------------------------
TESTED_DRUGS = [
    "DAMGO", "Fentanyl", "Morphine", "SR17018",
    "Mitragynine_Pseudoindoxyl", "Lofentanil", "Oliceridine",
    "PZM21", "Naloxone", "C6-guano",
]

# Map drug_name → safe column suffix
_DRUG_SUFFIX = {
    "DAMGO":                    "DAMGO",
    "Fentanyl":                 "Fentanyl",
    "Morphine":                 "Morphine",
    "SR17018":                  "SR17018",
    "Mitragynine_Pseudoindoxyl":"Mitragynine",
    "Lofentanil":               "Lofentanil",
    "Oliceridine":              "Oliceridine",
    "PZM21":                    "PZM21",
    "Naloxone":                 "Naloxone",
    "C6-guano":                 "C6guano",
}

_SHELL_ORDER = {"first": 0, "second": 1, "none": 2}


def _min_shell(series):
    """Return the closest shell category across structures."""
    vals = series.dropna().str.lower().unique()
    if len(vals) == 0:
        return "none"
    ordered = sorted(vals, key=lambda s: _SHELL_ORDER.get(s, 99))
    return ordered[0]


def load_ligand_distances_tested():
    fpath = PROC / "ligand_distances_tested_drugs" / "all_structures_combined.csv"
    df = pd.read_csv(fpath)

    # All unique positions
    all_rn = df["residue_number_human"].unique()
    result = pd.DataFrame({"residue_number_human": all_rn})

    for drug in TESTED_DRUGS:
        suffix = _DRUG_SUFFIX[drug]
        sub = df[df["drug_name"] == drug]
        if sub.empty:
            result[f"dist_{suffix}_min"]   = np.nan
            result[f"shell_{suffix}"]      = np.nan
            continue

        # min distance per position
        dist_agg = (sub.groupby("residue_number_human")["distance_sidechain_angstrom"]
                       .min()
                       .reset_index()
                       .rename(columns={"distance_sidechain_angstrom": f"dist_{suffix}_min"}))
        # min shell per position
        shell_agg = (sub.groupby("residue_number_human")["shell"]
                        .apply(_min_shell)
                        .reset_index()
                        .rename(columns={"shell": f"shell_{suffix}"}))

        result = result.merge(dist_agg,  on="residue_number_human", how="left")
        result = result.merge(shell_agg, on="residue_number_human", how="left")

    return result


# ---------------------------------------------------------------------------
# I. Chai distances – predicted drugs
# ---------------------------------------------------------------------------
CHAI_DRUGS = ["Naltrexone", "Nalbuphine", "Methadone", "Butorphanol"]


def load_chai_distances():
    fpath = PROC / "chai_distances" / "chai_predictions_combined.csv"
    df = pd.read_csv(fpath)

    all_rn = df["residue_number_human"].unique()
    result = pd.DataFrame({"residue_number_human": all_rn})

    for drug in CHAI_DRUGS:
        sub = df[df["drug_name"] == drug]
        if sub.empty:
            result[f"dist_{drug}_chai_min"] = np.nan
            continue

        dist_agg = (sub.groupby("residue_number_human")["distance_sidechain_angstrom"]
                       .min()
                       .reset_index()
                       .rename(columns={"distance_sidechain_angstrom": f"dist_{drug}_chai_min"}))
        result = result.merge(dist_agg, on="residue_number_human", how="left")

    return result


# ---------------------------------------------------------------------------
# Main merge
# ---------------------------------------------------------------------------
RECEPTOR_MIN = 1
RECEPTOR_MAX = 400  # OPRM1 is 400 aa; anything outside is G-protein / fusion tag


def build_backbone(sources):
    """Union of all residue_number_human values across all loaded DataFrames,
    filtered to receptor positions only (1–400)."""
    all_rn = set()
    for df in sources:
        if df is not None and "residue_number_human" in df.columns:
            rn = df["residue_number_human"].dropna().astype(int)
            all_rn.update(rn[(rn >= RECEPTOR_MIN) & (rn <= RECEPTOR_MAX)].tolist())
    backbone = pd.DataFrame({"residue_number_human": sorted(all_rn)})
    return backbone


# ---------------------------------------------------------------------------
# J. Full GPCRdb annotations
# ---------------------------------------------------------------------------
def load_gpcrdb_annotations():
    fpath = REPO_ROOT / "annotations" / "GPCRdb_OPRM1_table.csv"
    if not fpath.exists():
        print(f"  WARNING: GPCRdb annotations not found at {fpath}")
        return None
    df = pd.read_csv(fpath)
    # Columns: GPCRdb, pos, SSE, wt_aa
    df = df.rename(columns={"pos": "residue_number_human", "wt_aa": "wt_aa_gpcrdb"})
    df = df[["residue_number_human", "GPCRdb", "SSE", "wt_aa_gpcrdb"]].dropna(subset=["residue_number_human"])
    df["residue_number_human"] = df["residue_number_human"].astype(int)
    return df


def main():
    print("Loading source files...")

    ca_df         = load_ca_displacements()
    inter_df      = load_interactions()
    pi_df         = load_pi_interactions()          # may be None
    rot_df        = load_rotamers()
    bb_df         = load_backbone_angles()
    gp_df         = load_gprotein_distances()
    lig_cont_df   = load_ligand_contacts()
    lig_dist_df   = load_ligand_distances_tested()
    chai_dist_df  = load_chai_distances()
    gpcrdb_df     = load_gpcrdb_annotations()       # full GPCRdb table

    # Collect all non-None sources for backbone construction
    all_sources = [ca_df, inter_df, pi_df, rot_df, bb_df,
                   gp_df, lig_cont_df, lig_dist_df, chai_dist_df, gpcrdb_df]

    print("Building backbone of all positions...")
    master = build_backbone(all_sources)

    # ---- Merge A: CA displacements (also carries residue_name) ----
    master = master.merge(ca_df, on="residue_number_human", how="left")

    # ---- Merge E: Backbone angles (fill residue_name if still missing) ----
    bb_rn = bb_df[["residue_number_human", "residue_name",
                   "phi_mean", "phi_std",
                   "psi_mean", "psi_std",
                   "omega_mean", "omega_std"]]
    # Merge backbone angles; resolve residue_name conflict
    master = master.merge(
        bb_rn.rename(columns={"residue_name": "residue_name_bb"}),
        on="residue_number_human", how="left"
    )
    mask = master["residue_name"].isna()
    master.loc[mask, "residue_name"] = master.loc[mask, "residue_name_bb"]
    master = master.drop(columns=["residue_name_bb"])

    # ---- Merge B: Interactions ----
    master = master.merge(inter_df, on="residue_number_human", how="left")

    # ---- Merge C: Pi interactions (optional) ----
    if pi_df is not None:
        master = master.merge(pi_df, on="residue_number_human", how="left")

    # ---- Merge D: Rotamers ----
    master = master.merge(rot_df, on="residue_number_human", how="left")

    # ---- Merge F: G-protein distances ----
    master = master.merge(gp_df, on="residue_number_human", how="left")

    # ---- Merge G: Ligand contacts (drop GPCRdb/SSE — full table added below) ----
    lig_cont_df = lig_cont_df.drop(columns=["GPCRdb", "SSE"], errors="ignore")
    master = master.merge(lig_cont_df, on="residue_number_human", how="left")

    # ---- Merge H: Tested-drug distances ----
    master = master.merge(lig_dist_df, on="residue_number_human", how="left")

    # ---- Merge I: Chai distances ----
    master = master.merge(chai_dist_df, on="residue_number_human", how="left")

    # ---- Merge J: Full GPCRdb annotations ----
    if gpcrdb_df is not None:
        master = master.merge(gpcrdb_df, on="residue_number_human", how="left")

    # ---- Final sort ----
    master = master.sort_values("residue_number_human").reset_index(drop=True)

    # ---- Write output ----
    OUT.parent.mkdir(parents=True, exist_ok=True)
    master.to_csv(OUT, index=False)
    print(f"\nWrote {OUT}")

    # ---- Summary ----
    print(f"\nShape: {master.shape[0]} rows x {master.shape[1]} columns")
    print("\nColumn list:")
    for c in master.columns:
        print(f"  {c}")

    missing = (master.isna().mean() * 100).sort_values(ascending=False)
    missing_high = missing[missing > 10]
    if missing_high.empty:
        print("\nNo columns with >10% missing values.")
    else:
        print("\nColumns with >10% missing values:")
        for col, pct in missing_high.items():
            print(f"  {col:50s}  {pct:.1f}%")


if __name__ == "__main__":
    main()

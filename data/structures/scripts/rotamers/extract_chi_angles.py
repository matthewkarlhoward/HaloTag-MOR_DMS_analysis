#!/usr/bin/env python3
"""
Extract sidechain chi (χ) dihedral angles for all residues across all MOR structures.

For each structure:
  - Parse the verified single receptor chain
  - Calculate χ1, χ2, χ3, χ4 (where applicable) for each residue
  - Map to human residue numbering (mouse + 2)

Output: per-structure CSV files and a combined master CSV.

Chi angle definitions (IUPAC):
  χ1: N  - CA - CB - CG  (or CG1 for Val/Ile, SG for Cys, OG/OG1 for Ser/Thr)
  χ2: CA - CB - CG - CD  (or CD1, OD1, ND1, SD depending on residue)
  χ3: CB - CG - CD - NE  (or CE, OE1, NE2, SD depending on residue)
  χ4: CG - CD - NE - CZ  (or CE - SD - CG for Met)

Uses BioPython for structure parsing and numpy for dihedral calculation.
"""

import os
import csv
import math
import numpy as np
from collections import defaultdict

from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

from pathlib import Path
REPO_ROOT = Path(__file__).resolve().parents[3]
CIF_DIR = REPO_ROOT / "structures" / "raw" / "experimental"

# ── Structure definitions ────────────────────────────────────────────────────
INACTIVE_STRUCTURES = {
    "4dkl":  {"receptor_chain": "A", "species": "mouse"},
    "7ul4":  {"receptor_chain": "A", "species": "mouse"},
    "9PXU":  {"receptor_chain": "R", "species": "human"},
    "9bjk":  {"receptor_chain": "R", "species": "mouse"},
    "9MQI":  {"receptor_chain": "A", "species": "human"},
}

ACTIVE_STRUCTURES = {
    "5c1m":  {"receptor_chain": "A", "species": "mouse"},
    "6dde":  {"receptor_chain": "R", "species": "mouse"},
    "6ddf":  {"receptor_chain": "R", "species": "mouse"},
    "7sbf":  {"receptor_chain": "R", "species": "mouse"},
    "7scg":  {"receptor_chain": "D", "species": "mouse"},
    "7t2g":  {"receptor_chain": "R", "species": "mouse"},
    "7t2h":  {"receptor_chain": "D", "species": "mouse"},
    "7u2k":  {"receptor_chain": "D", "species": "mouse"},
    "7u2l":  {"receptor_chain": "D", "species": "mouse"},
    "8e0g":  {"receptor_chain": "A", "species": "mouse"},
    "8ef5":  {"receptor_chain": "R", "species": "human"},
    "8ef6":  {"receptor_chain": "R", "species": "human"},
    "8efb":  {"receptor_chain": "R", "species": "human"},
    "8efl":  {"receptor_chain": "R", "species": "human"},
    "8efo":  {"receptor_chain": "R", "species": "human"},
    "8efq":  {"receptor_chain": "R", "species": "human"},
    "8f7q":  {"receptor_chain": "R", "species": "human"},
    "8f7r":  {"receptor_chain": "R", "species": "human"},
    "8k9k":  {"receptor_chain": "R", "species": "human"},
    "8k9l":  {"receptor_chain": "R", "species": "human"},
    "8Y72":  {"receptor_chain": "R", "species": "human"},
    "8Y73":  {"receptor_chain": "R", "species": "human"},
    "9bqj":  {"receptor_chain": "D", "species": "mouse"},
    "9PY2":  {"receptor_chain": "R", "species": "human"},
    "9PY3":  {"receptor_chain": "R", "species": "human"},
    "9PY4":  {"receptor_chain": "R", "species": "human"},
    "9WST":  {"receptor_chain": "R", "species": "mouse"},
    "9WSW":  {"receptor_chain": "R", "species": "mouse"},
}


# ── Chi angle atom definitions per residue type ─────────────────────────────
# Each entry is a list of chi angles; each chi angle is a tuple of 4 atom names.
# Residue types without chi angles (Gly, Ala) are omitted.

CHI_DEFINITIONS = {
    "ARG": [
        ("N", "CA", "CB", "CG"),    # χ1
        ("CA", "CB", "CG", "CD"),   # χ2
        ("CB", "CG", "CD", "NE"),   # χ3
        ("CG", "CD", "NE", "CZ"),   # χ4
    ],
    "LYS": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD"),
        ("CB", "CG", "CD", "CE"),
        ("CG", "CD", "CE", "NZ"),
    ],
    "MET": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "SD"),
        ("CB", "CG", "SD", "CE"),
    ],
    "GLU": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD"),
        ("CB", "CG", "CD", "OE1"),
    ],
    "GLN": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD"),
        ("CB", "CG", "CD", "OE1"),
    ],
    "ASP": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "OD1"),
    ],
    "ASN": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "OD1"),
    ],
    "ILE": [
        ("N", "CA", "CB", "CG1"),
        ("CA", "CB", "CG1", "CD1"),
    ],
    "LEU": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD1"),
    ],
    "HIS": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "ND1"),
    ],
    "TRP": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD1"),
    ],
    "TYR": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD1"),
    ],
    "PHE": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD1"),
    ],
    "PRO": [
        ("N", "CA", "CB", "CG"),
        ("CA", "CB", "CG", "CD"),
    ],
    "THR": [
        ("N", "CA", "CB", "OG1"),
    ],
    "VAL": [
        ("N", "CA", "CB", "CG1"),
    ],
    "SER": [
        ("N", "CA", "CB", "OG"),
    ],
    "CYS": [
        ("N", "CA", "CB", "SG"),
    ],
}


def calc_dihedral(p0, p1, p2, p3):
    """
    Calculate dihedral angle (in degrees) from four 3D points.
    Uses the atan2 approach for numerical stability.
    Returns angle in range [-180, 180].
    """
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2

    # Normalize b1 so that it does not influence magnitude of vector
    b1 /= np.linalg.norm(b1)

    # v = projection of b0 onto plane perpendicular to b1
    v = b0 - np.dot(b0, b1) * b1
    # w = projection of b2 onto plane perpendicular to b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in the plane
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)

    return math.degrees(math.atan2(y, x))


def extract_chi_angles(structure, chain_id, species):
    """
    Extract chi angles for all residues in the specified chain.

    Returns:
        dict: {human_resnum: {"resname": str, "chi1": float|None, "chi2": float|None, ...}}
    """
    offset = 2 if species == "mouse" else 0
    results = {}

    model = structure[0]
    chain = model[chain_id]

    for residue in chain:
        if not is_aa(residue, standard=True):
            continue
        # Skip HETATMs (modified residues, etc.)
        het_flag = residue.get_id()[0]
        if het_flag != " ":
            continue

        auth_resnum = residue.get_id()[1]
        human_resnum = auth_resnum + offset
        resname = residue.get_resname().strip()

        entry = {
            "resname": resname,
            "chi1": None,
            "chi2": None,
            "chi3": None,
            "chi4": None,
        }

        if resname in CHI_DEFINITIONS:
            for i, atom_names in enumerate(CHI_DEFINITIONS[resname], start=1):
                try:
                    coords = []
                    for aname in atom_names:
                        atom = residue[aname]
                        coords.append(atom.get_vector().get_array())
                    angle = calc_dihedral(
                        np.array(coords[0]),
                        np.array(coords[1]),
                        np.array(coords[2]),
                        np.array(coords[3]),
                    )
                    entry[f"chi{i}"] = round(angle, 2)
                except KeyError:
                    # Atom missing (e.g., disordered residue, truncated sidechain)
                    entry[f"chi{i}"] = None

        results[human_resnum] = entry

    return results


def classify_rotamer_chi1(angle):
    """
    Classify χ1 angle into standard rotameric states.
    gauche+ (g+): roughly +60° (range ~0 to +120)
    trans (t):    roughly +180° (range ~120 to ~240 / equivalently -120 to +120 wrapping)
    gauche- (g-): roughly -60° (range ~-120 to 0)
    """
    if angle is None:
        return None
    # Normalize to [-180, 180]
    while angle > 180:
        angle -= 360
    while angle < -180:
        angle += 360

    if -120 <= angle < 0:
        return "g-"
    elif 0 <= angle < 120:
        return "g+"
    else:
        return "t"


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cif_dir = str(CIF_DIR)
    results_dir = str(REPO_ROOT / "structures" / "processed" / "rotamers")
    os.makedirs(results_dir, exist_ok=True)

    parser = MMCIFParser(QUIET=True)

    # Combine all structures
    all_structures = {}
    for pdb_id, info in INACTIVE_STRUCTURES.items():
        all_structures[pdb_id] = {**info, "state": "inactive"}
    for pdb_id, info in ACTIVE_STRUCTURES.items():
        all_structures[pdb_id] = {**info, "state": "active"}

    # Master data: {pdb_id: {human_resnum: {resname, chi1, chi2, chi3, chi4}}}
    master_data = {}
    all_residues = set()  # Track all human residue numbers seen

    print(f"Extracting chi angles from {len(all_structures)} structures...")
    print(f"  Inactive: {len(INACTIVE_STRUCTURES)}")
    print(f"  Active:   {len(ACTIVE_STRUCTURES)}")
    print()

    for pdb_id, info in sorted(all_structures.items()):
        chain_id = info["receptor_chain"]
        species = info["species"]
        state = info["state"]

        # Try common CIF filename patterns
        cif_file = None
        for pattern in [f"{pdb_id}.cif", f"{pdb_id.lower()}.cif", f"{pdb_id.upper()}.cif"]:
            candidate = os.path.join(cif_dir, pattern)
            if os.path.exists(candidate):
                cif_file = candidate
                break

        if cif_file is None:
            print(f"  WARNING: CIF file not found for {pdb_id}, skipping")
            continue

        try:
            structure = parser.get_structure(pdb_id, cif_file)
            chi_data = extract_chi_angles(structure, chain_id, species)
            master_data[pdb_id] = {"state": state, "data": chi_data}
            all_residues.update(chi_data.keys())

            # Count how many chi angles were extracted
            n_chi1 = sum(1 for v in chi_data.values() if v["chi1"] is not None)
            n_chi2 = sum(1 for v in chi_data.values() if v["chi2"] is not None)
            n_total = len(chi_data)
            print(f"  {pdb_id} ({state:8s}, chain {chain_id}): "
                  f"{n_total} residues, {n_chi1} χ1, {n_chi2} χ2")

        except Exception as e:
            print(f"  ERROR processing {pdb_id}: {e}")
            continue

    # ── Write per-structure CSV files ────────────────────────────────────────
    print(f"\nWriting per-structure CSV files to {results_dir}/")
    for pdb_id, entry in sorted(master_data.items()):
        state = entry["state"]
        chi_data = entry["data"]

        outfile = os.path.join(results_dir, f"{pdb_id}_chi_angles.csv")
        with open(outfile, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "residue_number_human", "residue_name",
                "chi1", "chi2", "chi3", "chi4",
                "chi1_rotamer", "state"
            ])
            for resnum in sorted(chi_data.keys()):
                d = chi_data[resnum]
                writer.writerow([
                    resnum,
                    d["resname"],
                    d["chi1"] if d["chi1"] is not None else "",
                    d["chi2"] if d["chi2"] is not None else "",
                    d["chi3"] if d["chi3"] is not None else "",
                    d["chi4"] if d["chi4"] is not None else "",
                    classify_rotamer_chi1(d["chi1"]) or "",
                    state,
                ])

    # ── Write master combined CSV ────────────────────────────────────────────
    sorted_residues = sorted(all_residues)
    sorted_pdbs = sorted(master_data.keys())

    master_file = os.path.join(results_dir, "all_chi_angles_combined.csv")
    print(f"Writing combined master CSV: {master_file}")
    with open(master_file, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["residue_number_human", "residue_name", "pdb_id", "state",
                   "chi1", "chi2", "chi3", "chi4", "chi1_rotamer"]
        writer.writerow(header)

        for resnum in sorted_residues:
            for pdb_id in sorted_pdbs:
                entry = master_data[pdb_id]
                state = entry["state"]
                chi_data = entry["data"]

                if resnum not in chi_data:
                    continue

                d = chi_data[resnum]
                writer.writerow([
                    resnum,
                    d["resname"],
                    pdb_id,
                    state,
                    d["chi1"] if d["chi1"] is not None else "",
                    d["chi2"] if d["chi2"] is not None else "",
                    d["chi3"] if d["chi3"] is not None else "",
                    d["chi4"] if d["chi4"] is not None else "",
                    classify_rotamer_chi1(d["chi1"]) or "",
                ])

    # ── Summary statistics ───────────────────────────────────────────────────
    summary_file = os.path.join(results_dir, "extraction_summary.txt")
    print(f"Writing summary: {summary_file}")
    with open(summary_file, "w") as f:
        f.write("Chi Angle Extraction Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Total structures processed: {len(master_data)}\n")
        f.write(f"  Inactive: {sum(1 for v in master_data.values() if v['state'] == 'inactive')}\n")
        f.write(f"  Active:   {sum(1 for v in master_data.values() if v['state'] == 'active')}\n")
        f.write(f"Total unique residue positions: {len(all_residues)}\n")
        f.write(f"Residue range: {min(all_residues)} - {max(all_residues)}\n\n")

        f.write("Per-structure counts:\n")
        f.write(f"{'PDB':>6s}  {'State':>8s}  {'Residues':>8s}  {'χ1':>5s}  {'χ2':>5s}  {'χ3':>5s}  {'χ4':>5s}\n")
        f.write("-" * 60 + "\n")
        for pdb_id in sorted(master_data.keys()):
            entry = master_data[pdb_id]
            chi_data = entry["data"]
            counts = {f"chi{i}": sum(1 for v in chi_data.values() if v[f"chi{i}"] is not None) for i in range(1, 5)}
            f.write(f"{pdb_id:>6s}  {entry['state']:>8s}  {len(chi_data):>8d}  "
                    f"{counts['chi1']:>5d}  {counts['chi2']:>5d}  {counts['chi3']:>5d}  {counts['chi4']:>5d}\n")

    print(f"\nDone! Results in {results_dir}/")
    print(f"  - {len(master_data)} per-structure CSV files")
    print(f"  - all_chi_angles_combined.csv (master)")
    print(f"  - extraction_summary.txt")


if __name__ == "__main__":
    main()

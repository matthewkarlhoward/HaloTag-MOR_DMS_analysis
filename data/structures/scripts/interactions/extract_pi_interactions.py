#!/usr/bin/env python3
"""
Extract pi interactions for all residues across all MOR structures.

Interaction types:
  1. Pi-pi stacking: aromatic ring centroid-centroid ≤5.5Å
     - Classified as parallel (dihedral <30°), t_shaped (60–90°), or other
  2. Pi-cation: aromatic centroid to cationic atom (ARG/LYS) ≤6.0Å
  3. Pi-ligand: aromatic centroid to any ligand heavy atom ≤5.5Å
     (only in structures that contain a bound ligand)

For each structure:
  - Parse verified single receptor chain
  - Find all intra-chain pi interactions
  - Map to human residue numbering

Output: combined master CSVs + per-residue counts + summary.
"""

import os
import csv
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

# Drug names associated with PDB IDs (for pi-ligand output labeling)
DRUG_NAMES = {
    "4dkl": "beta-FNA",
    "7ul4": "alvimopan",
    "9PXU": "JDTic",
    "9bjk": "naltrexone",
    "9MQI": "naloxone",
    "5c1m": "BU72",
    "6dde": "DAMGO",
    "6ddf": "DAMGO",
    "7sbf": "morphine",
    "7scg": "DAMGO",
    "7t2g": "fentanyl",
    "7t2h": "fentanyl",
    "7u2k": "lofentanil",
    "7u2l": "carfentanil",
    "8e0g": "endomorphin",
    "8ef5": "DAMGO",
    "8ef6": "leu-enkephalin",
    "8efb": "met-enkephalin",
    "8efl": "endomorphin-1",
    "8efo": "endorphin",
    "8efq": "endomorphin-2",
    "8f7q": "oxymorphone",
    "8f7r": "hydromorphone",
    "8k9k": "buprenorphine",
    "8k9l": "buprenorphine",
    "8Y72": "MP",
    "8Y73": "MP",
    "9bqj": "DAMGO",
    "9PY2": "fentanyl",
    "9PY3": "morphine",
    "9PY4": "nalbuphine",
    "9WST": "fentanyl",
    "9WSW": "morphine",
}

# ── Aromatic ring atom definitions ───────────────────────────────────────────
# Ring atoms used for centroid and normal-vector calculations
AROMATIC_RING_ATOMS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],  # 6-membered ring
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
}

AROMATIC_RESIDUES = set(AROMATIC_RING_ATOMS.keys())

# ── Cationic atom definitions ─────────────────────────────────────────────────
CATIONIC_ATOMS = {
    "ARG": ["NH1", "NH2", "NE"],
    "LYS": ["NZ"],
}

CATIONIC_RESIDUES = set(CATIONIC_ATOMS.keys())


def get_ring_centroid(residue, resname):
    """
    Compute the centroid of the aromatic ring for a given residue.
    Returns numpy array (3,) or None if atoms are missing.
    """
    ring_atom_names = AROMATIC_RING_ATOMS.get(resname, [])
    coords = []
    for aname in ring_atom_names:
        if aname in residue:
            coords.append(residue[aname].get_vector().get_array())
    if len(coords) < 3:
        return None
    return np.mean(coords, axis=0)


def get_ring_normal(residue, resname):
    """
    Compute the normal vector of the aromatic ring plane using the first
    three available ring atoms (cross product of two edge vectors).
    Returns normalised numpy array (3,) or None if atoms are missing.
    """
    ring_atom_names = AROMATIC_RING_ATOMS.get(resname, [])
    coords = []
    for aname in ring_atom_names:
        if aname in residue:
            coords.append(residue[aname].get_vector().get_array())
        if len(coords) == 3:
            break
    if len(coords) < 3:
        return None
    v1 = coords[1] - coords[0]
    v2 = coords[2] - coords[0]
    normal = np.cross(v1, v2)
    norm = np.linalg.norm(normal)
    if norm == 0:
        return None
    return normal / norm


def angle_between_normals(n1, n2):
    """
    Return angle (degrees) between two normal vectors.
    Takes the acute angle (0–90°) between the two planes.
    """
    cos_angle = np.clip(np.dot(n1, n2), -1.0, 1.0)
    angle = np.degrees(np.arccos(abs(cos_angle)))
    # Ensure 0–90 range (planes don't have an orientation)
    if angle > 90:
        angle = 180 - angle
    return angle


def classify_stacking(angle_deg):
    """
    Classify pi-pi stacking geometry from inter-plane angle.
      parallel:  angle < 30°
      t_shaped:  angle 60–90°
      other:     30–60°
    """
    if angle_deg < 30:
        return "parallel"
    elif angle_deg >= 60:
        return "t_shaped"
    else:
        return "other"


def extract_residues(structure, chain_id, species):
    """
    Extract standard amino acid residues with human numbering.
    Returns: list of (human_resnum, residue_object, resname)
    """
    offset = 2 if species == "mouse" else 0
    model = structure[0]
    chain = model[chain_id]
    residues = []

    for residue in chain:
        if not is_aa(residue, standard=True):
            continue
        het_flag = residue.get_id()[0]
        if het_flag != " ":
            continue
        auth_resnum = residue.get_id()[1]
        human_resnum = auth_resnum + offset
        resname = residue.get_resname().strip()
        residues.append((human_resnum, residue, resname))

    return residues


def extract_ligand_atoms(structure, chain_id):
    """
    Extract all HETATM heavy atoms NOT in the receptor chain.
    Returns: list of numpy arrays (3,) of coordinates.
    Only returns atoms from HETATM records (het_flag != " ") and
    excludes water (HOH/WAT).
    """
    model = structure[0]
    ligand_atoms = []
    for chain in model:
        for residue in chain:
            het_flag = residue.get_id()[0]
            # HETATM records start with "H_"
            if not het_flag.startswith("H_"):
                continue
            resname = residue.get_resname().strip()
            if resname in ("HOH", "WAT"):
                continue
            for atom in residue:
                # Skip hydrogens
                if atom.element == "H":
                    continue
                ligand_atoms.append(atom.get_vector().get_array())
    return ligand_atoms


def find_pi_pi(residues, distance_cutoff=5.5):
    """
    Find pi-pi stacking interactions between aromatic residue pairs.
    Criteria:
      - Both residues have defined aromatic rings
      - Centroid-centroid distance <= distance_cutoff
      - Separated by >= 2 residues in sequence

    Returns: list of dicts
    """
    pi_pi = []

    # Filter to aromatic residues with valid centroids
    aromatic = []
    for resnum, res, resname in residues:
        if resname not in AROMATIC_RESIDUES:
            continue
        centroid = get_ring_centroid(res, resname)
        normal = get_ring_normal(res, resname)
        if centroid is None or normal is None:
            continue
        aromatic.append((resnum, res, resname, centroid, normal))

    for i, (resnum_i, res_i, resname_i, centroid_i, normal_i) in enumerate(aromatic):
        for j, (resnum_j, res_j, resname_j, centroid_j, normal_j) in enumerate(aromatic):
            if j <= i:
                continue
            if abs(resnum_i - resnum_j) < 2:
                continue

            dist = float(np.linalg.norm(centroid_i - centroid_j))
            if dist > distance_cutoff:
                continue

            angle = angle_between_normals(normal_i, normal_j)
            stacking_type = classify_stacking(angle)

            pi_pi.append({
                "res1_num": resnum_i,
                "res1_name": resname_i,
                "res2_num": resnum_j,
                "res2_name": resname_j,
                "centroid_distance": round(dist, 3),
                "stacking_type": stacking_type,
            })

    return pi_pi


def find_pi_cation(residues, distance_cutoff=6.0):
    """
    Find pi-cation interactions between aromatic and cationic residues.
    Criteria:
      - Aromatic centroid to cationic heavy atom distance <= distance_cutoff
      - Report the closest cationic atom

    Returns: list of dicts
    """
    pi_cation = []

    aromatic = []
    for resnum, res, resname in residues:
        if resname not in AROMATIC_RESIDUES:
            continue
        centroid = get_ring_centroid(res, resname)
        if centroid is None:
            continue
        aromatic.append((resnum, res, resname, centroid))

    cationic = [(rn, res, rname) for rn, res, rname in residues if rname in CATIONIC_RESIDUES]

    for resnum_a, res_a, resname_a, centroid_a in aromatic:
        for resnum_c, res_c, resname_c in cationic:
            if resnum_a == resnum_c:
                continue

            cat_atom_names = CATIONIC_ATOMS[resname_c]
            min_dist = float("inf")
            best_atom = None
            for aname in cat_atom_names:
                if aname not in res_c:
                    continue
                atom_coord = res_c[aname].get_vector().get_array()
                dist = float(np.linalg.norm(centroid_a - atom_coord))
                if dist < min_dist:
                    min_dist = dist
                    best_atom = aname

            if min_dist <= distance_cutoff and best_atom is not None:
                pi_cation.append({
                    "aromatic_res_num": resnum_a,
                    "aromatic_res_name": resname_a,
                    "cation_res_num": resnum_c,
                    "cation_res_name": resname_c,
                    "cation_atom": best_atom,
                    "distance": round(min_dist, 3),
                })

    return pi_cation


def find_pi_ligand(residues, ligand_atoms, distance_cutoff=5.5):
    """
    Find pi-ligand interactions between aromatic receptor residues and
    ligand heavy atoms.
    Criteria:
      - Aromatic centroid to any ligand heavy atom distance <= distance_cutoff
      - Report minimum distance across all ligand atoms

    Returns: list of dicts
    """
    if not ligand_atoms:
        return []

    pi_ligand = []
    ligand_coords = np.array(ligand_atoms)  # shape (N, 3)

    for resnum, res, resname in residues:
        if resname not in AROMATIC_RESIDUES:
            continue
        centroid = get_ring_centroid(res, resname)
        if centroid is None:
            continue

        # Distances from centroid to all ligand atoms
        diffs = ligand_coords - centroid[np.newaxis, :]
        dists = np.linalg.norm(diffs, axis=1)
        min_dist = float(np.min(dists))

        if min_dist <= distance_cutoff:
            pi_ligand.append({
                "aromatic_res_num": resnum,
                "aromatic_res_name": resname,
                "min_ligand_distance": round(min_dist, 3),
            })

    return pi_ligand


def main():
    results_dir = str(REPO_ROOT / "structures" / "processed" / "pi_interactions")
    os.makedirs(results_dir, exist_ok=True)

    parser = MMCIFParser(QUIET=True)

    # Combine all structures
    all_structures = {}
    for pdb_id, info in INACTIVE_STRUCTURES.items():
        all_structures[pdb_id] = {**info, "state": "inactive"}
    for pdb_id, info in ACTIVE_STRUCTURES.items():
        all_structures[pdb_id] = {**info, "state": "active"}

    print(f"Extracting pi interactions from {len(all_structures)} structures...")
    print(f"  Inactive: {len(INACTIVE_STRUCTURES)}")
    print(f"  Active:   {len(ACTIVE_STRUCTURES)}")
    print()

    # Master interaction lists
    all_pi_pi = []
    all_pi_cation = []
    all_pi_ligand = []

    structures_processed = 0
    structures_skipped = 0

    for pdb_id, info in sorted(all_structures.items()):
        chain_id = info["receptor_chain"]
        species = info["species"]
        state = info["state"]
        drug_name = DRUG_NAMES.get(pdb_id, "unknown")

        # Find CIF file
        cif_file = None
        for pattern in [f"{pdb_id}.cif", f"{pdb_id.lower()}.cif", f"{pdb_id.upper()}.cif"]:
            candidate = os.path.join(str(CIF_DIR), pattern)
            if os.path.exists(candidate):
                cif_file = candidate
                break

        if cif_file is None:
            print(f"  WARNING: CIF file not found for {pdb_id}, skipping")
            structures_skipped += 1
            continue

        try:
            structure = parser.get_structure(pdb_id, cif_file)
            residues = extract_residues(structure, chain_id, species)
            ligand_atoms = extract_ligand_atoms(structure, chain_id)

            print(
                f"  {pdb_id} ({state:8s}, chain {chain_id}, "
                f"{len(residues)} res, {len(ligand_atoms)} ligand atoms)...",
                end="",
                flush=True,
            )

            # Find interactions
            pi_pi = find_pi_pi(residues)
            pi_cation = find_pi_cation(residues)
            pi_ligand = find_pi_ligand(residues, ligand_atoms)

            print(
                f" {len(pi_pi)} pi-pi, {len(pi_cation)} pi-cation, "
                f"{len(pi_ligand)} pi-ligand"
            )

            # Annotate with pdb_id and state
            for r in pi_pi:
                r["pdb_id"] = pdb_id
                r["state"] = state
            for r in pi_cation:
                r["pdb_id"] = pdb_id
                r["state"] = state
            for r in pi_ligand:
                r["pdb_id"] = pdb_id
                r["state"] = state
                r["drug_name"] = drug_name

            all_pi_pi.extend(pi_pi)
            all_pi_cation.extend(pi_cation)
            all_pi_ligand.extend(pi_ligand)

            structures_processed += 1

        except Exception as e:
            print(f" ERROR: {e}")
            structures_skipped += 1
            continue

    # ── Write combined master CSVs ───────────────────────────────────────────
    print(f"\nWriting combined CSVs...")

    pi_pi_file = os.path.join(results_dir, "all_pi_pi.csv")
    with open(pi_pi_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "pdb_id", "state",
            "res1_num", "res1_name",
            "res2_num", "res2_name",
            "centroid_distance", "stacking_type",
        ])
        writer.writeheader()
        for r in all_pi_pi:
            writer.writerow({k: r[k] for k in writer.fieldnames})
    print(f"  {pi_pi_file}: {len(all_pi_pi)} total pi-pi interactions")

    pi_cation_file = os.path.join(results_dir, "all_pi_cation.csv")
    with open(pi_cation_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "pdb_id", "state",
            "aromatic_res_num", "aromatic_res_name",
            "cation_res_num", "cation_res_name",
            "cation_atom", "distance",
        ])
        writer.writeheader()
        for r in all_pi_cation:
            writer.writerow({k: r[k] for k in writer.fieldnames})
    print(f"  {pi_cation_file}: {len(all_pi_cation)} total pi-cation interactions")

    pi_ligand_file = os.path.join(results_dir, "all_pi_ligand.csv")
    with open(pi_ligand_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "pdb_id", "state", "drug_name",
            "aromatic_res_num", "aromatic_res_name",
            "min_ligand_distance",
        ])
        writer.writeheader()
        for r in all_pi_ligand:
            writer.writerow({k: r[k] for k in writer.fieldnames})
    print(f"  {pi_ligand_file}: {len(all_pi_ligand)} total pi-ligand interactions")

    # ── Per-residue interaction count summary ────────────────────────────────
    print("\nGenerating per-residue pi interaction counts...")

    # residue_counts[pdb_id][human_resnum] = {pi_pi: N, pi_cation: N, pi_ligand: N}
    residue_counts = defaultdict(lambda: defaultdict(lambda: {
        "pi_pi_count": 0,
        "pi_cation_count": 0,
        "pi_ligand_count": 0,
    }))

    for r in all_pi_pi:
        pdb_id = r["pdb_id"]
        residue_counts[pdb_id][r["res1_num"]]["pi_pi_count"] += 1
        residue_counts[pdb_id][r["res2_num"]]["pi_pi_count"] += 1

    for r in all_pi_cation:
        pdb_id = r["pdb_id"]
        residue_counts[pdb_id][r["aromatic_res_num"]]["pi_cation_count"] += 1
        residue_counts[pdb_id][r["cation_res_num"]]["pi_cation_count"] += 1

    for r in all_pi_ligand:
        pdb_id = r["pdb_id"]
        residue_counts[pdb_id][r["aromatic_res_num"]]["pi_ligand_count"] += 1

    count_file = os.path.join(results_dir, "per_residue_pi_counts.csv")
    with open(count_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "pdb_id", "state", "residue_number_human",
            "pi_pi_count", "pi_cation_count", "pi_ligand_count",
        ])
        for pdb_id in sorted(residue_counts.keys()):
            state = all_structures[pdb_id]["state"]
            for resnum in sorted(residue_counts[pdb_id].keys()):
                c = residue_counts[pdb_id][resnum]
                writer.writerow([
                    pdb_id, state, resnum,
                    c["pi_pi_count"], c["pi_cation_count"], c["pi_ligand_count"],
                ])
    print(f"  {count_file}")

    # ── Summary ──────────────────────────────────────────────────────────────
    n_parallel = sum(1 for r in all_pi_pi if r["stacking_type"] == "parallel")
    n_t_shaped = sum(1 for r in all_pi_pi if r["stacking_type"] == "t_shaped")
    n_other    = sum(1 for r in all_pi_pi if r["stacking_type"] == "other")

    summary_file = os.path.join(results_dir, "extraction_summary.txt")
    with open(summary_file, "w") as f:
        f.write("Pi Interaction Extraction Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Structures processed: {structures_processed}\n")
        f.write(f"Structures skipped:   {structures_skipped}\n")
        f.write(f"  Inactive defined:   {len(INACTIVE_STRUCTURES)}\n")
        f.write(f"  Active defined:     {len(ACTIVE_STRUCTURES)}\n\n")

        f.write("Total pi interactions found:\n")
        f.write(f"  Pi-pi stacking:      {len(all_pi_pi):>8d}\n")
        f.write(f"    parallel:          {n_parallel:>8d}\n")
        f.write(f"    t-shaped:          {n_t_shaped:>8d}\n")
        f.write(f"    other:             {n_other:>8d}\n")
        f.write(f"  Pi-cation:           {len(all_pi_cation):>8d}\n")
        f.write(f"  Pi-ligand:           {len(all_pi_ligand):>8d}\n\n")

        f.write("Cutoffs used:\n")
        f.write("  Pi-pi centroid-centroid:  5.5 Å\n")
        f.write("  Pi-cation centroid-atom:  6.0 Å\n")
        f.write("  Pi-ligand centroid-atom:  5.5 Å\n\n")

        f.write("Stacking classification (inter-plane angle):\n")
        f.write("  parallel:  < 30°\n")
        f.write("  t_shaped:  60–90°\n")
        f.write("  other:     30–60°\n\n")

        f.write("Aromatic residues detected: PHE, TYR, TRP, HIS\n")
        f.write("Cationic residues detected: ARG, LYS\n")

    print(f"  {summary_file}")
    print(f"\nDone!")


if __name__ == "__main__":
    main()

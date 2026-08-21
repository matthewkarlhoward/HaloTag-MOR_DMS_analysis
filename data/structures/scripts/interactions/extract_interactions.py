#!/usr/bin/env python3
"""
Extract non-bonded interactions for all residues across all MOR structures.

Interaction types:
  1. Hydrogen bonds: donor-acceptor distance ≤3.5Å, D-H...A angle ≥120°
     - Backbone-backbone, backbone-sidechain, sidechain-sidechain
  2. Salt bridges: charged atom pairs (Asp/Glu COO⁻ ↔ Arg/Lys/His NH₃⁺) ≤4.0Å
  3. Hydrophobic contacts: carbon-carbon distance ≤4.5Å between non-polar sidechains
     (Ala, Val, Leu, Ile, Pro, Phe, Trp, Met)

For each structure:
  - Parse verified single receptor chain
  - Find all intra-chain interactions
  - Map to human residue numbering

Output: per-structure CSV + combined master CSV.
"""

import os
import csv
import numpy as np
from collections import defaultdict
from itertools import combinations

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


# ── Atom classification ──────────────────────────────────────────────────────

# Hydrogen bond donors and acceptors (heavy atoms only; H positions often missing)
# We use a distance-only criterion with donor/acceptor atom types
HBOND_DONORS = {
    # Backbone
    "N",
    # Sidechains
    "ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ",  # Nitrogen donors
    "OG", "OG1", "OH",  # Oxygen donors (Ser, Thr, Tyr)
}

HBOND_ACCEPTORS = {
    # Backbone
    "O",
    # Sidechains
    "OD1", "OD2", "OE1", "OE2",  # Asp, Glu
    "OG", "OG1", "OH",            # Ser, Thr, Tyr (also donors)
    "ND1", "NE2",                  # His (can be acceptor)
    "SD",                          # Met sulfur
    "SG",                          # Cys sulfur
}

# Salt bridge atoms
POSITIVE_ATOMS = {
    "ARG": ["NH1", "NH2", "NE"],
    "LYS": ["NZ"],
    "HIS": ["ND1", "NE2"],
}

NEGATIVE_ATOMS = {
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
}

# Hydrophobic residues and their sidechain carbon atoms
HYDROPHOBIC_RESIDUES = {"ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "TRP", "MET"}
BACKBONE_ATOMS = {"N", "CA", "C", "O"}


def get_distance(atom1, atom2):
    """Euclidean distance between two atoms."""
    return np.linalg.norm(atom1.get_vector().get_array() - atom2.get_vector().get_array())


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


def find_hbonds(residues, distance_cutoff=3.5):
    """
    Find hydrogen bonds between residue pairs.
    Uses distance-only criterion (heavy atom donor-acceptor ≤ cutoff).
    Only considers pairs separated by ≥ 2 residues in sequence.

    Returns: list of dicts with interaction details
    """
    hbonds = []

    for i, (resnum_i, res_i, resname_i) in enumerate(residues):
        for j, (resnum_j, res_j, resname_j) in enumerate(residues):
            if j <= i:
                continue
            if abs(resnum_i - resnum_j) < 2:
                continue  # Skip immediate neighbors

            # Check all donor-acceptor pairs
            for atom_d in res_i:
                if atom_d.get_name() not in HBOND_DONORS:
                    continue
                for atom_a in res_j:
                    if atom_a.get_name() not in HBOND_ACCEPTORS:
                        continue
                    dist = get_distance(atom_d, atom_a)
                    if dist <= distance_cutoff:
                        # Classify as bb-bb, bb-sc, sc-sc
                        d_is_bb = atom_d.get_name() in BACKBONE_ATOMS
                        a_is_bb = atom_a.get_name() in BACKBONE_ATOMS
                        if d_is_bb and a_is_bb:
                            hb_type = "bb-bb"
                        elif d_is_bb or a_is_bb:
                            hb_type = "bb-sc"
                        else:
                            hb_type = "sc-sc"

                        hbonds.append({
                            "res1_num": resnum_i,
                            "res1_name": resname_i,
                            "res1_atom": atom_d.get_name(),
                            "res2_num": resnum_j,
                            "res2_name": resname_j,
                            "res2_atom": atom_a.get_name(),
                            "distance": round(dist, 3),
                            "hb_type": hb_type,
                        })

            # Also check reverse direction (j as donor, i as acceptor)
            for atom_d in res_j:
                if atom_d.get_name() not in HBOND_DONORS:
                    continue
                for atom_a in res_i:
                    if atom_a.get_name() not in HBOND_ACCEPTORS:
                        continue
                    dist = get_distance(atom_d, atom_a)
                    if dist <= distance_cutoff:
                        d_is_bb = atom_d.get_name() in BACKBONE_ATOMS
                        a_is_bb = atom_a.get_name() in BACKBONE_ATOMS
                        if d_is_bb and a_is_bb:
                            hb_type = "bb-bb"
                        elif d_is_bb or a_is_bb:
                            hb_type = "bb-sc"
                        else:
                            hb_type = "sc-sc"

                        hbonds.append({
                            "res1_num": resnum_j,
                            "res1_name": resname_j,
                            "res1_atom": atom_d.get_name(),
                            "res2_num": resnum_i,
                            "res2_name": resname_i,
                            "res2_atom": atom_a.get_name(),
                            "distance": round(dist, 3),
                            "hb_type": hb_type,
                        })

    return hbonds


def find_salt_bridges(residues, distance_cutoff=4.0):
    """
    Find salt bridges: charged atom pairs within cutoff.
    Returns: list of dicts
    """
    salt_bridges = []

    # Collect positive and negative residues
    pos_residues = [(rn, res, rname) for rn, res, rname in residues if rname in POSITIVE_ATOMS]
    neg_residues = [(rn, res, rname) for rn, res, rname in residues if rname in NEGATIVE_ATOMS]

    for resnum_p, res_p, resname_p in pos_residues:
        pos_atoms_names = POSITIVE_ATOMS[resname_p]
        for resnum_n, res_n, resname_n in neg_residues:
            neg_atoms_names = NEGATIVE_ATOMS[resname_n]

            min_dist = float("inf")
            best_pair = None
            for pa_name in pos_atoms_names:
                if pa_name not in res_p:
                    continue
                pa = res_p[pa_name]
                for na_name in neg_atoms_names:
                    if na_name not in res_n:
                        continue
                    na = res_n[na_name]
                    dist = get_distance(pa, na)
                    if dist < min_dist:
                        min_dist = dist
                        best_pair = (pa_name, na_name)

            if min_dist <= distance_cutoff and best_pair:
                salt_bridges.append({
                    "res1_num": resnum_p,
                    "res1_name": resname_p,
                    "res1_atom": best_pair[0],
                    "res2_num": resnum_n,
                    "res2_name": resname_n,
                    "res2_atom": best_pair[1],
                    "distance": round(min_dist, 3),
                })

    return salt_bridges


def find_hydrophobic_contacts(residues, distance_cutoff=4.5):
    """
    Find hydrophobic packing contacts: sidechain carbon-carbon distances
    between hydrophobic residues within cutoff.
    Reports minimum C-C distance per residue pair.
    Only considers pairs separated by ≥ 3 residues.
    """
    contacts = []

    # Filter to hydrophobic residues
    hydro_residues = [(rn, res, rname) for rn, res, rname in residues if rname in HYDROPHOBIC_RESIDUES]

    for i, (resnum_i, res_i, resname_i) in enumerate(hydro_residues):
        for j, (resnum_j, res_j, resname_j) in enumerate(hydro_residues):
            if j <= i:
                continue
            if abs(resnum_i - resnum_j) < 3:
                continue

            # Get sidechain carbon atoms
            carbons_i = [a for a in res_i if a.element == "C" and a.get_name() not in BACKBONE_ATOMS]
            carbons_j = [a for a in res_j if a.element == "C" and a.get_name() not in BACKBONE_ATOMS]

            if not carbons_i or not carbons_j:
                continue

            min_dist = float("inf")
            best_pair = None
            for ci in carbons_i:
                for cj in carbons_j:
                    dist = get_distance(ci, cj)
                    if dist < min_dist:
                        min_dist = dist
                        best_pair = (ci.get_name(), cj.get_name())

            if min_dist <= distance_cutoff and best_pair:
                contacts.append({
                    "res1_num": resnum_i,
                    "res1_name": resname_i,
                    "res1_atom": best_pair[0],
                    "res2_num": resnum_j,
                    "res2_name": resname_j,
                    "res2_atom": best_pair[1],
                    "distance": round(min_dist, 3),
                })

    return contacts


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cif_dir = str(CIF_DIR)
    results_dir = str(REPO_ROOT / "structures" / "processed" / "interactions")
    os.makedirs(results_dir, exist_ok=True)

    parser = MMCIFParser(QUIET=True)

    # Combine all structures
    all_structures = {}
    for pdb_id, info in INACTIVE_STRUCTURES.items():
        all_structures[pdb_id] = {**info, "state": "inactive"}
    for pdb_id, info in ACTIVE_STRUCTURES.items():
        all_structures[pdb_id] = {**info, "state": "active"}

    print(f"Extracting interactions from {len(all_structures)} structures...")
    print(f"  Inactive: {len(INACTIVE_STRUCTURES)}")
    print(f"  Active:   {len(ACTIVE_STRUCTURES)}")
    print()

    # Master interaction lists
    all_hbonds = []
    all_salt_bridges = []
    all_hydrophobic = []

    for pdb_id, info in sorted(all_structures.items()):
        chain_id = info["receptor_chain"]
        species = info["species"]
        state = info["state"]

        # Find CIF file
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
            residues = extract_residues(structure, chain_id, species)

            print(f"  {pdb_id} ({state:8s}, chain {chain_id}, {len(residues)} res)...", end="", flush=True)

            # Find interactions
            hbonds = find_hbonds(residues)
            salt_bridges = find_salt_bridges(residues)
            hydrophobic = find_hydrophobic_contacts(residues)

            print(f" {len(hbonds)} H-bonds, {len(salt_bridges)} salt bridges, {len(hydrophobic)} hydrophobic")

            # Add pdb_id and state to each interaction
            for h in hbonds:
                h["pdb_id"] = pdb_id
                h["state"] = state
            for s in salt_bridges:
                s["pdb_id"] = pdb_id
                s["state"] = state
            for c in hydrophobic:
                c["pdb_id"] = pdb_id
                c["state"] = state

            all_hbonds.extend(hbonds)
            all_salt_bridges.extend(salt_bridges)
            all_hydrophobic.extend(hydrophobic)

            # Write per-structure CSVs
            per_struct_dir = os.path.join(results_dir, "per_structure")
            os.makedirs(per_struct_dir, exist_ok=True)

            # H-bonds
            hb_file = os.path.join(per_struct_dir, f"{pdb_id}_hbonds.csv")
            with open(hb_file, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=[
                    "res1_num", "res1_name", "res1_atom",
                    "res2_num", "res2_name", "res2_atom",
                    "distance", "hb_type"
                ])
                writer.writeheader()
                for h in hbonds:
                    writer.writerow({k: h[k] for k in writer.fieldnames})

            # Salt bridges
            sb_file = os.path.join(per_struct_dir, f"{pdb_id}_salt_bridges.csv")
            with open(sb_file, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=[
                    "res1_num", "res1_name", "res1_atom",
                    "res2_num", "res2_name", "res2_atom",
                    "distance"
                ])
                writer.writeheader()
                for s in salt_bridges:
                    writer.writerow({k: s[k] for k in writer.fieldnames})

            # Hydrophobic contacts
            hp_file = os.path.join(per_struct_dir, f"{pdb_id}_hydrophobic.csv")
            with open(hp_file, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=[
                    "res1_num", "res1_name", "res1_atom",
                    "res2_num", "res2_name", "res2_atom",
                    "distance"
                ])
                writer.writeheader()
                for c in hydrophobic:
                    writer.writerow({k: c[k] for k in writer.fieldnames})

        except Exception as e:
            print(f" ERROR: {e}")
            continue

    # ── Write combined master CSVs ───────────────────────────────────────────
    print(f"\nWriting combined CSVs...")

    hb_master = os.path.join(results_dir, "all_hbonds.csv")
    with open(hb_master, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "pdb_id", "state",
            "res1_num", "res1_name", "res1_atom",
            "res2_num", "res2_name", "res2_atom",
            "distance", "hb_type"
        ])
        writer.writeheader()
        for h in all_hbonds:
            writer.writerow({k: h[k] for k in writer.fieldnames})
    print(f"  {hb_master}: {len(all_hbonds)} total H-bonds")

    sb_master = os.path.join(results_dir, "all_salt_bridges.csv")
    with open(sb_master, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "pdb_id", "state",
            "res1_num", "res1_name", "res1_atom",
            "res2_num", "res2_name", "res2_atom",
            "distance"
        ])
        writer.writeheader()
        for s in all_salt_bridges:
            writer.writerow({k: s[k] for k in writer.fieldnames})
    print(f"  {sb_master}: {len(all_salt_bridges)} total salt bridges")

    hp_master = os.path.join(results_dir, "all_hydrophobic_contacts.csv")
    with open(hp_master, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "pdb_id", "state",
            "res1_num", "res1_name", "res1_atom",
            "res2_num", "res2_name", "res2_atom",
            "distance"
        ])
        writer.writeheader()
        for c in all_hydrophobic:
            writer.writerow({k: c[k] for k in writer.fieldnames})
    print(f"  {hp_master}: {len(all_hydrophobic)} total hydrophobic contacts")

    # ── Per-residue interaction count summary ────────────────────────────────
    # For each residue in each structure, count how many of each interaction type it participates in
    print("\nGenerating per-residue interaction counts...")

    # Collect per-residue counts per structure
    residue_counts = defaultdict(lambda: defaultdict(lambda: {
        "hbond_total": 0, "hbond_bb_bb": 0, "hbond_bb_sc": 0, "hbond_sc_sc": 0,
        "salt_bridge": 0, "hydrophobic": 0
    }))

    for h in all_hbonds:
        pdb_id = h["pdb_id"]
        residue_counts[pdb_id][h["res1_num"]]["hbond_total"] += 1
        residue_counts[pdb_id][h["res2_num"]]["hbond_total"] += 1
        residue_counts[pdb_id][h["res1_num"]][f"hbond_{h['hb_type'].replace('-', '_')}"] += 1
        residue_counts[pdb_id][h["res2_num"]][f"hbond_{h['hb_type'].replace('-', '_')}"] += 1

    for s in all_salt_bridges:
        pdb_id = s["pdb_id"]
        residue_counts[pdb_id][s["res1_num"]]["salt_bridge"] += 1
        residue_counts[pdb_id][s["res2_num"]]["salt_bridge"] += 1

    for c in all_hydrophobic:
        pdb_id = c["pdb_id"]
        residue_counts[pdb_id][c["res1_num"]]["hydrophobic"] += 1
        residue_counts[pdb_id][c["res2_num"]]["hydrophobic"] += 1

    # Write per-residue counts (combined across all structures)
    count_file = os.path.join(results_dir, "per_residue_interaction_counts.csv")
    with open(count_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "pdb_id", "state", "residue_number_human",
            "hbond_total", "hbond_bb_bb", "hbond_bb_sc", "hbond_sc_sc",
            "salt_bridge", "hydrophobic"
        ])
        for pdb_id in sorted(residue_counts.keys()):
            state = all_structures[pdb_id]["state"]
            for resnum in sorted(residue_counts[pdb_id].keys()):
                c = residue_counts[pdb_id][resnum]
                writer.writerow([
                    pdb_id, state, resnum,
                    c["hbond_total"], c["hbond_bb_bb"], c["hbond_bb_sc"], c["hbond_sc_sc"],
                    c["salt_bridge"], c["hydrophobic"]
                ])
    print(f"  {count_file}")

    # ── Summary ──────────────────────────────────────────────────────────────
    summary_file = os.path.join(results_dir, "extraction_summary.txt")
    with open(summary_file, "w") as f:
        f.write("Interaction Extraction Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Structures processed: {len(all_structures)}\n")
        f.write(f"  Inactive: {len(INACTIVE_STRUCTURES)}\n")
        f.write(f"  Active:   {len(ACTIVE_STRUCTURES)}\n\n")

        f.write(f"Total interactions found:\n")
        f.write(f"  Hydrogen bonds:      {len(all_hbonds):>8d}\n")
        n_bb = sum(1 for h in all_hbonds if h["hb_type"] == "bb-bb")
        n_bs = sum(1 for h in all_hbonds if h["hb_type"] == "bb-sc")
        n_ss = sum(1 for h in all_hbonds if h["hb_type"] == "sc-sc")
        f.write(f"    backbone-backbone: {n_bb:>8d}\n")
        f.write(f"    backbone-sidechain:{n_bs:>8d}\n")
        f.write(f"    sidechain-sidechain:{n_ss:>7d}\n")
        f.write(f"  Salt bridges:        {len(all_salt_bridges):>8d}\n")
        f.write(f"  Hydrophobic contacts:{len(all_hydrophobic):>8d}\n\n")

        f.write("Cutoffs used:\n")
        f.write("  H-bond donor-acceptor:    3.5 Å\n")
        f.write("  Salt bridge:              4.0 Å\n")
        f.write("  Hydrophobic C-C:          4.5 Å\n")
        f.write("  H-bond min seq separation:  2 residues\n")
        f.write("  Hydrophobic min seq sep:    3 residues\n")

    print(f"  {summary_file}")
    print(f"\nDone!")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Calculate per-residue CA (alpha carbon) displacement between active and inactive
MOR structures after structural superposition.

For each active-state structure:
  1. Superimpose onto each inactive-state structure (CA-based alignment)
  2. Calculate per-residue CA displacement (Å)
  3. Report distances to each inactive structure + mean across all inactive
Then calculate grand average per-residue displacement across all active structures.

Uses verified single receptor chains and human residue numbering alignment.
"""

import os
import csv
import numpy as np
from collections import defaultdict

from Bio.PDB import MMCIFParser, Superimposer
from Bio.PDB.Polypeptide import is_aa

from pathlib import Path
REPO_ROOT = Path(__file__).resolve().parents[3]
CIF_DIR = REPO_ROOT / "structures" / "raw" / "experimental"

# ── Structure definitions ────────────────────────────────────────────────────
# Classification:
#   inactive: antagonist-bound structures (receptor in inactive conformation)
#   active:   agonist-bound, typically G-protein-coupled
#
# Excluded entirely:
#   8qot  — apo (no orthosteric ligand)
#   9PXV, 9PXX, 9PXY — naloxone + G-protein outliers (intermediate conformation)
#   9PXW  — naloxone + G-protein intermediate
#   9WSV, 9WSX — arrestin-bound (not G-protein)

INACTIVE_STRUCTURES = {
    "4dkl":  {"receptor_chain": "A", "species": "mouse"},   # beta-FNA (irreversible antagonist)
    "7ul4":  {"receptor_chain": "A", "species": "mouse"},   # Naltrexone-derived (antagonist)
    "9PXU":  {"receptor_chain": "R", "species": "human"},   # Naloxone (antagonist, no G-protein)
    "9bjk":  {"receptor_chain": "R", "species": "mouse"},   # A1APU bitopic (naloxone + PAM, inactive)
    "9MQI":  {"receptor_chain": "A", "species": "human"},   # A1BNM (inactive state)
}

ACTIVE_STRUCTURES = {
    "5c1m":  {"receptor_chain": "A", "species": "mouse"},   # BU72
    "6dde":  {"receptor_chain": "R", "species": "mouse"},   # DAMGO
    "6ddf":  {"receptor_chain": "R", "species": "mouse"},   # DAMGO
    "7sbf":  {"receptor_chain": "R", "species": "mouse"},   # PZM21
    "7scg":  {"receptor_chain": "D", "species": "mouse"},   # FH210
    "7t2g":  {"receptor_chain": "R", "species": "mouse"},   # Mitragynine Pseudoindoxyl
    "7t2h":  {"receptor_chain": "D", "species": "mouse"},   # Lofentanil
    "7u2k":  {"receptor_chain": "D", "species": "mouse"},   # C6-guano
    "7u2l":  {"receptor_chain": "D", "species": "mouse"},   # C5-guano
    "8e0g":  {"receptor_chain": "A", "species": "mouse"},   # BU72 analog
    "8ef5":  {"receptor_chain": "R", "species": "human"},   # Fentanyl
    "8ef6":  {"receptor_chain": "R", "species": "human"},   # Morphine
    "8efb":  {"receptor_chain": "R", "species": "human"},   # Oliceridine
    "8efl":  {"receptor_chain": "R", "species": "human"},   # SR17018
    "8efo":  {"receptor_chain": "R", "species": "human"},   # PZM21
    "8efq":  {"receptor_chain": "R", "species": "human"},   # DAMGO
    "8f7q":  {"receptor_chain": "R", "species": "human"},   # Endomorphin-1 analog
    "8f7r":  {"receptor_chain": "R", "species": "human"},   # Endomorphin-2
    "8k9k":  {"receptor_chain": "R", "species": "human"},   # DAMGO
    "8k9l":  {"receptor_chain": "R", "species": "human"},   # DAMGO
    "8Y72":  {"receptor_chain": "R", "species": "human"},   # DAMGO
    "8Y73":  {"receptor_chain": "R", "species": "human"},   # DAMGO
    "9bqj":  {"receptor_chain": "D", "species": "mouse"},   # RO76
    "9PY2":  {"receptor_chain": "R", "species": "human"},   # A1CMV
    "9PY3":  {"receptor_chain": "R", "species": "human"},   # A1CMV
    "9PY4":  {"receptor_chain": "R", "species": "human"},   # A1CMV
    "9WST":  {"receptor_chain": "R", "species": "mouse"},   # DAMGO
    "9WSW":  {"receptor_chain": "R", "species": "mouse"},   # Endomorphin-2
}


def load_ca_coords(cif_file, pdb_id, chain_id, species):
    """
    Load CA atom coordinates from a structure, keyed by human residue number.

    Returns:
        dict: {human_resnum: (x, y, z, residue_name)}
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id, cif_file)
    model = structure[0]

    try:
        chain = model[chain_id]
    except KeyError:
        print(f"  WARNING: chain {chain_id} not found in {pdb_id}")
        return None

    human_offset = 2 if species == "mouse" else 0
    ca_coords = {}

    for residue in chain:
        if not is_aa(residue, standard=True):
            continue
        if 'CA' not in residue:
            continue

        res_num = residue.id[1]
        human_num = res_num + human_offset
        ca = residue['CA']
        ca_coords[human_num] = (
            ca.get_vector().get_array().copy(),
            residue.get_resname()
        )

    return ca_coords


def superimpose_and_calculate(active_coords, inactive_coords):
    """
    Superimpose active structure onto inactive using common CA atoms,
    then calculate per-residue CA displacement.

    Returns:
        dict: {human_resnum: displacement_angstrom}
        float: overall RMSD
        int: number of atoms used for alignment
    """
    # Find common residue positions
    common_positions = sorted(
        set(active_coords.keys()) & set(inactive_coords.keys())
    )

    if len(common_positions) < 50:
        return None, None, 0

    # Build coordinate arrays for superposition
    fixed_atoms = []  # inactive (reference)
    moving_atoms = []  # active (to be moved)

    for pos in common_positions:
        fixed_atoms.append(inactive_coords[pos][0])
        moving_atoms.append(active_coords[pos][0])

    fixed_array = np.array(fixed_atoms)
    moving_array = np.array(moving_atoms)

    # Calculate superposition using SVD (same as Bio.PDB.Superimposer)
    # Center both sets
    fixed_center = fixed_array.mean(axis=0)
    moving_center = moving_array.mean(axis=0)

    fixed_centered = fixed_array - fixed_center
    moving_centered = moving_array - moving_center

    # SVD to find optimal rotation
    H = moving_centered.T @ fixed_centered
    U, S, Vt = np.linalg.svd(H)

    # Handle reflection case
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.eye(3)
    if d < 0:
        sign_matrix[2, 2] = -1

    rotation = Vt.T @ sign_matrix @ U.T
    translation = fixed_center - rotation @ moving_center

    # Apply transformation to all active CA coordinates
    transformed = {}
    for pos, (coord, resname) in active_coords.items():
        new_coord = rotation @ coord + translation
        transformed[pos] = (new_coord, resname)

    # Calculate per-residue displacement at common positions
    displacements = {}
    for pos in common_positions:
        active_ca = transformed[pos][0]
        inactive_ca = inactive_coords[pos][0]
        dist = np.linalg.norm(active_ca - inactive_ca)
        displacements[pos] = dist

    # Overall RMSD
    rmsd = np.sqrt(np.mean([d**2 for d in displacements.values()]))

    return displacements, rmsd, len(common_positions)


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cif_dir = str(CIF_DIR)
    output_dir = str(REPO_ROOT / "structures" / "processed" / "ca_displacements")
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 80)
    print("MOR CA Distance Analysis: Active vs Inactive Structures")
    print("=" * 80)
    print(f"\nCIF directory: {cif_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Active structures:   {len(ACTIVE_STRUCTURES)}")
    print(f"Inactive structures: {len(INACTIVE_STRUCTURES)}")

    # ── Step 1: Load all structures ──────────────────────────────────────────
    print("\n--- Loading structures ---")

    inactive_data = {}
    for pdb_id, sdef in sorted(INACTIVE_STRUCTURES.items(), key=lambda x: x[0].lower()):
        cif_path = os.path.join(cif_dir, f"{pdb_id}.cif")
        if not os.path.exists(cif_path):
            print(f"  WARNING: {pdb_id}.cif not found, skipping")
            continue
        coords = load_ca_coords(cif_path, pdb_id, sdef["receptor_chain"], sdef["species"])
        if coords:
            inactive_data[pdb_id] = coords
            print(f"  {pdb_id}: {len(coords)} CA atoms (inactive, {sdef['species']})")

    active_data = {}
    for pdb_id, sdef in sorted(ACTIVE_STRUCTURES.items(), key=lambda x: x[0].lower()):
        cif_path = os.path.join(cif_dir, f"{pdb_id}.cif")
        if not os.path.exists(cif_path):
            print(f"  WARNING: {pdb_id}.cif not found, skipping")
            continue
        coords = load_ca_coords(cif_path, pdb_id, sdef["receptor_chain"], sdef["species"])
        if coords:
            active_data[pdb_id] = coords
            print(f"  {pdb_id}: {len(coords)} CA atoms (active, {sdef['species']})")

    print(f"\nLoaded: {len(inactive_data)} inactive, {len(active_data)} active")

    if not inactive_data or not active_data:
        print("ERROR: Need at least one active and one inactive structure.")
        return

    inactive_ids = sorted(inactive_data.keys(), key=str.lower)

    # ── Step 2: Per-active-structure analysis ────────────────────────────────
    print("\n--- Calculating CA displacements ---")

    # Accumulator for grand average across all active structures
    # grand_displacements[human_resnum] = list of displacement values
    grand_displacements = defaultdict(list)
    grand_residue_names = {}

    all_active_results = {}

    for active_id in sorted(active_data.keys(), key=str.lower):
        active_coords = active_data[active_id]
        active_species = ACTIVE_STRUCTURES[active_id]["species"]

        print(f"\n  {active_id} ({active_species}):")

        per_inactive = {}  # {inactive_id: {resnum: displacement}}
        rmsds = {}

        for inactive_id in inactive_ids:
            inactive_coords = inactive_data[inactive_id]

            displacements, rmsd, n_atoms = superimpose_and_calculate(
                active_coords, inactive_coords
            )

            if displacements is None:
                print(f"    vs {inactive_id}: SKIPPED (too few common residues)")
                continue

            per_inactive[inactive_id] = displacements
            rmsds[inactive_id] = rmsd
            print(f"    vs {inactive_id}: RMSD={rmsd:.2f}A, {n_atoms} common residues")

        if not per_inactive:
            print(f"    WARNING: No successful comparisons for {active_id}")
            continue

        # Collect all human positions across all inactive comparisons
        all_positions = set()
        for disps in per_inactive.values():
            all_positions.update(disps.keys())
        all_positions = sorted(all_positions)

        # Build results rows for this active structure
        rows = []
        for pos in all_positions:
            # Get residue name from active structure
            resname = active_coords[pos][1] if pos in active_coords else "UNK"
            grand_residue_names[pos] = resname

            row = {
                'residue_number_human': pos,
                'residue_name': resname,
            }

            values_for_mean = []
            for inactive_id in inactive_ids:
                if inactive_id in per_inactive and pos in per_inactive[inactive_id]:
                    val = per_inactive[inactive_id][pos]
                    row[f'dist_vs_{inactive_id}'] = round(val, 3)
                    values_for_mean.append(val)
                else:
                    row[f'dist_vs_{inactive_id}'] = None

            if values_for_mean:
                mean_val = np.mean(values_for_mean)
                row['mean_displacement'] = round(mean_val, 3)
                grand_displacements[pos].append(mean_val)
            else:
                row['mean_displacement'] = None

            rows.append(row)

        all_active_results[active_id] = rows

        # Write individual CSV
        outfile = os.path.join(output_dir, f"{active_id}_ca_displacements.csv")
        fieldnames = ['residue_number_human', 'residue_name']
        fieldnames += [f'dist_vs_{iid}' for iid in inactive_ids]
        fieldnames += ['mean_displacement']

        with open(outfile, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        mean_rmsd = np.mean(list(rmsds.values()))
        print(f"    Mean RMSD across inactive structures: {mean_rmsd:.2f}A")
        print(f"    Saved: {os.path.basename(outfile)}")

    # ── Step 3: Grand average across all active structures ───────────────────
    print("\n--- Calculating grand average ---")

    grand_positions = sorted(grand_displacements.keys())
    grand_rows = []

    for pos in grand_positions:
        values = grand_displacements[pos]
        resname = grand_residue_names.get(pos, "UNK")
        grand_rows.append({
            'residue_number_human': pos,
            'residue_name': resname,
            'n_active_structures': len(values),
            'mean_ca_displacement': round(np.mean(values), 3),
            'std_ca_displacement': round(np.std(values), 3),
            'min_ca_displacement': round(np.min(values), 3),
            'max_ca_displacement': round(np.max(values), 3),
        })

    grand_file = os.path.join(output_dir, "grand_average_ca_displacements.csv")
    with open(grand_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'residue_number_human', 'residue_name', 'n_active_structures',
            'mean_ca_displacement', 'std_ca_displacement',
            'min_ca_displacement', 'max_ca_displacement'
        ])
        writer.writeheader()
        writer.writerows(grand_rows)

    print(f"  Grand average saved: {grand_file}")
    print(f"  Positions covered: {len(grand_positions)}")

    # ── Step 4: Summary report ───────────────────────────────────────────────
    summary_file = os.path.join(output_dir, "summary_report.txt")
    with open(summary_file, 'w') as f:
        f.write("MOR CA Distance Analysis: Active vs Inactive Structures\n")
        f.write("=" * 80 + "\n\n")

        f.write("Method:\n")
        f.write("  For each active structure, superimpose onto each inactive structure\n")
        f.write("  using SVD-based CA alignment on common residue positions (human numbering).\n")
        f.write("  Then calculate per-residue CA displacement after superposition.\n")
        f.write("  Mouse structures use +2 offset for human numbering alignment.\n\n")

        f.write(f"Inactive structures ({len(inactive_data)}):\n")
        for pid in inactive_ids:
            sp = INACTIVE_STRUCTURES[pid]["species"]
            f.write(f"  {pid} ({sp})\n")

        f.write(f"\nActive structures ({len(active_data)}):\n")
        for pid in sorted(active_data.keys(), key=str.lower):
            sp = ACTIVE_STRUCTURES[pid]["species"]
            f.write(f"  {pid} ({sp})\n")

        f.write(f"\n{'='*80}\n")
        f.write("Top 30 residues with largest mean CA displacement (grand average):\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Position':<10} {'Residue':<8} {'Mean(A)':<10} {'Std(A)':<10} {'Min(A)':<10} {'Max(A)':<10} {'N':<5}\n")
        f.write("-" * 80 + "\n")

        top_movers = sorted(grand_rows, key=lambda x: x['mean_ca_displacement'], reverse=True)
        for row in top_movers[:30]:
            f.write(f"{row['residue_number_human']:<10} "
                    f"{row['residue_name']:<8} "
                    f"{row['mean_ca_displacement']:<10.3f} "
                    f"{row['std_ca_displacement']:<10.3f} "
                    f"{row['min_ca_displacement']:<10.3f} "
                    f"{row['max_ca_displacement']:<10.3f} "
                    f"{row['n_active_structures']:<5}\n")

        f.write(f"\n{'='*80}\n")
        f.write("Per-active-structure mean displacement:\n")
        f.write("-" * 80 + "\n")

        for active_id in sorted(all_active_results.keys(), key=str.lower):
            rows = all_active_results[active_id]
            mean_vals = [r['mean_displacement'] for r in rows if r['mean_displacement'] is not None]
            if mean_vals:
                overall_mean = np.mean(mean_vals)
                overall_max = np.max(mean_vals)
                f.write(f"  {active_id}: mean={overall_mean:.2f}A, max={overall_max:.2f}A, "
                        f"{len(mean_vals)} residues\n")

    print(f"  Summary saved: {summary_file}")

    print(f"\n{'='*80}")
    print("DONE!")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()

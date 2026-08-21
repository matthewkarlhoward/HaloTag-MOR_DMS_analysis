#!/usr/bin/env python3
"""
Analyze MOR structures: G-alpha protein distances + shell annotation.

Calculates minimum all-atom and sidechain distances from every receptor residue
to the G-alpha protein, then classifies residues into:
  First shell:  sidechain atoms ≤ 5.0 Å from any G-alpha atom
  Second shell: sidechain atoms ≤ 4.5 Å from any first-shell residue's sidechain
                (not already first shell)

Uses manually verified receptor chains (single chain per structure) and
verified G-alpha chains. Excludes 8qot and structures without G-alpha.
"""

import os
import numpy as np
import csv
from pathlib import Path

from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

REPO_ROOT = Path(__file__).resolve().parents[3]
CIF_DIR = REPO_ROOT / "structures" / "raw" / "experimental"


# ── Verified structure definitions ───────────────────────────────────────────
# receptor_chain: single verified chain closest to ligand
# galpha_chain: verified G-alpha subunit chain
# species: mouse (+2 offset) or human (no offset)
#
# Excluded: 8qot (apo), 4dkl, 5c1m, 7ul4, 8e0g, 9bjk, 9MQI, 9PXU, 9WSV, 9WSX
#   (no G-alpha protein in complex)
STRUCTURE_DEFS = {
    "6dde": {"receptor_chain": "R", "galpha_chain": "A", "species": "mouse"},
    "6ddf": {"receptor_chain": "R", "galpha_chain": "A", "species": "mouse"},
    "7sbf": {"receptor_chain": "R", "galpha_chain": "A", "species": "mouse"},
    "7scg": {"receptor_chain": "D", "galpha_chain": "A", "species": "mouse"},
    "7t2g": {"receptor_chain": "R", "galpha_chain": "A", "species": "mouse"},
    "7t2h": {"receptor_chain": "D", "galpha_chain": "A", "species": "mouse"},
    "7u2k": {"receptor_chain": "D", "galpha_chain": "A", "species": "mouse"},
    "7u2l": {"receptor_chain": "D", "galpha_chain": "A", "species": "mouse"},
    "8ef5": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "8ef6": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "8efb": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "8efl": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "8efo": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "8efq": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "8f7q": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "8f7r": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "8k9k": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "8k9l": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "8Y72": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "8Y73": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "9bqj": {"receptor_chain": "D", "galpha_chain": "A", "species": "mouse"},
    "9PXV": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "9PXW": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "9PXX": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "9PXY": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "9PY2": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "9PY3": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "9PY4": {"receptor_chain": "R", "galpha_chain": "A", "species": "human"},
    "9WST": {"receptor_chain": "R", "galpha_chain": "A", "species": "mouse"},
    "9WSW": {"receptor_chain": "R", "galpha_chain": "A", "species": "mouse"},
}

# Shell definition thresholds (same as ligand distance analysis)
FIRST_SHELL_CUTOFF = 5.0    # sidechain ≤5Å from G-alpha
SECOND_SHELL_CUTOFF = 4.5   # sidechain ≤4.5Å from first-shell sidechain

BACKBONE_ATOMS = {'N', 'CA', 'C', 'O'}


# ── Distance calculation functions ───────────────────────────────────────────

def get_sidechain_atoms(residue):
    """Get sidechain atoms for a residue. For GLY, returns all atoms."""
    if residue.get_resname() == 'GLY':
        return list(residue.get_atoms())
    return [a for a in residue.get_atoms() if a.get_name() not in BACKBONE_ATOMS]


def calculate_min_distance(residue, target_atoms):
    """Minimum distance between any atom in residue and any target atom."""
    res_atoms = list(residue.get_atoms())
    if not res_atoms or not target_atoms:
        return None
    min_dist = float('inf')
    for ra in res_atoms:
        for ta in target_atoms:
            d = ra - ta
            if d < min_dist:
                min_dist = d
    return min_dist


def calculate_min_sidechain_distance(residue, target_atoms):
    """Minimum distance between sidechain atoms of residue and target atoms."""
    sc_atoms = get_sidechain_atoms(residue)
    if not sc_atoms or not target_atoms:
        return None
    min_dist = float('inf')
    for sa in sc_atoms:
        for ta in target_atoms:
            d = sa - ta
            if d < min_dist:
                min_dist = d
    return min_dist


def calculate_min_sidechain_sidechain_distance(res_a, res_b):
    """Minimum sidechain-sidechain distance between two residues."""
    atoms_a = get_sidechain_atoms(res_a)
    atoms_b = get_sidechain_atoms(res_b)
    if not atoms_a or not atoms_b:
        return float('inf')
    min_dist = float('inf')
    for aa in atoms_a:
        for ab in atoms_b:
            d = aa - ab
            if d < min_dist:
                min_dist = d
    return min_dist


# ── Shell classification ─────────────────────────────────────────────────────

def classify_shells(receptor_residues, galpha_atoms):
    """
    Classify receptor residues into G-protein binding shells.

    First shell:  sidechain atoms ≤ 5.0 Å from any G-alpha atom
    Second shell: NOT first shell, but sidechain atoms ≤ 4.5 Å from any
                  first-shell residue's sidechain atoms

    Returns dict: {(chain_id, resnum): shell_label}
    """
    first_shell_residues = []
    residue_shell = {}

    for res in receptor_residues:
        sc_dist = calculate_min_sidechain_distance(res, galpha_atoms)
        rid = (res.get_parent().id, res.id[1])
        if sc_dist is not None and sc_dist <= FIRST_SHELL_CUTOFF:
            residue_shell[rid] = 'first_shell'
            first_shell_residues.append(res)
        else:
            residue_shell[rid] = 'none'

    print(f"  First shell ({FIRST_SHELL_CUTOFF}A sidechain->G-alpha): {len(first_shell_residues)} residues")

    second_shell_count = 0
    for res in receptor_residues:
        rid = (res.get_parent().id, res.id[1])
        if residue_shell[rid] == 'first_shell':
            continue
        for fs_res in first_shell_residues:
            sc_sc_dist = calculate_min_sidechain_sidechain_distance(res, fs_res)
            if sc_sc_dist <= SECOND_SHELL_CUTOFF:
                residue_shell[rid] = 'second_shell'
                second_shell_count += 1
                break

    print(f"  Second shell ({SECOND_SHELL_CUTOFF}A sidechain->first-shell sidechain): {second_shell_count} residues")

    return residue_shell


# ── Main analysis ────────────────────────────────────────────────────────────

def analyze_structure(cif_file, pdb_id, sdef, output_dir):
    """Analyze one structure for receptor-to-G-alpha distances with shell classification."""
    parser = MMCIFParser(QUIET=True)

    print(f"\n{'='*70}")
    print(f"Analyzing {pdb_id}")
    print(f"{'='*70}")

    try:
        structure = parser.get_structure(pdb_id, cif_file)
    except Exception as e:
        print(f"  ERROR parsing {pdb_id}: {e}")
        return None

    model = structure[0]
    receptor_chain_id = sdef["receptor_chain"]
    galpha_chain_id = sdef["galpha_chain"]
    species = sdef["species"]
    human_offset = 2 if species == "mouse" else 0

    # Validate chains exist
    try:
        receptor_chain = model[receptor_chain_id]
    except KeyError:
        print(f"  WARNING: Receptor chain {receptor_chain_id} not found")
        return None
    try:
        galpha_chain = model[galpha_chain_id]
    except KeyError:
        print(f"  WARNING: G-alpha chain {galpha_chain_id} not found")
        return None

    # Collect all G-alpha atoms (all atoms from all standard AA residues)
    galpha_residues = [r for r in galpha_chain if is_aa(r, standard=True)]
    galpha_atoms = []
    for res in galpha_residues:
        galpha_atoms.extend(list(res.get_atoms()))

    # Collect receptor residues
    receptor_residues = [r for r in receptor_chain if is_aa(r, standard=True)]

    print(f"  Receptor chain {receptor_chain_id}: {len(receptor_residues)} residues")
    print(f"  G-alpha chain {galpha_chain_id}: {len(galpha_residues)} residues, {len(galpha_atoms)} atoms")
    print(f"  Species: {species.upper()} (human offset: +{human_offset})")

    # Classify shells
    shell_map = classify_shells(receptor_residues, galpha_atoms)

    # Calculate distances and build results
    results = []
    for residue in receptor_residues:
        res_num = residue.id[1]
        res_name = residue.get_resname()
        rid = (receptor_chain_id, res_num)

        min_dist = calculate_min_distance(residue, galpha_atoms)
        min_sc_dist = calculate_min_sidechain_distance(residue, galpha_atoms)

        if min_dist is not None:
            results.append({
                'pdb_id': pdb_id,
                'species': species,
                'receptor_chain': receptor_chain_id,
                'residue_number': res_num,
                'residue_number_human': res_num + human_offset,
                'residue_name': res_name,
                'galpha_chain': galpha_chain_id,
                'distance_angstrom': round(min_dist, 3),
                'distance_sidechain_angstrom': round(min_sc_dist, 3) if min_sc_dist else None,
                'shell': shell_map.get(rid, 'none'),
            })

    # Save individual CSV
    if results:
        outfile = os.path.join(output_dir, f"{pdb_id}_gprotein_distances.csv")
        fieldnames = [
            'pdb_id', 'species', 'receptor_chain',
            'residue_number', 'residue_number_human', 'residue_name',
            'galpha_chain', 'distance_angstrom', 'distance_sidechain_angstrom',
            'shell'
        ]
        with open(outfile, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)

        n_first = sum(1 for r in results if r['shell'] == 'first_shell')
        n_second = sum(1 for r in results if r['shell'] == 'second_shell')
        print(f"  Wrote {len(results)} measurements -> {Path(outfile).name}")
        print(f"  Shell summary: {n_first} first, {n_second} second, {len(results)-n_first-n_second} none")
    else:
        print(f"  WARNING: No distance measurements for {pdb_id}")

    return results


def generate_summary(all_results, output_dir):
    """Write combined CSV and text summary."""
    fieldnames = [
        'pdb_id', 'species', 'receptor_chain',
        'residue_number', 'residue_number_human', 'residue_name',
        'galpha_chain', 'distance_angstrom', 'distance_sidechain_angstrom',
        'shell'
    ]

    # Combined CSV
    combined_file = os.path.join(output_dir, "all_structures_combined.csv")
    with open(combined_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for results in all_results:
            writer.writerows(results)
    print(f"\nCombined CSV: {combined_file}")

    # Summary report
    summary_file = os.path.join(output_dir, "summary_report.txt")
    with open(summary_file, 'w') as f:
        f.write("MOR G-alpha Protein Distance Analysis — Experimental Structures\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Total structures analyzed: {len(all_results)}\n")
        total = sum(len(r) for r in all_results)
        f.write(f"Total distance measurements: {total:,}\n\n")

        f.write("Shell definitions:\n")
        f.write(f"  First shell:  sidechain atoms <= {FIRST_SHELL_CUTOFF} A from G-alpha\n")
        f.write(f"  Second shell: sidechain atoms <= {SECOND_SHELL_CUTOFF} A from first-shell sidechain\n")
        f.write(f"                (not already first shell)\n\n")

        f.write("Excluded structures:\n")
        f.write("  8qot — apo (no orthosteric ligand)\n")
        f.write("  4dkl, 5c1m, 7ul4, 8e0g, 9bjk, 9MQI, 9PXU, 9WSV, 9WSX — no G-alpha protein\n\n")

        # Species counts
        mouse_count = sum(1 for r in all_results if r and r[0]['species'] == 'mouse')
        human_count = sum(1 for r in all_results if r and r[0]['species'] == 'human')
        f.write(f"Mouse structures: {mouse_count}\n")
        f.write(f"Human structures: {human_count}\n\n")
        f.write("-" * 80 + "\n")

        # Per-structure details
        for results in all_results:
            pdb_id = results[0]['pdb_id']
            species = results[0]['species']
            galpha = results[0]['galpha_chain']

            distances = [r['distance_angstrom'] for r in results]
            first_shell = sorted(
                [r for r in results if r['shell'] == 'first_shell'],
                key=lambda x: x['distance_sidechain_angstrom'] or 999
            )
            second_shell = sorted(
                [r for r in results if r['shell'] == 'second_shell'],
                key=lambda x: x['distance_sidechain_angstrom'] or 999
            )

            f.write(f"\n{pdb_id}:\n")
            f.write(f"  Species: {species.upper()}\n")
            f.write(f"  G-alpha chain: {galpha}\n")
            f.write(f"  Distance range: {min(distances):.2f} -- {max(distances):.2f} A\n")
            f.write(f"  Mean distance: {np.mean(distances):.2f} A\n")
            f.write(f"  First shell:  {len(first_shell)} residues\n")
            f.write(f"  Second shell: {len(second_shell)} residues\n")

            if first_shell:
                res_str = ", ".join(
                    f"{r['residue_name']}{r['residue_number_human']}({r['distance_sidechain_angstrom']:.1f}A)"
                    for r in first_shell
                )
                f.write(f"  First shell residues:  {res_str}\n")

            if second_shell:
                res_str = ", ".join(
                    f"{r['residue_name']}{r['residue_number_human']}"
                    for r in second_shell
                )
                f.write(f"  Second shell residues: {res_str}\n")

    print(f"Summary report: {summary_file}")


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cif_dir = str(CIF_DIR)
    output_dir = str(REPO_ROOT / "structures" / "processed" / "gprotein_distances")
    os.makedirs(output_dir, exist_ok=True)

    print(f"CIF directory:       {cif_dir}")
    print(f"Output directory:    {output_dir}")
    print(f"Structures to analyze: {len(STRUCTURE_DEFS)}")
    print(f"First shell cutoff:  {FIRST_SHELL_CUTOFF} A (sidechain -> G-alpha)")
    print(f"Second shell cutoff: {SECOND_SHELL_CUTOFF} A (sidechain -> first-shell sidechain)")

    all_results = []
    failed = []

    for pdb_id, sdef in sorted(STRUCTURE_DEFS.items(), key=lambda x: x[0].lower()):
        cif_path = os.path.join(cif_dir, f"{pdb_id}.cif")
        if not os.path.exists(cif_path):
            print(f"\n  WARNING: {cif_path} not found, skipping {pdb_id}")
            failed.append(pdb_id)
            continue

        results = analyze_structure(cif_path, pdb_id, sdef, output_dir)
        if results:
            all_results.append(results)
        else:
            failed.append(pdb_id)

    if all_results:
        generate_summary(all_results, output_dir)

    print(f"\n{'='*70}")
    print(f"DONE: {len(all_results)} structures analyzed, {len(failed)} failed")
    if failed:
        print(f"Failed: {', '.join(failed)}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()

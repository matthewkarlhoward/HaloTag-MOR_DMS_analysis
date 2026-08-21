#!/usr/bin/env python3
"""
Analyze ChAI-predicted MOR structures: ligand distances + shell annotation.

Same methodology as experimental structure analysis:
  First shell:  sidechain atoms ≤ 5.0 Å from ligand
  Second shell: sidechain-sidechain ≤ 4.5 Å from first-shell residue (not already first shell)

ChAI predictions use residue numbering 1-400 for the receptor (chain A).
User requests: ignore receptor residues 1-60 (disordered N-terminus).
Human MOR numbering is used directly (no offset needed - these are predictions).
"""

import os
import numpy as np
import csv
from pathlib import Path

from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

REPO_ROOT = Path(__file__).resolve().parents[3]
CHAI_DIR = REPO_ROOT / "structures" / "raw" / "chai_predicted"

# ── Structure definitions ────────────────────────────────────────────────────
STRUCTURE_DEFS = {
    "MOR-Naltrexone_rank1": {
        "drug": "Naltrexone",
        "cif_path": "MOR-Naltrexone/pred.rank_1.cif",
        "ligand_code": "LIG2",
        "receptor_chain": "A",
        "min_residue": 61,  # skip 1-60
    },
    "MOR-Nalbuphine_rank0": {
        "drug": "Nalbuphine",
        "cif_path": "MOR-GNAI-Nalbuphine/pred.rank_0_nalbuphine.cif",
        "ligand_code": "LIG3",
        "receptor_chain": "A",
        "min_residue": 61,
    },
    "MOR-Methadone_rank0": {
        "drug": "Methadone",
        "cif_path": "MOR-GNAI-Methadone/pred.rank_0_methadone.cif",
        "ligand_code": "LIG3",
        "receptor_chain": "A",
        "min_residue": 61,
    },
    "MOR-Butorphanol_rank0": {
        "drug": "Butorphanol",
        "cif_path": "MOR-GNAI-Butorphanol/pred.rank_0_butorphanol.cif",
        "ligand_code": "LIG3",
        "receptor_chain": "A",
        "min_residue": 61,
    },
}

# Shell thresholds (same as experimental)
FIRST_SHELL_CUTOFF = 5.0
SECOND_SHELL_CUTOFF = 4.5

BACKBONE_ATOMS = {'N', 'CA', 'C', 'O'}


def get_sidechain_atoms(residue):
    if residue.get_resname() == 'GLY':
        return list(residue.get_atoms())
    return [a for a in residue.get_atoms() if a.get_name() not in BACKBONE_ATOMS]


def calculate_min_distance(residue, ligand_atoms):
    res_atoms = list(residue.get_atoms())
    if not res_atoms or not ligand_atoms:
        return None
    min_dist = float('inf')
    for ra in res_atoms:
        for la in ligand_atoms:
            d = ra - la
            if d < min_dist:
                min_dist = d
    return min_dist


def calculate_min_sidechain_distance(residue, ligand_atoms):
    res_atoms = get_sidechain_atoms(residue)
    if not res_atoms or not ligand_atoms:
        return None
    min_dist = float('inf')
    for ra in res_atoms:
        for la in ligand_atoms:
            d = ra - la
            if d < min_dist:
                min_dist = d
    return min_dist


def calculate_min_sidechain_sidechain_distance(res_a, res_b):
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


def classify_shells(receptor_residues, ligand_atoms):
    first_shell_residues = []
    residue_shell = {}

    for res in receptor_residues:
        sc_dist = calculate_min_sidechain_distance(res, ligand_atoms)
        rid = (res.get_parent().id, res.id[1])
        if sc_dist is not None and sc_dist <= FIRST_SHELL_CUTOFF:
            residue_shell[rid] = 'first_shell'
            first_shell_residues.append(res)
        else:
            residue_shell[rid] = 'none'

    print(f"  First shell ({FIRST_SHELL_CUTOFF}A): {len(first_shell_residues)} residues")

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

    print(f"  Second shell ({SECOND_SHELL_CUTOFF}A): {second_shell_count} residues")
    return residue_shell


def identify_gprotein_chains(structure, exclude_chains):
    model = structure[0]
    gprotein_chains = []
    exclude = set(exclude_chains) if exclude_chains else set()
    for chain in model:
        cid = chain.get_id()
        if cid in exclude:
            continue
        aa_count = sum(1 for r in chain if is_aa(r, standard=True))
        if aa_count > 50:
            gprotein_chains.append(cid)
    return gprotein_chains


def analyze_structure(struct_id, sdef, base_dir, output_dir):
    parser = MMCIFParser(QUIET=True)
    cif_file = os.path.join(base_dir, sdef['cif_path'])

    print(f"\n{'='*70}")
    print(f"Analyzing {struct_id}  --  {sdef['drug']}")
    print(f"{'='*70}")

    try:
        structure = parser.get_structure(struct_id, cif_file)
    except Exception as e:
        print(f"  ERROR: {e}")
        return None

    model = structure[0]
    receptor_chain_id = sdef['receptor_chain']
    ligand_code = sdef['ligand_code']
    min_res = sdef.get('min_residue', 1)

    # Collect ligand atoms
    ligand_atoms = []
    for chain in model:
        for residue in chain:
            if residue.get_resname() == ligand_code:
                ligand_atoms.extend(list(residue.get_atoms()))

    if not ligand_atoms:
        print(f"  WARNING: No atoms found for ligand {ligand_code}")
        return None
    print(f"  Found {len(ligand_atoms)} atoms for ligand {ligand_code}")

    # G-protein chains
    gprotein_chains = identify_gprotein_chains(structure, [receptor_chain_id])
    print(f"  Receptor chain: {receptor_chain_id}, G-protein chains: {gprotein_chains}")
    print(f"  Skipping receptor residues < {min_res}")

    # Collect receptor residues (skip 1-60)
    receptor_residues = [
        res for res in model[receptor_chain_id]
        if is_aa(res, standard=True) and res.id[1] >= min_res
    ]
    print(f"  Receptor residues: {len(receptor_residues)} (from {min_res})")

    # Classify shells
    shell_map = classify_shells(receptor_residues, ligand_atoms)

    # Calculate distances
    results = []
    for residue in receptor_residues:
        res_num = residue.id[1]
        res_name = residue.get_resname()
        rid = (receptor_chain_id, res_num)

        min_dist = calculate_min_distance(residue, ligand_atoms)
        min_sc_dist = calculate_min_sidechain_distance(residue, ligand_atoms)

        if min_dist is not None:
            results.append({
                'pdb_id': struct_id,
                'species': 'human_prediction',
                'drug_name': sdef['drug'],
                'receptor_chain': receptor_chain_id,
                'residue_number': res_num,
                'residue_number_human': res_num,  # ChAI uses human numbering directly
                'residue_name': res_name,
                'ligand_name': ligand_code,
                'distance_angstrom': round(min_dist, 3),
                'distance_sidechain_angstrom': round(min_sc_dist, 3) if min_sc_dist else None,
                'shell': shell_map.get(rid, 'none'),
                'gprotein_chains': ','.join(gprotein_chains),
            })

    # Save CSV
    if results:
        outfile = os.path.join(output_dir, f"{struct_id}_distances.csv")
        fieldnames = [
            'pdb_id', 'species', 'drug_name', 'receptor_chain',
            'residue_number', 'residue_number_human', 'residue_name',
            'ligand_name', 'distance_angstrom', 'distance_sidechain_angstrom',
            'shell', 'gprotein_chains'
        ]
        with open(outfile, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)

        n_first = sum(1 for r in results if r['shell'] == 'first_shell')
        n_second = sum(1 for r in results if r['shell'] == 'second_shell')
        print(f"  Wrote {len(results)} measurements -> {Path(outfile).name}")
        print(f"  Shell summary: {n_first} first, {n_second} second, {len(results)-n_first-n_second} none")

    return results


def generate_summary(all_results, output_dir):
    fieldnames = [
        'pdb_id', 'species', 'drug_name', 'receptor_chain',
        'residue_number', 'residue_number_human', 'residue_name',
        'ligand_name', 'distance_angstrom', 'distance_sidechain_angstrom',
        'shell', 'gprotein_chains'
    ]

    combined_file = os.path.join(output_dir, "chai_predictions_combined.csv")
    with open(combined_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for results in all_results:
            writer.writerows(results)
    print(f"\nCombined CSV: {combined_file}")

    summary_file = os.path.join(output_dir, "summary_report.txt")
    with open(summary_file, 'w') as f:
        f.write("MOR Ligand Distance Analysis — ChAI Predictions\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Total structures analyzed: {len(all_results)}\n")
        total = sum(len(r) for r in all_results)
        f.write(f"Total distance measurements: {total:,}\n\n")
        f.write("Shell definitions:\n")
        f.write(f"  First shell:  sidechain atoms <= {FIRST_SHELL_CUTOFF} A from ligand\n")
        f.write(f"  Second shell: sidechain atoms <= {SECOND_SHELL_CUTOFF} A from first-shell sidechain\n")
        f.write(f"                (not already first shell)\n\n")
        f.write("Note: Receptor residues 1-60 excluded (disordered N-terminus in predictions)\n")
        f.write("Note: ChAI predictions use human MOR numbering directly\n\n")
        f.write("-" * 80 + "\n")

        for results in all_results:
            sid = results[0]['pdb_id']
            drug = results[0]['drug_name']
            gp = results[0]['gprotein_chains']
            distances = [r['distance_angstrom'] for r in results]
            first_shell = sorted(
                [r for r in results if r['shell'] == 'first_shell'],
                key=lambda x: x['distance_sidechain_angstrom'] or 999
            )
            second_shell = sorted(
                [r for r in results if r['shell'] == 'second_shell'],
                key=lambda x: x['distance_sidechain_angstrom'] or 999
            )

            f.write(f"\n{sid}:\n")
            f.write(f"  Drug: {drug}\n")
            f.write(f"  G-protein chains: {gp or 'none'}\n")
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
    base_dir = str(CHAI_DIR)
    output_dir = str(REPO_ROOT / "structures" / "processed" / "chai_distances")
    os.makedirs(output_dir, exist_ok=True)

    print(f"Base directory: {base_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Structures to analyze: {len(STRUCTURE_DEFS)}")

    all_results = []
    failed = []

    for struct_id, sdef in STRUCTURE_DEFS.items():
        results = analyze_structure(struct_id, sdef, base_dir, output_dir)
        if results:
            all_results.append(results)
        else:
            failed.append(struct_id)

    if all_results:
        generate_summary(all_results, output_dir)

    print(f"\n{'='*70}")
    print(f"DONE: {len(all_results)} analyzed, {len(failed)} failed")
    if failed:
        print(f"Failed: {', '.join(failed)}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()

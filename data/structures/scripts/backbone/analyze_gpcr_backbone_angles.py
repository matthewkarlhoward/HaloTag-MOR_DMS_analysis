#!/usr/bin/env python3
"""
Analyze GPCR structures to calculate backbone dihedral angles (phi, psi, omega)
for each receptor residue.

Uses manually verified single receptor chains (from STRUCTURE_DEFS) to avoid
picking up duplicate chains in structures with multiple receptor copies.
"""

import os
import numpy as np
import csv
from pathlib import Path

from Bio.PDB import MMCIFParser, calc_dihedral
from Bio.PDB.Polypeptide import is_aa

REPO_ROOT = Path(__file__).resolve().parents[3]
CIF_DIR = REPO_ROOT / "structures" / "raw" / "experimental"

# ── Verified structure definitions ───────────────────────────────────────────
# receptor_chain: single verified chain closest to the ligand
# species: mouse (+2 offset) or human (no offset)
# Excludes 8qot per user request.
STRUCTURE_DEFS = {
    "4dkl": {"receptor_chain": "A", "species": "mouse"},
    "5c1m": {"receptor_chain": "A", "species": "mouse"},
    "6dde": {"receptor_chain": "R", "species": "mouse"},
    "6ddf": {"receptor_chain": "R", "species": "mouse"},
    "7sbf": {"receptor_chain": "R", "species": "mouse"},
    "7scg": {"receptor_chain": "D", "species": "mouse"},
    "7t2g": {"receptor_chain": "R", "species": "mouse"},
    "7t2h": {"receptor_chain": "D", "species": "mouse"},
    "7u2k": {"receptor_chain": "D", "species": "mouse"},
    "7u2l": {"receptor_chain": "D", "species": "mouse"},
    "7ul4": {"receptor_chain": "A", "species": "mouse"},
    "8e0g": {"receptor_chain": "A", "species": "mouse"},
    "8ef5": {"receptor_chain": "R", "species": "human"},
    "8ef6": {"receptor_chain": "R", "species": "human"},
    "8efb": {"receptor_chain": "R", "species": "human"},
    "8efl": {"receptor_chain": "R", "species": "human"},
    "8efo": {"receptor_chain": "R", "species": "human"},
    "8efq": {"receptor_chain": "R", "species": "human"},
    "8f7q": {"receptor_chain": "R", "species": "human"},
    "8f7r": {"receptor_chain": "R", "species": "human"},
    "8k9k": {"receptor_chain": "R", "species": "human"},
    "8k9l": {"receptor_chain": "R", "species": "human"},
    "8Y72": {"receptor_chain": "R", "species": "human"},
    "8Y73": {"receptor_chain": "R", "species": "human"},
    "9bjk": {"receptor_chain": "R", "species": "mouse"},
    "9bqj": {"receptor_chain": "D", "species": "mouse"},
    "9MQI": {"receptor_chain": "A", "species": "human"},
    "9PXU": {"receptor_chain": "R", "species": "human"},
    "9PXV": {"receptor_chain": "R", "species": "human"},
    "9PXW": {"receptor_chain": "R", "species": "human"},
    "9PXX": {"receptor_chain": "R", "species": "human"},
    "9PXY": {"receptor_chain": "R", "species": "human"},
    "9PY2": {"receptor_chain": "R", "species": "human"},
    "9PY3": {"receptor_chain": "R", "species": "human"},
    "9PY4": {"receptor_chain": "R", "species": "human"},
    "9WST": {"receptor_chain": "R", "species": "mouse"},
    "9WSV": {"receptor_chain": "R", "species": "mouse"},
    "9WSW": {"receptor_chain": "R", "species": "mouse"},
    "9WSX": {"receptor_chain": "R", "species": "mouse"},
}


def calculate_backbone_angles(residue_list, residue_index):
    """
    Calculate phi, psi, and omega angles for a given residue.

    Phi:   C(i-1) - N(i) - CA(i) - C(i)
    Psi:   N(i) - CA(i) - C(i) - N(i+1)
    Omega: CA(i) - C(i) - N(i+1) - CA(i+1)

    Returns (phi, psi, omega) in degrees, or None for missing atoms.
    """
    current = residue_list[residue_index]
    phi = psi = omega = None

    # PHI (requires previous residue)
    if residue_index > 0:
        prev = residue_list[residue_index - 1]
        try:
            phi = np.degrees(calc_dihedral(
                prev['C'].get_vector(),
                current['N'].get_vector(),
                current['CA'].get_vector(),
                current['C'].get_vector()
            ))
        except KeyError:
            pass

    # PSI (requires next residue)
    if residue_index < len(residue_list) - 1:
        next_res = residue_list[residue_index + 1]
        try:
            psi = np.degrees(calc_dihedral(
                current['N'].get_vector(),
                current['CA'].get_vector(),
                current['C'].get_vector(),
                next_res['N'].get_vector()
            ))
        except KeyError:
            pass

    # OMEGA (requires next residue)
    if residue_index < len(residue_list) - 1:
        next_res = residue_list[residue_index + 1]
        try:
            omega = np.degrees(calc_dihedral(
                current['CA'].get_vector(),
                current['C'].get_vector(),
                next_res['N'].get_vector(),
                next_res['CA'].get_vector()
            ))
        except KeyError:
            pass

    return phi, psi, omega


def analyze_structure(cif_file, pdb_id, sdef, output_dir):
    """Analyze a single structure using its verified chain definition."""
    parser = MMCIFParser(QUIET=True)

    print(f"\nAnalyzing {pdb_id}...")

    try:
        structure = parser.get_structure(pdb_id, cif_file)
    except Exception as e:
        print(f"  Error parsing {pdb_id}: {e}")
        return None

    model = structure[0]
    chain_id = sdef["receptor_chain"]
    species = sdef["species"]
    human_offset = 2 if species == "mouse" else 0

    try:
        chain = model[chain_id]
    except KeyError:
        print(f"  WARNING: Chain {chain_id} not found in {pdb_id}")
        return None

    aa_residues = [r for r in chain if is_aa(r, standard=True)]
    print(f"  Chain {chain_id}, species={species.upper()}, {len(aa_residues)} residues, offset=+{human_offset}")

    results = []
    for idx, residue in enumerate(aa_residues):
        res_num = residue.id[1]
        phi, psi, omega = calculate_backbone_angles(aa_residues, idx)

        results.append({
            'pdb_id': pdb_id,
            'species': species,
            'receptor_chain': chain_id,
            'residue_number': res_num,
            'residue_number_human': res_num + human_offset,
            'residue_name': residue.get_resname(),
            'phi': round(phi, 3) if phi is not None else None,
            'psi': round(psi, 3) if psi is not None else None,
            'omega': round(omega, 3) if omega is not None else None,
        })

    if results:
        output_file = os.path.join(output_dir, f"{pdb_id}_backbone_angles.csv")
        fieldnames = ['pdb_id', 'species', 'receptor_chain', 'residue_number',
                      'residue_number_human', 'residue_name', 'phi', 'psi', 'omega']
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        print(f"  Saved {len(results)} angle measurements")

    return results


def generate_summary_report(all_results, output_dir):
    """Generate a summary report of all structures analyzed."""
    if not all_results:
        print("No results to summarize.")
        return

    # Combined CSV
    combined_file = os.path.join(output_dir, "all_structures_combined.csv")
    fieldnames = ['pdb_id', 'species', 'receptor_chain', 'residue_number',
                  'residue_number_human', 'residue_name', 'phi', 'psi', 'omega']
    with open(combined_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for results in all_results:
            writer.writerows(results)
    print(f"\nCombined results saved to {combined_file}")

    # Summary report
    summary_file = os.path.join(output_dir, "summary_report.txt")
    with open(summary_file, 'w') as f:
        f.write("GPCR μOR - Backbone Dihedral Angles Analysis - Summary Report\n")
        f.write("=" * 80 + "\n\n")
        f.write("Uses verified single receptor chains from STRUCTURE_DEFS.\n")
        f.write("8qot excluded from analysis.\n\n")

        f.write(f"Total structures analyzed: {len(all_results)}\n")
        total_measurements = sum(len(r) for r in all_results)
        f.write(f"Total residue angle measurements: {total_measurements}\n\n")

        mouse_count = sum(1 for r in all_results if r and r[0]['species'] == 'mouse')
        human_count = sum(1 for r in all_results if r and r[0]['species'] == 'human')
        f.write(f"Mouse structures: {mouse_count}\n")
        f.write(f"Human structures: {human_count}\n\n")

        f.write("Per-Structure Summary:\n")
        f.write("-" * 80 + "\n")

        for results in all_results:
            if not results:
                continue
            pdb_id = results[0]['pdb_id']
            species = results[0]['species']
            chain = results[0]['receptor_chain']

            phi_count = sum(1 for r in results if r['phi'] is not None)
            psi_count = sum(1 for r in results if r['psi'] is not None)
            omega_count = sum(1 for r in results if r['omega'] is not None)

            phi_values = [r['phi'] for r in results if r['phi'] is not None]
            psi_values = [r['psi'] for r in results if r['psi'] is not None]
            omega_values = [r['omega'] for r in results if r['omega'] is not None]

            f.write(f"\n{pdb_id}:\n")
            f.write(f"  Species: {species.upper()}\n")
            f.write(f"  Receptor chain: {chain}\n")
            f.write(f"  Total residues: {len(results)}\n")
            f.write(f"  Phi angles: {phi_count}\n")
            f.write(f"  Psi angles: {psi_count}\n")
            f.write(f"  Omega angles: {omega_count}\n")
            if phi_values:
                f.write(f"  Phi range: {min(phi_values):.1f}° to {max(phi_values):.1f}°\n")
            if psi_values:
                f.write(f"  Psi range: {min(psi_values):.1f}° to {max(psi_values):.1f}°\n")
            if omega_values:
                f.write(f"  Omega range: {min(omega_values):.1f}° to {max(omega_values):.1f}°\n")

    print(f"Summary report saved to {summary_file}")


def main():
    """Main analysis pipeline."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cif_dir = str(CIF_DIR)
    output_dir = str(REPO_ROOT / "structures" / "processed" / "backbone_angles")
    os.makedirs(output_dir, exist_ok=True)

    print(f"CIF directory: {cif_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Structures to analyze: {len(STRUCTURE_DEFS)}")

    all_results = []
    for pdb_id, sdef in sorted(STRUCTURE_DEFS.items(), key=lambda x: x[0].lower()):
        # Try case-sensitive filename first, then lowercase
        cif_path = os.path.join(cif_dir, f"{pdb_id}.cif")
        if not os.path.exists(cif_path):
            cif_path = os.path.join(cif_dir, f"{pdb_id.lower()}.cif")
        if not os.path.exists(cif_path):
            cif_path = os.path.join(cif_dir, f"{pdb_id.upper()}.cif")
        if not os.path.exists(cif_path):
            print(f"\n  WARNING: CIF file not found for {pdb_id}")
            continue

        results = analyze_structure(cif_path, pdb_id, sdef, output_dir)
        if results:
            all_results.append(results)

    generate_summary_report(all_results, output_dir)

    print("\n" + "=" * 80)
    print(f"Analysis complete! {len(all_results)} structures processed.")
    print(f"Results saved in: {output_dir}")


if __name__ == "__main__":
    main()

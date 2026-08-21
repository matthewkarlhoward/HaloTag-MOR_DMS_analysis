#!/usr/bin/env python3
"""
Analyze ALL experimental MOR structures: ligand distances + shell annotation.

Calculates minimum all-atom and sidechain distances from every receptor residue
to the bound ligand, then classifies residues into:
  First shell:  sidechain atoms ≤ 5.0 Å from ligand
  Second shell: sidechain atoms ≤ 4.5 Å from any first-shell residue's sidechain
                (not already first shell)

Handles both small-molecule and peptide ligands (DAMGO, endomorphins).
Uses manually verified ligand codes, ligand chains, and receptor chains.

Skipped structures:
  8qot — apo (antibody-bound, no orthosteric ligand)
"""

import os
import numpy as np
import csv
from pathlib import Path

from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

REPO_ROOT = Path(__file__).resolve().parents[3]
CIF_DIR = REPO_ROOT / "structures" / "raw" / "experimental"

# ── Structure definitions ────────────────────────────────────────────────────
# Each entry specifies:
#   drug:           common drug name
#   type:           "small_mol" or "peptide"
#   ligand_code:    PDB 3-letter code for small molecules
#   ligand_chains:  list of chain IDs containing the peptide ligand
#   receptor_chain: single chain ID for the receptor (verified closest to ligand)

STRUCTURE_DEFS = {
    # --- beta-FNA (covalent antagonist) ---
    "4dkl": {"drug": "beta-FNA", "type": "small_mol", "ligand_code": "BF0", "receptor_chain": "A"},

    # --- BU72 (agonist) ---
    "5c1m": {"drug": "BU72", "type": "small_mol", "ligand_code": "VF1", "receptor_chain": "A"},

    # --- DAMGO (peptide agonist: Tyr-D-Ala-Gly-MePhe-Gly-ol) ---
    "6dde": {"drug": "DAMGO", "type": "peptide", "ligand_chains": ["D"], "receptor_chain": "R"},
    "6ddf": {"drug": "DAMGO", "type": "peptide", "ligand_chains": ["D"], "receptor_chain": "R"},
    "8efq": {"drug": "DAMGO", "type": "peptide", "ligand_chains": ["P"], "receptor_chain": "R"},
    "8k9k": {"drug": "DAMGO", "type": "peptide", "ligand_chains": ["S"], "receptor_chain": "R"},
    "8k9l": {"drug": "DAMGO", "type": "peptide", "ligand_chains": ["S"], "receptor_chain": "R"},
    "8Y72": {"drug": "DAMGO", "type": "peptide", "ligand_chains": ["D"], "receptor_chain": "R"},
    "8Y73": {"drug": "DAMGO", "type": "peptide", "ligand_chains": ["D"], "receptor_chain": "R"},
    "9WST": {"drug": "DAMGO", "type": "peptide", "ligand_chains": ["P"], "receptor_chain": "R"},
    "9WSV": {"drug": "DAMGO", "type": "peptide", "ligand_chains": ["P"], "receptor_chain": "R"},

    # --- Endomorphin-1/2 and analogs (peptide ligands, all standard AA + NH2 cap) ---
    "8f7q": {"drug": "Endomorphin-1_analog", "type": "peptide", "ligand_chains": ["P"], "receptor_chain": "R"},
    "8f7r": {"drug": "Endomorphin-2", "type": "peptide", "ligand_chains": ["P"], "receptor_chain": "R"},
    "9WSW": {"drug": "Endomorphin-2", "type": "peptide", "ligand_chains": ["P"], "receptor_chain": "R"},
    "9WSX": {"drug": "Endomorphin-2", "type": "peptide", "ligand_chains": ["P"], "receptor_chain": "R"},

    # --- Fentanyl ---
    "8ef5": {"drug": "Fentanyl", "type": "small_mol", "ligand_code": "7V7", "receptor_chain": "R"},

    # --- Morphine ---
    "8ef6": {"drug": "Morphine", "type": "small_mol", "ligand_code": "MOI", "receptor_chain": "R"},

    # --- Oliceridine (TRV130) ---
    "8efb": {"drug": "Oliceridine", "type": "small_mol", "ligand_code": "WH2", "receptor_chain": "R"},

    # --- SR17018 ---
    "8efl": {"drug": "SR17018", "type": "small_mol", "ligand_code": "WH9", "receptor_chain": "R"},

    # --- PZM21 ---
    "7sbf": {"drug": "PZM21", "type": "small_mol", "ligand_code": "8QY", "receptor_chain": "R"},
    "8efo": {"drug": "PZM21", "type": "small_mol", "ligand_code": "8QY", "receptor_chain": "R"},

    # --- Mitragynine Pseudoindoxyl ---
    "7t2g": {"drug": "Mitragynine_Pseudoindoxyl", "type": "small_mol", "ligand_code": "EIG", "receptor_chain": "R"},

    # --- Lofentanil ---
    "7t2h": {"drug": "Lofentanil", "type": "small_mol", "ligand_code": "EID", "receptor_chain": "D"},

    # --- C6-guano ---
    "7u2k": {"drug": "C6-guano", "type": "small_mol", "ligand_code": "KZR", "receptor_chain": "D"},

    # --- C5-guano ---
    "7u2l": {"drug": "C5-guano", "type": "small_mol", "ligand_code": "L0X", "receptor_chain": "D"},

    # --- FH210 ---
    "7scg": {"drug": "FH210", "type": "small_mol", "ligand_code": "8RI", "receptor_chain": "D"},

    # --- Naltrexone-derived ---
    "7ul4": {"drug": "Naltrexone-derived", "type": "small_mol", "ligand_code": "NG0", "receptor_chain": "A"},

    # --- BU72 analog ---
    "8e0g": {"drug": "BU72_analog", "type": "small_mol", "ligand_code": "A1A7R", "receptor_chain": "A"},

    # --- Naloxone ---
    "9PXU": {"drug": "Naloxone", "type": "small_mol", "ligand_code": "A1APV", "receptor_chain": "R"},
    "9PXV": {"drug": "Naloxone", "type": "small_mol", "ligand_code": "A1APV", "receptor_chain": "R"},
    "9PXW": {"drug": "Naloxone", "type": "small_mol", "ligand_code": "A1APV", "receptor_chain": "R"},
    "9PXX": {"drug": "Naloxone", "type": "small_mol", "ligand_code": "A1APV", "receptor_chain": "R"},
    "9PXY": {"drug": "Naloxone", "type": "small_mol", "ligand_code": "A1APV", "receptor_chain": "R"},

    # --- 9bjk: Naloxone + bitopic compound A1APU (use A1APU as primary) ---
    "9bjk": {"drug": "A1APU_bitopic", "type": "small_mol", "ligand_code": "A1APU", "receptor_chain": "R"},

    # --- RO76 ---
    "9bqj": {"drug": "RO76", "type": "small_mol", "ligand_code": "A1AQ2", "receptor_chain": "D"},

    # --- 9MQI: unknown ligand A1BNM ---
    "9MQI": {"drug": "A1BNM", "type": "small_mol", "ligand_code": "A1BNM", "receptor_chain": "A"},

    # --- 9PY2/3/4: compound A1CMV ---
    "9PY2": {"drug": "A1CMV", "type": "small_mol", "ligand_code": "A1CMV", "receptor_chain": "R"},
    "9PY3": {"drug": "A1CMV", "type": "small_mol", "ligand_code": "A1CMV", "receptor_chain": "R"},
    "9PY4": {"drug": "A1CMV", "type": "small_mol", "ligand_code": "A1CMV", "receptor_chain": "R"},

    # SKIPPED: 8qot (apo, no orthosteric ligand)
}

# Shell definition thresholds
FIRST_SHELL_CUTOFF = 5.0    # sidechain ≤5Å from ligand
SECOND_SHELL_CUTOFF = 4.5   # sidechain ≤4.5Å from first-shell sidechain


# ── Distance calculation functions ───────────────────────────────────────────

BACKBONE_ATOMS = {'N', 'CA', 'C', 'O'}


def get_sidechain_atoms(residue):
    """Get sidechain atoms for a residue. For GLY, returns all atoms."""
    if residue.get_resname() == 'GLY':
        return list(residue.get_atoms())
    return [a for a in residue.get_atoms() if a.get_name() not in BACKBONE_ATOMS]


def calculate_min_distance(residue, ligand_atoms):
    """Minimum distance between any atom in residue and any ligand atom."""
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
    """Minimum distance between sidechain atoms and ligand atoms."""
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


# ── Species detection ────────────────────────────────────────────────────────

def determine_species(structure, receptor_chain_id):
    """Mouse MOR has ASP at position 147 (human = 149). Returns species and offset."""
    model = structure[0]
    has_asp147 = False
    has_asp149 = False
    for res in model[receptor_chain_id]:
        if res.get_resname() == 'ASP':
            rnum = res.id[1]
            if rnum == 147:
                has_asp147 = True
            elif rnum == 149:
                has_asp149 = True
    if has_asp147 and not has_asp149:
        return 'mouse', 2
    elif has_asp149:
        return 'human', 0
    return 'unknown', 0


# ── Chain identification ─────────────────────────────────────────────────────

def identify_gprotein_chains(structure, exclude_chains):
    """Identify G-protein chains (excluding receptor and ligand chains)."""
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


# ── Ligand atom collection ───────────────────────────────────────────────────

def collect_ligand_atoms_small_mol(structure, ligand_code):
    """Collect all atoms from HETATM residues matching ligand_code."""
    model = structure[0]
    atoms = []
    for chain in model:
        for residue in chain:
            if residue.get_resname() == ligand_code:
                atoms.extend(list(residue.get_atoms()))
    if not atoms:
        print(f"  WARNING: No atoms found for ligand code {ligand_code}")
    else:
        print(f"  Found {len(atoms)} atoms for ligand {ligand_code}")
    return atoms


def collect_ligand_atoms_peptide(structure, ligand_chain_ids):
    """Collect ALL atoms from specified chains (for peptide ligands)."""
    model = structure[0]
    atoms = []
    for cid in ligand_chain_ids:
        try:
            chain = model[cid]
            for residue in chain:
                atoms.extend(list(residue.get_atoms()))
        except KeyError:
            print(f"  WARNING: Chain {cid} not found")
    if not atoms:
        print(f"  WARNING: No atoms found on chains {ligand_chain_ids}")
    else:
        print(f"  Found {len(atoms)} atoms on peptide chains {ligand_chain_ids}")
    return atoms


# ── Shell classification ─────────────────────────────────────────────────────

def classify_shells(receptor_residues, ligand_atoms):
    """
    Classify receptor residues into binding shells.

    First shell:  sidechain atoms ≤ 5.0 Å from any ligand atom
    Second shell: NOT first shell, but sidechain atoms ≤ 4.5 Å from any
                  first-shell residue's sidechain atoms

    Returns dict: {(chain_id, resnum): shell_label}
    """
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

    print(f"  First shell ({FIRST_SHELL_CUTOFF}A sidechain->ligand): {len(first_shell_residues)} residues")

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

def analyze_structure(cif_file, sdef, output_dir):
    """Analyze one structure with specified ligand definition."""
    parser = MMCIFParser(QUIET=True)
    pdb_id = Path(cif_file).stem

    print(f"\n{'='*70}")
    print(f"Analyzing {pdb_id}  --  {sdef['drug']}")
    print(f"{'='*70}")

    try:
        structure = parser.get_structure(pdb_id, cif_file)
    except Exception as e:
        print(f"  ERROR parsing {pdb_id}: {e}")
        return None

    # Collect ligand atoms
    if sdef['type'] == 'peptide':
        ligand_chain_ids = sdef['ligand_chains']
        ligand_atoms = collect_ligand_atoms_peptide(structure, ligand_chain_ids)
        ligand_label = sdef['drug']
    else:
        ligand_chain_ids = None
        ligand_atoms = collect_ligand_atoms_small_mol(structure, sdef['ligand_code'])
        ligand_label = sdef['ligand_code']

    if not ligand_atoms:
        return None

    # Use specified receptor chain
    receptor_chain_id = sdef['receptor_chain']

    # Identify G-protein chains for reporting
    exclude = [receptor_chain_id] + (sdef.get('ligand_chains') or [])
    gprotein_chains = identify_gprotein_chains(structure, exclude)
    print(f"  Receptor chain: {receptor_chain_id}")
    print(f"  G-protein chains: {gprotein_chains}")

    if receptor_chain_id not in [c.id for c in structure[0]]:
        print(f"  WARNING: Receptor chain {receptor_chain_id} not found in {pdb_id}")
        return None

    # Determine species
    species, human_offset = determine_species(structure, receptor_chain_id)
    print(f"  Species: {species.upper()} (human offset: +{human_offset})")

    # Collect all standard AA residues on receptor chain
    model = structure[0]
    receptor_residues = [
        res for res in model[receptor_chain_id]
        if is_aa(res, standard=True)
    ]

    # Classify shells
    shell_map = classify_shells(receptor_residues, ligand_atoms)

    # Calculate distances and build results
    results = []
    for residue in receptor_residues:
        res_num = residue.id[1]
        res_name = residue.get_resname()
        rid = (receptor_chain_id, res_num)

        min_dist = calculate_min_distance(residue, ligand_atoms)
        min_sc_dist = calculate_min_sidechain_distance(residue, ligand_atoms)

        if min_dist is not None:
            results.append({
                'pdb_id': pdb_id,
                'species': species,
                'drug_name': sdef['drug'],
                'receptor_chain': receptor_chain_id,
                'residue_number': res_num,
                'residue_number_human': res_num + human_offset,
                'residue_name': res_name,
                'ligand_name': ligand_label,
                'distance_angstrom': round(min_dist, 3),
                'distance_sidechain_angstrom': round(min_sc_dist, 3) if min_sc_dist else None,
                'shell': shell_map.get(rid, 'none'),
                'gprotein_chains': ','.join(gprotein_chains),
            })

    # Save individual CSV
    if results:
        outfile = os.path.join(output_dir, f"{pdb_id}_distances.csv")
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
    else:
        print(f"  WARNING: No distance measurements for {pdb_id}")

    return results


def generate_summary(all_results, output_dir):
    """Write combined CSV and text summary."""
    fieldnames = [
        'pdb_id', 'species', 'drug_name', 'receptor_chain',
        'residue_number', 'residue_number_human', 'residue_name',
        'ligand_name', 'distance_angstrom', 'distance_sidechain_angstrom',
        'shell', 'gprotein_chains'
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
        f.write("MOR Ligand Distance Analysis — ALL Experimental Structures\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Total structures analyzed: {len(all_results)}\n")
        total = sum(len(r) for r in all_results)
        f.write(f"Total distance measurements: {total:,}\n\n")

        f.write("Shell definitions:\n")
        f.write(f"  First shell:  sidechain atoms <= {FIRST_SHELL_CUTOFF} A from ligand\n")
        f.write(f"  Second shell: sidechain atoms <= {SECOND_SHELL_CUTOFF} A from first-shell sidechain\n")
        f.write(f"                (not already first shell)\n\n")

        f.write(f"Skipped structures:\n")
        f.write(f"  8qot — apo (antibody-bound, no orthosteric ligand)\n\n")

        # Drug summary
        drugs = {}
        for results in all_results:
            drug = results[0]['drug_name']
            drugs.setdefault(drug, []).append(results[0]['pdb_id'])
        f.write("Structures per drug:\n")
        for drug, pdbs in sorted(drugs.items()):
            f.write(f"  {drug}: {', '.join(pdbs)}\n")
        f.write("\n" + "-" * 80 + "\n")

        # Per-structure details
        for results in all_results:
            pdb_id = results[0]['pdb_id']
            species = results[0]['species']
            drug = results[0]['drug_name']
            ligand = results[0]['ligand_name']
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

            f.write(f"\n{pdb_id}:\n")
            f.write(f"  Drug: {drug}\n")
            f.write(f"  Species: {species.upper()}\n")
            f.write(f"  Ligand code: {ligand}\n")
            f.write(f"  G-protein chains: {gp}\n")
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
    struct_dir = str(CIF_DIR)
    output_dir = str(REPO_ROOT / "structures" / "processed" / "ligand_distances")
    os.makedirs(output_dir, exist_ok=True)

    print(f"Structure directory: {struct_dir}")
    print(f"Output directory:    {output_dir}")
    print(f"Structures to analyze: {len(STRUCTURE_DEFS)}")
    print(f"First shell cutoff:  {FIRST_SHELL_CUTOFF} A (sidechain -> ligand)")
    print(f"Second shell cutoff: {SECOND_SHELL_CUTOFF} A (sidechain -> first-shell sidechain)")

    all_results = []
    failed = []

    for pdb_id, sdef in STRUCTURE_DEFS.items():
        cif_file = os.path.join(struct_dir, f"{pdb_id}.cif")
        if not os.path.exists(cif_file):
            print(f"\n  WARNING: {cif_file} not found, skipping {pdb_id}")
            failed.append(pdb_id)
            continue

        results = analyze_structure(cif_file, sdef, output_dir)
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

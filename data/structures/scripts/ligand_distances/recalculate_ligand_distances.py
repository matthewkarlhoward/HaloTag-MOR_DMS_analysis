#!/usr/bin/env python3
"""
Recalculate ligand distances for specific MOR structures with correct ligand assignments.
Annotates residues as first shell (sidechain ≤5Å from ligand) or second shell
(sidechain ≤4.5Å from any first-shell residue's sidechain, but not first shell itself).

Structures and ligands:
  DAMGO:       8efq (chain P)
  Fentanyl:    8ef5 (7V7)
  Morphine:    8ef6 (MOI)
  SR17018:     8efl (WH9)
  C6-guano:    7u2k (KZR)
  PZM21:       8efo (8QY)
  MP:          7t2g (EIG)
  Lofentanil:  7t2h (EID)
  Oliceridine: 8efb (WH2)
  Naloxone:    9PXU (A1APV)
"""

import os
import numpy as np
import csv
from pathlib import Path

from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa

REPO_ROOT = Path(__file__).resolve().parents[3]
CIF_DIR = REPO_ROOT / "structures" / "raw" / "experimental"

# ── Structure definitions ──────────────────────────────────────────────────────

STRUCTURE_DEFS = {
    # DAMGO structure - peptide ligand spanning entire chain(s)
    # NOTE: BioPython uses auth_asym_id (PDB chain), not label_asym_id (mmCIF chain)
    # receptor_chain: the single chain closest to the ligand (verified empirically)
    "8efq": {"drug": "DAMGO", "type": "peptide", "ligand_chains": ["P"], "receptor_chain": "R"},
    # Small-molecule ligands
    "8ef5": {"drug": "Fentanyl",  "type": "small_mol", "ligand_code": "7V7", "receptor_chain": "R"},
    "8ef6": {"drug": "Morphine",  "type": "small_mol", "ligand_code": "MOI", "receptor_chain": "R"},
    "8efl": {"drug": "SR17018",   "type": "small_mol", "ligand_code": "WH9", "receptor_chain": "R"},
    "7u2k": {"drug": "C6-guano",  "type": "small_mol", "ligand_code": "KZR", "receptor_chain": "D"},
    "8efo": {"drug": "PZM21",     "type": "small_mol", "ligand_code": "8QY", "receptor_chain": "R"},
    "7t2g": {"drug": "Mitragynine_Pseudoindoxyl", "type": "small_mol", "ligand_code": "EIG", "receptor_chain": "R"},
    "7t2h": {"drug": "Lofentanil", "type": "small_mol", "ligand_code": "EID", "receptor_chain": "D"},
    "8efb": {"drug": "Oliceridine", "type": "small_mol", "ligand_code": "WH2", "receptor_chain": "R"},
    "9PXU": {"drug": "Naloxone",  "type": "small_mol", "ligand_code": "A1APV", "receptor_chain": "R"},
}

# Shell definition thresholds
FIRST_SHELL_CUTOFF = 5.0    # sidechain ≤5Å from ligand
SECOND_SHELL_CUTOFF = 4.5   # sidechain ≤4.5Å from first-shell sidechain


# ── Distance calculation functions ─────────────────────────────────────────────

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


# ── Species detection ──────────────────────────────────────────────────────────

def determine_species(structure, receptor_chain_ids):
    """Mouse MOR has ASP at position 147 (human = 149). Returns species and offset."""
    model = structure[0]
    has_asp147 = False
    has_asp149 = False
    for cid in receptor_chain_ids:
        for res in model[cid]:
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


# ── Chain identification ───────────────────────────────────────────────────────

def identify_receptor_and_gprotein_chains(structure, ligand_chain_ids=None):
    """Identify receptor and G-protein chains (excluding known ligand chains)."""
    model = structure[0]
    receptor_chains = []
    gprotein_chains = []
    exclude = set(ligand_chain_ids) if ligand_chain_ids else set()
    for chain in model:
        cid = chain.get_id()
        if cid in exclude:
            continue
        aa_count = sum(1 for r in chain if is_aa(r, standard=True))
        if 200 < aa_count < 500:
            receptor_chains.append(cid)
        elif aa_count > 50:
            gprotein_chains.append(cid)
    return receptor_chains, gprotein_chains


# ── Ligand atom collection ─────────────────────────────────────────────────────

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
    """Collect ALL atoms from specified chains (for peptide ligands like DAMGO)."""
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


# ── Shell classification ──────────────────────────────────────────────────────

def classify_shells(receptor_residues, ligand_atoms):
    """
    Classify receptor residues into binding shells.

    First shell:  sidechain atoms ≤ 5.0 Å from any ligand atom
    Second shell: NOT first shell, but sidechain atoms ≤ 4.5 Å from any
                  first-shell residue's sidechain atoms

    Parameters
    ----------
    receptor_residues : list of Bio.PDB.Residue
        All standard AA residues on the receptor chain
    ligand_atoms : list of Bio.PDB.Atom
        All atoms belonging to the ligand

    Returns
    -------
    dict : {residue_id: 'first_shell', 'second_shell', or 'none'}
    """
    # Pass 1: identify first-shell residues
    first_shell_residues = []
    residue_shell = {}

    for res in receptor_residues:
        sc_dist = calculate_min_sidechain_distance(res, ligand_atoms)
        rid = (res.get_parent().id, res.id[1])  # (chain, resnum)
        if sc_dist is not None and sc_dist <= FIRST_SHELL_CUTOFF:
            residue_shell[rid] = 'first_shell'
            first_shell_residues.append(res)
        else:
            residue_shell[rid] = 'none'

    print(f"  First shell ({FIRST_SHELL_CUTOFF}Å sidechain→ligand): {len(first_shell_residues)} residues")

    # Pass 2: identify second-shell residues
    # For each non-first-shell residue, check sidechain-sidechain distance
    # to every first-shell residue
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
                break  # only need one contact

    print(f"  Second shell ({SECOND_SHELL_CUTOFF}Å sidechain→first-shell sidechain): {second_shell_count} residues")

    return residue_shell


# ── Main analysis ──────────────────────────────────────────────────────────────

def analyze_structure(cif_file, sdef, output_dir):
    """Analyze one structure with specified ligand definition."""
    parser = MMCIFParser(QUIET=True)
    pdb_id = Path(cif_file).stem

    print(f"\n{'='*70}")
    print(f"Analyzing {pdb_id}  —  {sdef['drug']}")
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
        ligand_label = "DAMGO"
    else:
        ligand_chain_ids = None
        ligand_atoms = collect_ligand_atoms_small_mol(structure, sdef['ligand_code'])
        ligand_label = sdef['ligand_code']

    if not ligand_atoms:
        return None

    # Use specified receptor chain
    receptor_chain_id = sdef['receptor_chain']

    # Identify other chains for reporting
    exclude_chains = sdef.get('ligand_chains', [])
    all_receptor_chains, gprotein_chains = identify_receptor_and_gprotein_chains(
        structure, ligand_chain_ids=exclude_chains
    )
    print(f"  Using receptor chain: {receptor_chain_id}  (all candidates: {all_receptor_chains})")
    print(f"  G-protein chains: {gprotein_chains}")

    if receptor_chain_id not in [c.id for c in structure[0]]:
        print(f"  WARNING: Receptor chain {receptor_chain_id} not found in {pdb_id}")
        return None

    # Determine species
    species, human_offset = determine_species(structure, [receptor_chain_id])
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
        print(f"  Wrote {len(results)} measurements → {outfile}")
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
        f.write("MOR Ligand Distance Analysis — Corrected Ligand Assignments\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Total structures analyzed: {len(all_results)}\n")
        total = sum(len(r) for r in all_results)
        f.write(f"Total distance measurements: {total:,}\n\n")

        f.write("Shell definitions:\n")
        f.write(f"  First shell:  sidechain atoms ≤ {FIRST_SHELL_CUTOFF} Å from ligand\n")
        f.write(f"  Second shell: sidechain atoms ≤ {SECOND_SHELL_CUTOFF} Å from first-shell sidechain\n")
        f.write(f"                (not already first shell)\n\n")

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
            f.write(f"  Distance range: {min(distances):.2f} – {max(distances):.2f} Å\n")
            f.write(f"  Mean distance: {np.mean(distances):.2f} Å\n")
            f.write(f"  First shell:  {len(first_shell)} residues\n")
            f.write(f"  Second shell: {len(second_shell)} residues\n")

            if first_shell:
                res_str = ", ".join(
                    f"{r['residue_name']}{r['residue_number_human']}({r['distance_sidechain_angstrom']:.1f}Å)"
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
    output_dir = str(REPO_ROOT / "structures" / "processed" / "ligand_distances_tested_drugs")
    os.makedirs(output_dir, exist_ok=True)

    print(f"Structure directory: {struct_dir}")
    print(f"Output directory:    {output_dir}")
    print(f"Structures to analyze: {len(STRUCTURE_DEFS)}")
    print(f"First shell cutoff:  {FIRST_SHELL_CUTOFF} Å (sidechain → ligand)")
    print(f"Second shell cutoff: {SECOND_SHELL_CUTOFF} Å (sidechain → first-shell sidechain)")

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

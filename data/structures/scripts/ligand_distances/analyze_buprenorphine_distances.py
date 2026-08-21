#!/usr/bin/env python3
"""
Ligand distance analysis for the buprenorphine PDB structure.
Applies the same shell-classification logic as recalculate_ligand_distances.py
but reads a PDB file instead of mmCIF.

Structure:
  Buprenorphine: buprenorphine.pdb
    receptor chain : R
    ligand residue : LIG (residue 701, HETATM, chain R)
    species        : human (ASP149 present → offset = 0)

Shell definitions (identical to the CIF workflow):
  First shell  : sidechain ≤ 5.0 Å from any ligand atom
  Second shell : NOT first shell AND sidechain ≤ 4.5 Å from any
                 first-shell residue's sidechain atoms
"""

import os
import csv
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
CIF_DIR = REPO_ROOT / "structures" / "raw" / "experimental"

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

# ── Constants ─────────────────────────────────────────────────────────────────

PDB_FILE      = str(CIF_DIR / "buprenorphine.pdb")
RECEPTOR_CHAIN = "R"
LIGAND_CODE    = "LIG"
DRUG_NAME      = "Buprenorphine"
OUTPUT_DIR     = str(REPO_ROOT / "structures" / "processed" / "ligand_distances_tested_drugs")

FIRST_SHELL_CUTOFF  = 5.0   # sidechain → ligand (Å)
SECOND_SHELL_CUTOFF = 4.5   # sidechain → first-shell sidechain (Å)

BACKBONE_ATOMS = {'N', 'CA', 'C', 'O'}


# ── Atom-selection helpers ────────────────────────────────────────────────────

def get_sidechain_atoms(residue):
    """Sidechain atoms; for GLY return all atoms."""
    if residue.get_resname() == 'GLY':
        return list(residue.get_atoms())
    return [a for a in residue.get_atoms() if a.get_name() not in BACKBONE_ATOMS]


def min_dist(atoms_a, atoms_b):
    """Minimum distance between two atom lists; returns inf if either is empty."""
    if not atoms_a or not atoms_b:
        return float('inf')
    best = float('inf')
    for a in atoms_a:
        for b in atoms_b:
            d = a - b
            if d < best:
                best = d
    return best


def calculate_min_distance(residue, ligand_atoms):
    """Minimum distance between any residue atom and any ligand atom."""
    return min_dist(list(residue.get_atoms()), ligand_atoms)


def calculate_min_sidechain_distance(residue, ligand_atoms):
    """Minimum sidechain–ligand distance."""
    return min_dist(get_sidechain_atoms(residue), ligand_atoms)


def calculate_min_sidechain_sidechain_distance(res_a, res_b):
    """Minimum sidechain–sidechain distance between two residues."""
    return min_dist(get_sidechain_atoms(res_a), get_sidechain_atoms(res_b))


# ── Species detection ─────────────────────────────────────────────────────────

def determine_species(model, receptor_chain_id):
    """
    Mouse MOR: ASP147 (human MOR: ASP149).
    Returns (species_str, human_offset).
    """
    chain = model[receptor_chain_id]
    has_asp147 = any(r.get_resname() == 'ASP' and r.id[1] == 147 for r in chain)
    has_asp149 = any(r.get_resname() == 'ASP' and r.id[1] == 149 for r in chain)
    if has_asp147 and not has_asp149:
        return 'mouse', 2
    elif has_asp149:
        return 'human', 0
    return 'unknown', 0


# ── Shell classification ───────────────────────────────────────────────────────

def classify_shells(receptor_residues, ligand_atoms):
    """
    Returns dict {(chain_id, res_num): 'first_shell'|'second_shell'|'none'}.
    """
    first_shell_residues = []
    residue_shell = {}

    # Pass 1 – first shell
    for res in receptor_residues:
        rid = (res.get_parent().id, res.id[1])
        sc_dist = calculate_min_sidechain_distance(res, ligand_atoms)
        if sc_dist <= FIRST_SHELL_CUTOFF:
            residue_shell[rid] = 'first_shell'
            first_shell_residues.append(res)
        else:
            residue_shell[rid] = 'none'

    print(f"  First shell ({FIRST_SHELL_CUTOFF} Å sidechain→ligand): "
          f"{len(first_shell_residues)} residues")

    # Pass 2 – second shell
    second_count = 0
    for res in receptor_residues:
        rid = (res.get_parent().id, res.id[1])
        if residue_shell[rid] == 'first_shell':
            continue
        for fs_res in first_shell_residues:
            if calculate_min_sidechain_sidechain_distance(res, fs_res) <= SECOND_SHELL_CUTOFF:
                residue_shell[rid] = 'second_shell'
                second_count += 1
                break

    print(f"  Second shell ({SECOND_SHELL_CUTOFF} Å sidechain→first-shell sidechain): "
          f"{second_count} residues")

    return residue_shell


# ── Ligand atom collection ─────────────────────────────────────────────────────

def collect_ligand_atoms(model, ligand_code):
    """
    Collect HETATM atoms whose residue name matches ligand_code.
    In BioPython PDB parsing, HETATM residues have id[0] == 'H_<resname>'
    or simply a het flag. We search all chains for a matching resname.
    """
    atoms = []
    for chain in model:
        for residue in chain:
            if residue.get_resname().strip() == ligand_code:
                atoms.extend(list(residue.get_atoms()))
    if not atoms:
        print(f"  WARNING: No atoms found for ligand '{ligand_code}'")
    else:
        print(f"  Found {len(atoms)} atoms for ligand '{ligand_code}'")
    return atoms


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    pdb_id = Path(PDB_FILE).stem   # 'buprenorphine'

    print(f"\n{'='*70}")
    print(f"Analyzing {pdb_id}  —  {DRUG_NAME}")
    print(f"{'='*70}")
    print(f"  PDB file:       {PDB_FILE}")
    print(f"  Receptor chain: {RECEPTOR_CHAIN}")
    print(f"  Ligand code:    {LIGAND_CODE}")
    print(f"  First shell:    ≤ {FIRST_SHELL_CUTOFF} Å (sidechain → ligand)")
    print(f"  Second shell:   ≤ {SECOND_SHELL_CUTOFF} Å (sidechain → first-shell sidechain)")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, PDB_FILE)
    model = structure[0]

    # Collect ligand atoms
    ligand_atoms = collect_ligand_atoms(model, LIGAND_CODE)
    if not ligand_atoms:
        print("ERROR: Could not find ligand atoms. Aborting.")
        return

    # Species
    species, human_offset = determine_species(model, RECEPTOR_CHAIN)
    print(f"  Species: {species.upper()} (human_offset: +{human_offset})")

    # Receptor residues (standard amino acids on chain R)
    receptor_residues = [
        res for res in model[RECEPTOR_CHAIN]
        if is_aa(res, standard=True)
    ]
    print(f"  Receptor residues (standard AA): {len(receptor_residues)}")

    # Shell classification
    shell_map = classify_shells(receptor_residues, ligand_atoms)

    # Build results
    results = []
    for residue in receptor_residues:
        res_num  = residue.id[1]
        res_name = residue.get_resname()
        rid      = (RECEPTOR_CHAIN, res_num)

        all_dist = calculate_min_distance(residue, ligand_atoms)
        sc_dist  = calculate_min_sidechain_distance(residue, ligand_atoms)

        if all_dist == float('inf'):
            continue

        results.append({
            'pdb_id':                      pdb_id,
            'species':                     species,
            'drug_name':                   DRUG_NAME,
            'receptor_chain':              RECEPTOR_CHAIN,
            'residue_number':              res_num,
            'residue_number_human':        res_num + human_offset,
            'residue_name':                res_name,
            'ligand_name':                 LIGAND_CODE,
            'distance_angstrom':           round(all_dist, 3),
            'distance_sidechain_angstrom': round(sc_dist, 3) if sc_dist != float('inf') else None,
            'shell':                       shell_map.get(rid, 'none'),
            'gprotein_chains':             '',   # none present in this structure
        })

    # Save CSV
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    outfile = os.path.join(OUTPUT_DIR, f"{pdb_id}_distances.csv")
    fieldnames = [
        'pdb_id', 'species', 'drug_name', 'receptor_chain',
        'residue_number', 'residue_number_human', 'residue_name',
        'ligand_name', 'distance_angstrom', 'distance_sidechain_angstrom',
        'shell', 'gprotein_chains',
    ]
    with open(outfile, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    n_first  = sum(1 for r in results if r['shell'] == 'first_shell')
    n_second = sum(1 for r in results if r['shell'] == 'second_shell')
    n_none   = len(results) - n_first - n_second

    print(f"\n  Wrote {len(results)} measurements → {outfile}")
    print(f"  Shell summary: {n_first} first, {n_second} second, {n_none} none")

    # Print first-shell residues sorted by sidechain distance
    first_shell = sorted(
        [r for r in results if r['shell'] == 'first_shell'],
        key=lambda x: x['distance_sidechain_angstrom'] or 999,
    )
    if first_shell:
        print("\n  First-shell residues (sorted by sidechain distance):")
        for r in first_shell:
            print(f"    {r['residue_name']}{r['residue_number_human']:>4}  "
                  f"sc={r['distance_sidechain_angstrom']:.2f} Å  "
                  f"all={r['distance_angstrom']:.2f} Å")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()

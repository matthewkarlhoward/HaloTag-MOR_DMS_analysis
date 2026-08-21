#!/usr/bin/env python3
"""
Analyze MOR structures: G-alpha protein distances + shell annotation.

Extended version of analyze_gpcr_gprotein_distances.py that additionally
records WHICH G-alpha residue is closest to each receptor residue.

Adds columns:
  galpha_residue_number       — sequence number of the closest G-alpha residue
  galpha_residue_name         — amino acid name of the closest G-alpha residue
  galpha_residue_number_sc    — same, but based on sidechain-to-all distance
  galpha_residue_name_sc      — same
"""

import os
import numpy as np
import csv
from pathlib import Path

from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

REPO_ROOT = Path(__file__).resolve().parents[3]
CIF_DIR = REPO_ROOT / "structures" / "raw" / "experimental"
OUTPUT_DIR = Path("/Users/mkh/GitHub/mor_dms_analysis/structures/distances/gprotein")

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

FIRST_SHELL_CUTOFF  = 5.0
SECOND_SHELL_CUTOFF = 4.5
BACKBONE_ATOMS = {'N', 'CA', 'C', 'O'}


def get_sidechain_atoms(residue):
    if residue.get_resname() == 'GLY':
        return list(residue.get_atoms())
    return [a for a in residue.get_atoms() if a.get_name() not in BACKBONE_ATOMS]


def find_closest_galpha_residue(receptor_res, galpha_residues):
    """
    For a receptor residue, find:
      - minimum all-atom distance to any G-alpha atom, and which G-alpha residue provides it
      - minimum sidechain distance to any G-alpha atom, and which G-alpha residue provides it

    Returns:
      (min_allatom_dist, closest_galpha_res_allatom,
       min_sc_dist,      closest_galpha_res_sc)
    """
    all_atoms = list(receptor_res.get_atoms())
    sc_atoms  = get_sidechain_atoms(receptor_res)

    min_all_dist = float('inf')
    min_sc_dist  = float('inf')
    closest_all  = None
    closest_sc   = None

    for galpha_res in galpha_residues:
        ga_atoms = list(galpha_res.get_atoms())
        if not ga_atoms:
            continue

        # All-atom
        for ra in all_atoms:
            for ga in ga_atoms:
                d = ra - ga
                if d < min_all_dist:
                    min_all_dist = d
                    closest_all  = galpha_res

        # Sidechain
        for sa in sc_atoms:
            for ga in ga_atoms:
                d = sa - ga
                if d < min_sc_dist:
                    min_sc_dist = d
                    closest_sc  = galpha_res

    return (
        min_all_dist if min_all_dist < float('inf') else None, closest_all,
        min_sc_dist  if min_sc_dist  < float('inf') else None, closest_sc,
    )


def classify_shells(receptor_residues, galpha_residues):
    """Same shell logic as original, using per-residue iteration."""
    galpha_atoms = [a for res in galpha_residues for a in res.get_atoms()]

    first_shell_residues = []
    residue_shell = {}

    for res in receptor_residues:
        sc_atoms = get_sidechain_atoms(res)
        rid = (res.get_parent().id, res.id[1])
        min_sc = min((sa - ga for sa in sc_atoms for ga in galpha_atoms),
                     default=float('inf'))
        if min_sc <= FIRST_SHELL_CUTOFF:
            residue_shell[rid] = 'first_shell'
            first_shell_residues.append(res)
        else:
            residue_shell[rid] = 'none'

    for res in receptor_residues:
        rid = (res.get_parent().id, res.id[1])
        if residue_shell[rid] == 'first_shell':
            continue
        sc_a = get_sidechain_atoms(res)
        for fs_res in first_shell_residues:
            sc_b = get_sidechain_atoms(fs_res)
            d = min((a - b for a in sc_a for b in sc_b), default=float('inf'))
            if d <= SECOND_SHELL_CUTOFF:
                residue_shell[rid] = 'second_shell'
                break

    return residue_shell


def analyze_structure(cif_file, pdb_id, sdef):
    parser = MMCIFParser(QUIET=True)
    print(f"\n{'='*60}\nAnalyzing {pdb_id}\n{'='*60}")

    try:
        structure = parser.get_structure(pdb_id, cif_file)
    except Exception as e:
        print(f"  ERROR: {e}")
        return None

    model = structure[0]
    species      = sdef["species"]
    human_offset = 2 if species == "mouse" else 0

    try:
        receptor_chain = model[sdef["receptor_chain"]]
        galpha_chain   = model[sdef["galpha_chain"]]
    except KeyError as e:
        print(f"  WARNING: chain not found: {e}")
        return None

    galpha_residues  = [r for r in galpha_chain   if is_aa(r, standard=True)]
    receptor_residues = [r for r in receptor_chain if is_aa(r, standard=True)]

    print(f"  Receptor: {len(receptor_residues)} residues | Gα: {len(galpha_residues)} residues")

    shell_map = classify_shells(receptor_residues, galpha_residues)

    results = []
    for res in receptor_residues:
        res_num  = res.id[1]
        res_name = res.get_resname()
        rid      = (sdef["receptor_chain"], res_num)

        min_all, cl_all, min_sc, cl_sc = find_closest_galpha_residue(res, galpha_residues)

        if min_all is None:
            continue

        results.append({
            'pdb_id':                    pdb_id,
            'species':                   species,
            'receptor_chain':            sdef["receptor_chain"],
            'residue_number':            res_num,
            'residue_number_human':      res_num + human_offset,
            'residue_name':              res_name,
            'galpha_chain':              sdef["galpha_chain"],
            'distance_angstrom':         round(min_all, 3),
            'distance_sidechain_angstrom': round(min_sc, 3) if min_sc is not None else None,
            'shell':                     shell_map.get(rid, 'none'),
            'galpha_residue_number':     cl_all.id[1]        if cl_all else None,
            'galpha_residue_name':       cl_all.get_resname() if cl_all else None,
            'galpha_residue_number_sc':  cl_sc.id[1]         if cl_sc  else None,
            'galpha_residue_name_sc':    cl_sc.get_resname()  if cl_sc  else None,
        })

    return results


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        'pdb_id', 'species', 'receptor_chain',
        'residue_number', 'residue_number_human', 'residue_name',
        'galpha_chain',
        'distance_angstrom', 'distance_sidechain_angstrom', 'shell',
        'galpha_residue_number', 'galpha_residue_name',
        'galpha_residue_number_sc', 'galpha_residue_name_sc',
    ]

    all_results = []
    failed = []

    for pdb_id, sdef in sorted(STRUCTURE_DEFS.items(), key=lambda x: x[0].lower()):
        cif_path = CIF_DIR / f"{pdb_id}.cif"
        if not cif_path.exists():
            print(f"  WARNING: {cif_path} not found, skipping")
            failed.append(pdb_id)
            continue

        results = analyze_structure(str(cif_path), pdb_id, sdef)
        if not results:
            failed.append(pdb_id)
            continue

        # Per-structure CSV
        out_csv = OUTPUT_DIR / f"{pdb_id}_gprotein_distances.csv"
        with open(out_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)

        n1 = sum(1 for r in results if r['shell'] == 'first_shell')
        n2 = sum(1 for r in results if r['shell'] == 'second_shell')
        print(f"  -> {out_csv.name}  ({n1} first shell, {n2} second shell)")
        all_results.extend(results)

    # Combined CSV
    combined = OUTPUT_DIR / "all_structures_combined.csv"
    with open(combined, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)

    print(f"\nDONE: {len(STRUCTURE_DEFS) - len(failed)} structures | combined -> {combined}")
    if failed:
        print(f"Failed/skipped: {', '.join(failed)}")


if __name__ == "__main__":
    main()

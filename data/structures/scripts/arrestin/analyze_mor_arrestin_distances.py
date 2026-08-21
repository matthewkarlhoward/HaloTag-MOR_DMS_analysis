#!/usr/bin/env python3
"""
MOR - beta-arrestin distance / shell analysis for 9WSV and 9WSX.

Mirrors analyze_gpcr_gprotein_distances_with_galpha_residue.py but for
MOR-arrestin complexes. Chain R = receptor (MOR), chain C = beta-arrestin 1.
"""

import csv
from pathlib import Path

from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

REPO_ROOT = Path(__file__).resolve().parents[3]
CIF_DIR = REPO_ROOT / "structures" / "raw" / "experimental"
OUTPUT_DIR = Path(
    "/Users/mkh/GitHub/mor_dms_analysis/structures/distances/arrestin"
)

STRUCTURE_DEFS = {
    "9WSV": {"receptor_chain": "R", "arrestin_chain": "C", "species": "mouse",
             "arrestin_isoform": "barr1"},
    "9WSX": {"receptor_chain": "R", "arrestin_chain": "C", "species": "mouse",
             "arrestin_isoform": "barr1"},
}

FIRST_SHELL_CUTOFF = 5.0
SECOND_SHELL_CUTOFF = 4.5
BACKBONE_ATOMS = {"N", "CA", "C", "O"}


def get_sidechain_atoms(residue):
    if residue.get_resname() == "GLY":
        return list(residue.get_atoms())
    return [a for a in residue.get_atoms() if a.get_name() not in BACKBONE_ATOMS]


def find_closest_arrestin_residue(receptor_res, arr_residues):
    all_atoms = list(receptor_res.get_atoms())
    sc_atoms = get_sidechain_atoms(receptor_res)

    min_all = float("inf")
    min_sc = float("inf")
    closest_all = None
    closest_sc = None

    for arr in arr_residues:
        a_atoms = list(arr.get_atoms())
        if not a_atoms:
            continue
        for ra in all_atoms:
            for aa in a_atoms:
                d = ra - aa
                if d < min_all:
                    min_all = d
                    closest_all = arr
        for sa in sc_atoms:
            for aa in a_atoms:
                d = sa - aa
                if d < min_sc:
                    min_sc = d
                    closest_sc = arr

    return (
        min_all if min_all < float("inf") else None, closest_all,
        min_sc if min_sc < float("inf") else None, closest_sc,
    )


def classify_shells(receptor_residues, arr_residues):
    arr_atoms = [a for res in arr_residues for a in res.get_atoms()]

    first_shell = []
    shell = {}

    for res in receptor_residues:
        sc = get_sidechain_atoms(res)
        rid = (res.get_parent().id, res.id[1])
        min_sc = min((s - a for s in sc for a in arr_atoms), default=float("inf"))
        if min_sc <= FIRST_SHELL_CUTOFF:
            shell[rid] = "first_shell"
            first_shell.append(res)
        else:
            shell[rid] = "none"

    for res in receptor_residues:
        rid = (res.get_parent().id, res.id[1])
        if shell[rid] == "first_shell":
            continue
        sc_a = get_sidechain_atoms(res)
        for fs in first_shell:
            sc_b = get_sidechain_atoms(fs)
            d = min((a - b for a in sc_a for b in sc_b), default=float("inf"))
            if d <= SECOND_SHELL_CUTOFF:
                shell[rid] = "second_shell"
                break

    return shell


def analyze(cif_file, pdb_id, sdef):
    parser = MMCIFParser(QUIET=True)
    print(f"\n{'='*60}\nAnalyzing {pdb_id}\n{'='*60}")
    try:
        s = parser.get_structure(pdb_id, cif_file)
    except Exception as e:
        print(f"  ERROR: {e}")
        return None

    model = s[0]
    species = sdef["species"]
    human_offset = 2 if species == "mouse" else 0

    try:
        receptor_chain = model[sdef["receptor_chain"]]
        arr_chain = model[sdef["arrestin_chain"]]
    except KeyError as e:
        print(f"  WARNING: chain not found: {e}")
        return None

    arr_residues = [r for r in arr_chain if is_aa(r, standard=True)]
    rec_residues = [r for r in receptor_chain if is_aa(r, standard=True)]
    print(f"  Receptor: {len(rec_residues)} | Arrestin: {len(arr_residues)}")

    shell_map = classify_shells(rec_residues, arr_residues)

    out = []
    for res in rec_residues:
        res_num = res.id[1]
        rid = (sdef["receptor_chain"], res_num)
        min_all, cl_all, min_sc, cl_sc = find_closest_arrestin_residue(res, arr_residues)
        if min_all is None:
            continue
        out.append({
            "pdb_id": pdb_id,
            "species": species,
            "arrestin_isoform": sdef["arrestin_isoform"],
            "receptor_chain": sdef["receptor_chain"],
            "residue_number": res_num,
            "residue_number_human": res_num + human_offset,
            "residue_name": res.get_resname(),
            "arrestin_chain": sdef["arrestin_chain"],
            "distance_angstrom": round(min_all, 3),
            "distance_sidechain_angstrom": round(min_sc, 3) if min_sc is not None else None,
            "shell": shell_map.get(rid, "none"),
            "arrestin_residue_number": cl_all.id[1] if cl_all else None,
            "arrestin_residue_name": cl_all.get_resname() if cl_all else None,
            "arrestin_residue_number_sc": cl_sc.id[1] if cl_sc else None,
            "arrestin_residue_name_sc": cl_sc.get_resname() if cl_sc else None,
        })
    return out


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "pdb_id", "species", "arrestin_isoform", "receptor_chain",
        "residue_number", "residue_number_human", "residue_name",
        "arrestin_chain",
        "distance_angstrom", "distance_sidechain_angstrom", "shell",
        "arrestin_residue_number", "arrestin_residue_name",
        "arrestin_residue_number_sc", "arrestin_residue_name_sc",
    ]

    all_rows = []
    for pdb_id, sdef in sorted(STRUCTURE_DEFS.items()):
        cif = CIF_DIR / f"{pdb_id}.cif"
        if not cif.exists():
            print(f"  WARNING: {cif} not found")
            continue
        rows = analyze(str(cif), pdb_id, sdef)
        if not rows:
            continue
        out_csv = OUTPUT_DIR / f"{pdb_id}_arrestin_distances.csv"
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)
        n1 = sum(1 for r in rows if r["shell"] == "first_shell")
        n2 = sum(1 for r in rows if r["shell"] == "second_shell")
        print(f"  -> {out_csv.name} ({n1} first shell, {n2} second shell)")
        all_rows.extend(rows)

    combined = OUTPUT_DIR / "all_structures_combined.csv"
    with open(combined, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(all_rows)
    print(f"\nDONE: combined -> {combined}")


if __name__ == "__main__":
    main()

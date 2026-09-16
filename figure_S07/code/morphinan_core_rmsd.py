#!/usr/bin/env python3
"""
Morphinan-core RMSD of the Chai-1 predicted poses to the experimental naloxone pose.

Answers: "for the morphinan cores, the RMSD to naloxone is *** A".

Method
------
1. Superpose each Chai-1 model's receptor onto the naloxone reference receptor
   (9PXU) by Kabsch on shared CA atoms, and apply that transform to the ligand.
   This is the same receptor-frame superposition used for the overlay panel in
   supplemental figure 7b; RMSD is therefore measured in the pocket frame, not
   after independently best-fitting the two ligands, which would flatter the
   agreement.
2. Perceive each ligand from its coordinates and find the maximum common
   substructure with naloxone (rdkit FindMCS, ring-matching enforced). For the
   4,5-epoxymorphinans this recovers the morphinan core; for butorphanol, which
   lacks the epoxy bridge, it recovers the shared phenanthrene core.
3. RMSD over the matched core atoms, minimised over symmetry-equivalent
   matches only (no further superposition).

Usage:  python3 morphinan_core_rmsd.py
"""
from pathlib import Path
import itertools
import re
import sys

import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import rdFMCS, rdDetermineBonds
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError:
    sys.exit("rdkit is required:  pip install rdkit")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
STRUCT = ROOT / "data" / "structures"
# Two naloxone references. 9PXU is the inactive, Nb6-bound structure used for
# the ligand-contact analysis; 9PXY is the active, Gi1-bound one. The Chai-1
# models are Gi-bound and active-like, so both are reported.
# (label, path, receptor chain, naloxone component id)
REFS = [("9PXU (inactive, Nb6)", STRUCT / "raw" / "experimental" / "9PXU.cif", "R", "A1APV"),
        ("9PXY (active, Gi1)",   STRUCT / "raw" / "experimental" / "9PXY.cif", "R", "A1APV")]

# Chai-1 models, and the rank used for contact assignment in this study.
MODELS = [
    ("Naltrexone",  "MOR-Naltrexone",       "pred.rank_1.cif"),
    ("Nalbuphine",  "MOR-GNAI-Nalbuphine",  "pred.rank_0_nalbuphine.cif"),
    ("Butorphanol", "MOR-GNAI-Butorphanol", "pred.rank_0_butorphanol.cif"),
    ("Methadone",   "MOR-GNAI-Methadone",   "pred.rank_0_methadone.cif"),
]

# Superposition is restricted to the structured receptor (TM1-H8). The Chai-1
# models are full length; 9PXU is a construct with disordered termini, so
# including them lets flexible tails dominate the fit.
POS_MIN, POS_MAX = 65, 355


# ---------------------------------------------------------------- parsing ---
def parse_cif_atoms(path):
    """Minimal mmCIF atom_site reader -> list of dicts. Handles the column
    order of both RCSB entries and Chai-1 output by reading the loop header."""
    cols, rows, in_loop = [], [], False
    for line in open(path, errors="ignore"):
        s = line.strip()
        if s.startswith("_atom_site."):
            cols.append(s.split(".", 1)[1]); in_loop = True; continue
        if in_loop and (s.startswith(("ATOM", "HETATM"))):
            f = s.split()
            # Chai-1 omits the trailing pdbx_PDB_model_num on some rows; every
            # column up to B_iso is still aligned, so zip what is present.
            if len(f) >= 13:
                rows.append(dict(zip(cols, f)))
        elif in_loop and s.startswith("#"):
            in_loop = False
    return cols, rows


def get(d, *names, default=None):
    for n in names:
        if n in d:
            return d[n]
    return default


def receptor_ca(rows, chain=None):
    """{residue number: xyz} for protein CA atoms. If chain is given, use it;
    otherwise take the longest polymer chain."""
    from collections import defaultdict
    per_chain = defaultdict(dict)
    for r in rows:
        if r.get("group_PDB") != "ATOM":
            continue
        if get(r, "label_atom_id", "auth_atom_id") != "CA":
            continue
        ch = get(r, "auth_asym_id", "label_asym_id")
        num = get(r, "auth_seq_id", "label_seq_id")
        try:
            num = int(num)
        except (TypeError, ValueError):
            continue
        xyz = [float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])]
        per_chain[ch][num] = xyz
    if not per_chain:
        return {}
    if chain is not None and chain in per_chain:
        return per_chain[chain]
    return max(per_chain.values(), key=len)


def ligand_atoms(rows, comp=None):
    """Heavy-atom (element, xyz) for a HETATM ligand. If comp is None, take the
    largest non-solvent HETATM component."""
    from collections import defaultdict
    groups = defaultdict(list)
    SKIP = {"HOH", "NA", "CL", "ZN", "MG", "SO4", "PO4", "GOL", "EDO", "CLR",
            "OLA", "OLC", "PEG", "ACT", "GDP", "GTP", "K"}
    for r in rows:
        if r.get("group_PDB") != "HETATM":
            continue
        cid = get(r, "label_comp_id", "auth_comp_id")
        if cid in SKIP:
            continue
        el = get(r, "type_symbol", default="C")
        if el == "H":
            continue
        ch = get(r, "auth_asym_id", "label_asym_id")
        key = (cid, ch)
        groups[key].append((el, [float(r["Cartn_x"]), float(r["Cartn_y"]),
                                float(r["Cartn_z"])]))
    if comp is not None:
        for (cid, ch), v in groups.items():
            if cid == comp:
                return v
    return max(groups.values(), key=len) if groups else []


# ------------------------------------------------------------- geometry ---
def kabsch(P, Q):
    """Rotation+translation mapping P onto Q (both N x 3)."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, Qc - R @ Pc


def to_mol(atoms):
    """Build an rdkit Mol with bonds perceived from coordinates."""
    m = Chem.RWMol()
    conf = Chem.Conformer(len(atoms))
    for i, (el, xyz) in enumerate(atoms):
        a = Chem.Atom(el.capitalize())
        a.SetNoImplicit(True)
        m.AddAtom(a)
        conf.SetAtomPosition(i, tuple(float(c) for c in xyz))
    mol = m.GetMol()
    mol.AddConformer(conf)
    try:
        rdDetermineBonds.DetermineConnectivity(mol)
    except Exception:
        return None
    # ring perception must be initialised before substructure matching
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_SYMMRINGS |
                          Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
                     catchErrors=True)
    Chem.FastFindRings(mol)
    return mol


def core_rmsd(mol_a, mol_b):
    """MCS between two coordinate-perceived mols; RMSD over matched atoms,
    minimised across symmetry-equivalent mappings. No re-superposition."""
    res = rdFMCS.FindMCS([mol_a, mol_b], ringMatchesRingOnly=True,
                         completeRingsOnly=True, timeout=30,
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         bondCompare=rdFMCS.BondCompare.CompareAny)
    if res.canceled or res.numAtoms < 6:
        return None, 0
    patt = Chem.MolFromSmarts(res.smartsString)
    ma = mol_a.GetSubstructMatches(patt, uniquify=False, maxMatches=2000)
    mb = mol_b.GetSubstructMatches(patt, uniquify=False, maxMatches=2000)
    if not ma or not mb:
        return None, 0
    ca = mol_a.GetConformer().GetPositions()
    cb = mol_b.GetConformer().GetPositions()
    best = None
    for A in ma:
        for B in mb:
            d = ca[list(A)] - cb[list(B)]
            r = float(np.sqrt((d ** 2).sum(1).mean()))
            best = r if best is None or r < best else best
    return best, res.numAtoms


# ------------------------------------------------------------------ main ---
ALL = []


def main():
  for ref_label, ref_path, ref_chain, ref_comp in REFS:
    _, ref_rows = parse_cif_atoms(ref_path)
    ref_ca = receptor_ca(ref_rows, ref_chain)
    ref_lig = ligand_atoms(ref_rows, ref_comp)
    print(f"\n=== reference: naloxone in {ref_label} — "
          f"receptor CA {len(ref_ca)}, ligand heavy atoms {len(ref_lig)} ===")
    ref_mol = to_mol(ref_lig)

    print(f"{'ligand':<14}{'model':<34}{'CA used':>8}{'CA RMSD':>9}"
          f"{'core':>6}{'core RMSD':>11}")
    out = []
    for name, folder, fn in MODELS:
        path = STRUCT / "raw" / "chai_predicted" / folder / fn
        if not path.exists():
            print(f"  {name}: missing {path}")
            continue
        _, rows = parse_cif_atoms(path)
        ca = receptor_ca(rows)
        shared = sorted(i for i in (set(ca) & set(ref_ca))
                        if POS_MIN <= i <= POS_MAX)
        if len(shared) < 50:
            print(f"  {name}: only {len(shared)} shared CA, skipping")
            continue
        P = np.array([ca[i] for i in shared])
        Q = np.array([ref_ca[i] for i in shared])
        R, t = kabsch(P, Q)
        ca_rmsd = float(np.sqrt((((P @ R.T) + t - Q) ** 2).sum(1).mean()))

        lig = ligand_atoms(rows)
        moved = [(el, list((R @ np.array(xyz)) + t)) for el, xyz in lig]
        mol = to_mol(moved)
        if mol is None or ref_mol is None:
            print(f"  {name}: bond perception failed")
            continue
        rms, n = core_rmsd(mol, ref_mol)
        print(f"{name:<14}{folder+'/'+fn.replace('pred.','').replace('.cif',''):<34}"
              f"{len(shared):>8}{ca_rmsd:>9.2f}{n:>6}{(f'{rms:.2f}' if rms else 'n/a'):>11}")
        out.append((name, len(shared), ca_rmsd, n, rms))
        ALL.append(dict(reference=ref_label, ligand=name, model=f"{folder}/{fn}",
                        n_ca=len(shared), ca_rmsd=round(ca_rmsd, 2),
                        core_atoms=n, core_rmsd=round(rms, 2) if rms else None))

    morph = [o for o in out if o[0] in ("Naltrexone", "Nalbuphine", "Butorphanol")]
    if morph:
        vals = [o[4] for o in morph if o[4]]
        if not vals:
            continue
        print(f"\nmorphinan-core RMSD to naloxone: "
              f"{min(vals):.2f}-{max(vals):.2f} A (mean {np.mean(vals):.2f}) "
              f"across {', '.join(o[0] for o in morph)}")
        print("Methadone is not a morphinan and is reported separately.")


if __name__ == "__main__":
    main()
    import pandas as pd
    out = HERE / "morphinan_core_rmsd.csv"
    pd.DataFrame(ALL).to_csv(out, index=False)
    print(f"\nwrote {out}")

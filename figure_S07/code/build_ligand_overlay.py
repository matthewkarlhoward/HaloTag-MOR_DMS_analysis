#!/usr/bin/env python3
"""
Generate a ChimeraX .cxc script that:
  1. Opens all included MOR experimental structures.
  2. Strips each to receptor chain (+ peptide ligand chain if any), removing
     waters, ions, lipids, sugars, cryoprotectants, and G-protein-bound
     nucleotides.
  3. Matchmaker-aligns every structure to the 8efq reference receptor.
  4. Deletes receptor protein from non-reference models, leaving only the
     ligand poses overlaid in the 8efq orthosteric pocket.
  5. Colors each ligand a unique color and renders cleanly.
"""

from pathlib import Path

CIF_DIR = "/Users/mkh/GitHub/mor_dms_analysis/structures/raw/experimental"
OUT_CXC = "/Users/mkh/GitHub/mor_dms_analysis/structures/scripts/chimerax/ligand_overlay.cxc"

REFERENCE_PDB = "8efq"

STRUCTS = [
    ("8efq", "R", "P"),
    ("4dkl", "A", None),
    ("5c1m", "A", None),
    ("6dde", "R", "D"),
    ("6ddf", "R", "D"),
    ("7sbf", "R", None),
    ("7scg", "D", None),
    ("7t2g", "R", None),
    ("7t2h", "D", None),
    ("7u2k", "D", None),
    ("7u2l", "D", None),
    ("7ul4", "A", None),
    ("8Y72", "R", "D"),
    ("8Y73", "R", "D"),
    ("8e0g", "A", None),
    ("8ef5", "R", None),
    ("8ef6", "R", None),
    ("8efb", "R", None),
    ("8efl", "R", None),
    ("8efo", "R", None),
    ("8f7q", "R", "P"),
    ("8f7r", "R", "P"),
    ("8k9k", "R", "S"),
    ("8k9l", "R", "S"),
    ("9MQI", "A", None),
    ("9ODL", "R", None),
    ("9PXU", "R", None),
    ("9PXW", "R", None),
    ("9PY2", "R", None),
    ("9PY3", "R", None),
    ("9PY4", "R", None),
    ("9WST", "R", "P"),
    ("9WSV", "R", "P"),
    ("9WSW", "R", "P"),
    ("9WSX", "R", "P"),
    ("9bjk", "R", None),
    ("9bqj", "D", None),
]

DROP_RESIDUES = [
    "HOH", "WAT", "NA", "CL", "MG", "K", "CA", "ZN",
    "SO4", "PO4", "GOL", "EDO", "PEG", "P6G", "PG4", "1PE", "MPG",
    "CLR", "OLC", "OLA", "PLM", "VSN",
    "NAG", "FUC", "MAN", "BMA",
    "GTP", "GDP", "GNP",
]

# 37 visually distinct colors (kelly-like palette, hand-picked for contrast on white)
LIGAND_COLORS = [
    "#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231",
    "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE",
    "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#800000",
    "#AAFFC3", "#808000", "#FFD8B1", "#000075", "#808080",
    "#FF1493", "#00CED1", "#FF8C00", "#8A2BE2", "#00FA9A",
    "#DC143C", "#1E90FF", "#FFD700", "#7B68EE", "#20B2AA",
    "#FF6347", "#4682B4", "#D2691E", "#9ACD32", "#BA55D3",
    "#CD5C5C", "#5F9EA0",
]


def main():
    assert STRUCTS[0][0] == REFERENCE_PDB, "Reference must be first in STRUCTS"
    lines = []

    lines.append("# ─── MOR experimental structures: ligand overlay ─────────────")
    lines.append("# Reference receptor: 8efq (human MOR / DAMGO / Gi1)")
    lines.append("")
    lines.append("close session")
    lines.append("set bgColor white")
    lines.append("")

    # Open all structures and strip each
    lines.append("# ─── Open + strip each structure ────────────────────────────")
    for i, (pdb, rec_chain, lig_chain) in enumerate(STRUCTS, start=1):
        chains_to_keep = rec_chain if lig_chain is None else f"{rec_chain},{lig_chain}"
        lines.append(f"# [{i}] {pdb}  receptor=/{rec_chain}" + (f"  peptide_ligand=/{lig_chain}" if lig_chain else ""))
        lines.append(f"open {CIF_DIR}/{pdb}.cif")
        lines.append(f"delete #{i} & ~/{chains_to_keep}")
        drop_sel = ",".join(DROP_RESIDUES)
        lines.append(f"delete #{i}:{drop_sel}")
        lines.append("")

    # Align all to reference (model #1 is 8efq receptor chain R)
    lines.append("# ─── Matchmaker alignment to 8efq receptor ──────────────────")
    nonref_range = f"#2-{len(STRUCTS)}"
    lines.append(f"matchmaker {nonref_range} to #1/{STRUCTS[0][1]} alg sw matrix BLOSUM-62 ssFraction 0.3")
    lines.append("")

    # Strip non-reference receptor protein, keep ligands
    lines.append("# ─── Delete non-reference receptor protein, keep ligands ─────")
    for i, (pdb, rec_chain, lig_chain) in enumerate(STRUCTS, start=1):
        if i == 1:
            continue
        # Delete protein on the receptor chain (this drops the receptor backbone+sidechains)
        # but leaves HETATM ligands on that chain, plus the peptide ligand chain untouched.
        lines.append(f"delete #{i} & protein & /{rec_chain}")
    lines.append("")

    # Render
    lines.append("# ─── Render ─────────────────────────────────────────────────")
    lines.append("# Reference receptor as cartoon")
    lines.append("hide #1 atoms")
    lines.append("show #1 cartoon")
    lines.append("color #1 #DDDDDD")
    lines.append("transparency #1 60 cartoon")
    lines.append("")
    lines.append("# Hide cartoon on non-reference (only ligands remain)")
    lines.append(f"hide #2-{len(STRUCTS)} cartoon")
    lines.append(f"hide #2-{len(STRUCTS)} ribbon")
    lines.append("")
    lines.append("# Show every remaining atom (ligands) as sticks")
    lines.append(f"show #2-{len(STRUCTS)} atoms")
    lines.append(f"style #2-{len(STRUCTS)} stick")
    lines.append("")
    lines.append("# Peptide ligand on reference (DAMGO) also as sticks")
    lines.append(f"show #1/{STRUCTS[0][2]} atoms")
    lines.append(f"style #1/{STRUCTS[0][2]} stick")
    lines.append("")

    # Color each ligand uniquely
    lines.append("# ─── Color each ligand uniquely ─────────────────────────────")
    # Reference DAMGO peptide ligand colored too
    lines.append(f"color #1/{STRUCTS[0][2]} {LIGAND_COLORS[0]}")
    for i, (pdb, rec_chain, lig_chain) in enumerate(STRUCTS, start=1):
        if i == 1:
            continue
        color = LIGAND_COLORS[(i - 1) % len(LIGAND_COLORS)]
        # carbon-only color so heteroatoms keep their element color
        lines.append(f"color #{i} {color}  # {pdb}")
    lines.append("")
    lines.append("# Recolor heteroatoms by element for clarity")
    lines.append("color byhetero")
    lines.append("")

    lines.append("# ─── View ───────────────────────────────────────────────────")
    lines.append("view")
    lines.append(f"view #1/{STRUCTS[0][1]}:50-340  # frame on MOR transmembrane core")
    lines.append("")
    lines.append("# ─── Done ───────────────────────────────────────────────────")
    lines.append("echo Loaded and aligned " + str(len(STRUCTS)) + " MOR structures; non-reference receptors stripped, ligands overlaid.")

    Path(OUT_CXC).write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT_CXC}")
    print(f"  {len(STRUCTS)} structures, reference = {REFERENCE_PDB}")


if __name__ == "__main__":
    main()

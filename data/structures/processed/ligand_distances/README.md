# ligand_distances/

Per-residue distance from each receptor position to the bound orthosteric ligand,
across all 40 ligand-bound MOR PDBs (8qot apo is excluded).

## Producer

- `structures/scripts/ligand_distances/analyze_all_structures.py`
  computes the per-PDB tables.
- `structures/scripts/ligand_distances/harmonize_ligand_distances.py`
  merges the legacy directories into the canonical
  `all_structures_ligand_distances.csv` (receptor-chain-only, 45,864 rows).

## Files

- `all_structures_ligand_distances.csv` — canonical master.
- `all_structures_combined.csv`, `experimental_and_chai_combined.csv` — combined
  legacy/extended views.
- `per_pdb_legacy/` — original per-PDB CSVs kept for provenance.
- `summary_report.txt` — text summary written by the producer.

## Key columns

| column | unit | meaning |
|---|---|---|
| `pdb_id` | — | PDB code |
| `species` | — | mouse / human |
| `drug_name` | — | manually verified orthosteric ligand |
| `receptor_chain` | — | single verified receptor chain |
| `residue_number` | — | PDB residue number |
| `residue_number_human` | — | human-aligned numbering (mouse + 2) |
| `residue_name` | — | three-letter aa |
| `distance_angstrom` | A | min all-atom distance to any ligand atom |
| `distance_sidechain_angstrom` | A | min sidechain-atom distance to ligand (GLY uses all atoms) |
| `shell` | — | `first` (sc <= 5.0 A from ligand), `second` (sc <= 4.5 A from a first-shell sc; not first), or `none` |

# gprotein_distances/

Per-receptor-residue minimum distances to the G-alpha subunit, plus shell
classification, for the MOR-G-protein complexes.

## Producer

- `structures/scripts/gprotein/analyze_gpcr_gprotein_distances.py`
  — base table.
- `structures/scripts/gprotein/analyze_gpcr_gprotein_distances_with_galpha_residue.py`
  — extended version that also records the closest G-alpha residue per
  receptor residue.

## Files

- `<pdb_id>_gprotein_distances.csv` — per-PDB distance + shell table.
- `all_structures_combined.csv` — long-format combined master.
- `summary_report.txt`.

## Key columns

| column | unit | meaning |
|---|---|---|
| `pdb_id`, `residue_number`, `residue_number_human`, `residue_name` | — | receptor residue id |
| `distance_angstrom` | A | min all-atom distance to any G-alpha atom |
| `distance_sidechain_angstrom` | A | min sidechain-atom distance to G-alpha |
| `shell` | — | `first` (sc <= 5.0 A from G-alpha), `second` (sc <= 4.5 A from a first-shell sc), or `none` |
| `galpha_residue_number`, `galpha_residue_name` | — | closest G-alpha residue (extended script only) |
| `galpha_residue_number_sc`, `galpha_residue_name_sc` | — | same, sidechain-based (extended script only) |

8qot (apo) and structures without a G-alpha chain are excluded. Numbering uses
`residue_number_human` (mouse PDBs +2).

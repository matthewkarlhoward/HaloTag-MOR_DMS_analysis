# pi_interactions/

Pi-pi, pi-cation, and pi-ligand interactions for every MOR structure, plus
state comparisons and pair conservation.

## Producers

All under `structures/scripts/interactions/`:

- `extract_pi_interactions.py` — per-structure pi interactions
  (`all_pi_pi.csv`, `all_pi_cation.csv`, `all_pi_ligand.csv`,
  `per_residue_pi_counts.csv`, `extraction_summary.txt`).
- `compare_pi_interactions.py` — per-residue active-vs-inactive deltas
  (`active_vs_inactive_pi_interactions.csv`).
- `pi_pair_conservation.py` — per-pair active / inactive conservation CSVs.

## Files

- `all_pi_pi.csv`, `all_pi_cation.csv`, `all_pi_ligand.csv` — long-format
  per-interaction tables.
- `per_residue_pi_counts.csv` — per-(pdb, residue) counts by interaction type.
- `active_vs_inactive_pi_interactions.csv` — per-residue mean count by state
  and `delta = active - inactive`.
- `extraction_summary.txt`.

## Geometry cutoffs

- Pi-pi stacking: aromatic ring centroid distance <= 5.5 A. Classified as
  `parallel` (dihedral < 30 deg), `t_shaped` (60-90 deg), or `other`.
- Pi-cation: aromatic centroid to ARG/LYS cation atom <= 6.0 A.
- Pi-ligand: aromatic centroid to any ligand heavy atom <= 5.5 A (only in
  ligand-bound structures).

## Key columns

Pair files: `pdb_id`, `state`, `res1_num`, `res1_name`, `res2_num`,
`res2_name`, `interaction_subtype`, `distance_angstrom`. Numbers are
`residue_number_human` (mouse PDBs +2). State from
`structures/metadata/pdb_state_table.csv`.

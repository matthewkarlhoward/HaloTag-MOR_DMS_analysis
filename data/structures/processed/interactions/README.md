# interactions/

Non-bonded intra-receptor interactions (H-bonds, salt bridges, hydrophobic
contacts, water-mediated bridges) for every MOR structure, plus state
comparisons and pair-level conservation tables.

## Producers

All under `structures/scripts/interactions/`:

- `extract_interactions.py` — per-structure H-bonds / salt bridges /
  hydrophobic contacts (`all_hbonds.csv`, `all_salt_bridges.csv`,
  `all_hydrophobic_contacts.csv`, per-residue counts).
- `extract_water_bridges.py` — water-mediated H-bond bridges
  (`all_water_bridges.csv`, `water_bridge_pair_conservation.csv`).
- `compare_interactions.py` — active-vs-inactive per-residue comparison and
  within-state variance (`active_vs_inactive_interactions.csv`,
  `active_within_state_interactions.csv`,
  `inactive_within_state_interactions.csv`, `comparison_summary.txt`).
- `hbond_bb_pair_conservation.py` — sc-bb and bb-bb hbond pair conservation
  CSVs (matching the format of `hbond_sc_sc_pair_conservation.csv`).

## Files

- `all_hbonds.csv`, `all_salt_bridges.csv`, `all_hydrophobic_contacts.csv`,
  `all_water_bridges.csv` — long-format per-pair tables.
- `per_residue_interaction_counts.csv` — per-(pdb, residue) counts.
- `hbond_sc_sc_pair_conservation.csv`, `hbond_sc_bb_pair_conservation.csv`,
  `hbond_bb_bb_pair_conservation.csv`, `salt_bridge_pair_conservation.csv`,
  `hydrophobic_pair_conservation.csv`, `cation_pi_pair_conservation.csv` —
  per-pair active / inactive conservation fractions.
- `active_vs_inactive_interactions.csv`,
  `active_within_state_interactions.csv`,
  `inactive_within_state_interactions.csv`, `comparison_summary.txt`,
  `extraction_summary.txt`.

## Geometry cutoffs

- H-bond: donor-acceptor distance <= 3.5 A, D-H...A angle >= 120 deg.
- Salt bridge: charged atoms (Asp/Glu COO- and Arg/Lys/His NH3+) <= 4.0 A.
- Hydrophobic: C-C distance <= 4.5 A between non-polar sidechains
  (Ala, Val, Leu, Ile, Pro, Phe, Trp, Met).
- Water bridge: receptor polar atom <= 3.5 A from a shared water O atom.

## Key columns

Pair files: `pdb_id`, `state`, `res1_num`, `res1_name`, `res2_num`, `res2_name`
(numbers are `residue_number_human`, mouse PDBs +2).

Conservation files: `res1_num`, `res1_name`, `res2_num`, `res2_name`,
`n_active`, `n_inactive`, `frac_active`, `frac_inactive`.

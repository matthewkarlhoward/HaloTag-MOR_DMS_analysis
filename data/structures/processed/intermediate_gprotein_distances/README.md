# intermediate_gprotein_distances/

Per-receptor-residue contacts to G-protein for intermediate-state MOR
structures (e.g. 9ODL, 9PXW). Used to characterise the
intermediate-vs-active/inactive conformations of the cytoplasmic interface.

## Producer

Intermediate-state contact script in the same family as the gprotein scripts in
`structures/scripts/gprotein/`. The producer for the
intermediate variant is not in the tracked `06_structural/` tree; the schema
mirrors `gprotein_distances/`.

## Files

- `<pdb_id>_contacts.csv` — per-PDB contact / distance table.
- `intermediate_contacts_combined.csv` — combined master across intermediate
  structures.

## Key columns

Same schema as `gprotein_distances/`:
`pdb_id`, `residue_number`, `residue_number_human`, `residue_name`,
`distance_angstrom`, `distance_sidechain_angstrom`, `shell`.

State assignment from `structures/metadata/pdb_state_table.csv`. Numbering
uses `residue_number_human` (mouse PDBs +2).

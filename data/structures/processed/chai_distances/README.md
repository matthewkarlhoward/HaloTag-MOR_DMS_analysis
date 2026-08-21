# chai_distances/

Per-residue ligand distances for the four ChAI-1 predicted MOR-ligand
structures (Naltrexone, Nalbuphine, Methadone, Butorphanol) — drugs with no
experimental MOR structure.

## Producer

`structures/scripts/chai/analyze_chai_predictions.py`

## Files

- `MOR-<drug>_rank<k>_distances.csv` — per-PDB distance + shell table.
- `chai_predictions_combined.csv` — concatenated long-format master.
- `summary_report.txt` — text summary from the producer.

## Key columns

Identical to `ligand_distances/`: `pdb_id`, `drug_name`, `receptor_chain`,
`residue_number`, `residue_number_human`, `residue_name`, `distance_angstrom`,
`distance_sidechain_angstrom` (A), `shell`.

ChAI predictions use human MOR numbering directly (no offset). Receptor
residues 1-60 are excluded (disordered N-terminus in predictions).

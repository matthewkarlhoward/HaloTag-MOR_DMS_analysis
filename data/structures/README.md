# structures/

Structural data and per-residue summaries for 41 mu-opioid receptor (MOR) PDBs
plus 4 ChAI predictions. All numbering is harmonised to **human MOR**: mouse
PDB residue numbers are offset **+2** (e.g. mouse ASP147 = human ASP149). The
column `residue_number_human` is the canonical key across the repo.

## Layout

```
structures/
  raw/                   immutable inputs
    experimental/        mmCIF files for the 41 PDBs (and 9ODL intermediate)
    chai_predicted/      ChAI-1 predictions for ligands without experimental MOR
  processed/             per-residue tables produced by 06_structural scripts
    ligand_distances/
    chai_distances/
    state_rrcs/
    pocket_metrics/
    rotamers/
    backbone_angles/
    ca_displacements/
    ca_displacements_between_state/
    ca_displacements_within_state_variance/
    gprotein_distances/
    arrestin_distances/
    intermediate_gprotein_distances/
    interactions/
    pi_interactions/
    structural_master.csv      one row per residue_number_human (canonical join)
  interactions/          third-party / contributed contact data
    gpcrdb/              GPCRdb non-bonded interaction Excel exports + helpers
    rrcs/                RRCS heatmaps + clustering output
  metadata/
    pdb_state_table.csv  canonical state classification (see below)
```

## Canonical files

- **`metadata/pdb_state_table.csv`** — canonical state classification for every
  structure (3-tier: `active` / `intermediate` / `inactive`). Covers 41 PDBs
  with 4 explicitly excluded. All downstream state-stratified analyses
  (pocket metrics, state RRCS, interactions, ca-displacements) read this file.
- **`processed/ligand_distances/all_structures_ligand_distances.csv`** —
  canonical per-residue ligand-distance master, **receptor-chain-only**,
  45,864 rows across 40 PDBs. Built by
  `structures/scripts/ligand_distances/harmonize_ligand_distances.py`.
- **`processed/structural_master.csv`** — per-position master keyed on
  `residue_number_human`, joining ligand distances, RRCS, interactions,
  rotamers, backbone angles, pocket metrics, ca displacements and gprotein
  distances. Built by
  `structures/scripts/build_structural_master.py`.

## Producer scripts

All Python producers live under
`structures/scripts/`:

| processed/ subdir                              | producer script |
|------------------------------------------------|-----------------|
| `ligand_distances/`                            | `ligand_distances/analyze_all_structures.py` + `ligand_distances/harmonize_ligand_distances.py` |
| `chai_distances/`                              | `chai/analyze_chai_predictions.py` |
| `state_rrcs/`                                  | `rrcs/state_rrcs.py` (consumes `rrcs/calculate_rrcs.py` output) |
| `pocket_metrics/`                              | `pocket/pocket_metrics.py` |
| `rotamers/`                                    | `rotamers/extract_chi_angles.py` (+ `rotamers/compare_chi_angles.py`) |
| `backbone_angles/`                             | `backbone/analyze_gpcr_backbone_angles.py` (+ cluster scripts) |
| `ca_displacements/`                            | `ca_displacements/calculate_ca_distances.py` |
| `ca_displacements_between_state/`              | `ca_displacements/calculate_ca_distances.py` |
| `ca_displacements_within_state_variance/`      | `ca_displacements/calculate_ca_distances.py` (variance pass) |
| `gprotein_distances/`                          | `gprotein/analyze_gpcr_gprotein_distances.py` and `..._with_galpha_residue.py` |
| `arrestin_distances/`                          | arrestin distance script (lives outside `06_structural/`; same logic as gprotein) |
| `intermediate_gprotein_distances/`             | intermediate-state contact script (same family as gprotein) |
| `interactions/`                                | `interactions/extract_interactions.py`, `extract_water_bridges.py`, `compare_interactions.py`, `hbond_bb_pair_conservation.py` |
| `pi_interactions/`                             | `interactions/extract_pi_interactions.py`, `compare_pi_interactions.py`, `pi_pair_conservation.py` |

The two shell helpers in `raw/gpcrdb/` (`create_master.sh`,
`add_columns.sh`) only operate on the GPCRdb Excel exports in that directory.

## Conventions

- **Numbering:** `residue_number` is the PDB number; `residue_number_human` is
  the human-aligned number used as the join key. Mouse PDBs add +2.
- **Receptor chain:** every per-residue table is filtered to a single, manually
  verified receptor chain per PDB (no duplicate ASU copies, no G-protein /
  arrestin / nanobody chains).
- **State:** `active` / `intermediate` / `inactive` come from
  `metadata/pdb_state_table.csv`. State-stratified means/deltas use
  `delta = active - inactive` unless stated otherwise.
- **Distances** are in angstroms; **angles** in degrees.

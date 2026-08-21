# Figure 6 — Efficacy determinants in the orthosteric pocket

**a** Workflow for assigning first- and second-shell contacts, with the
μOR–buprenorphine–miniGi/Gβ/Gγ model shown as the example complex.

First shell = receptor residues with any heavy atom within **4.5 Å** of the ligand;
second shell = residues within **5 Å** of a first-shell residue. Ten of the fifteen
ligands have experimental structures; buprenorphine was determined here
(`data/structures/raw/experimental/buprenorphine.pdb`), carfentanil uses lofentanil
(7T2H) as a close analog, and methadone, naltrexone, nalbuphine and butorphanol use
Chai-1 predictions (`data/structures/raw/chai_predicted/`).

| Script | Original path | Writes |
|---|---|---|
| `code/analyze_all_structures.py` | `structures/scripts/ligand_distances/` | `data/structures/processed/ligand_distances/all_structures_ligand_distances.csv` |
| `code/analyze_chai_predictions.py` | `structures/scripts/chai/` | `data/structures/processed/chai_distances/` |

The canonical merged table used downstream is
`data/structures/processed/ligand_distances/experimental_and_chai_combined.csv`.

**b** Positional average missense score heatmap over every first-shell contact across the
15 ligands; black dots mark actual contacts for that ligand.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/orthosteric_contact_map.R` | `plots/interfaces/orthosteric/` | `data/dms_scores/composite_dms_scores.csv`, `data/structures/processed/ligand_distances/experimental_and_chai_combined.csv` (+ `per_pdb_legacy/buprenorphine_distances.csv`), `data/annotations/GPCRdb_OPRM1_table.csv` | `panels/orthosteric_contact_map_first_shell.pdf` (second-shell output is supplemental figure 7d) |

**c** Radar plots of the average positional loss-of-function score for orthosteric
residues in each ligand class, with TM/ECL arcs outside the position ring.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/radars_arcs.R` | `plots/radars/` | `data/dms_scores/composite_dms_scores.csv` | `panels/circular_lof_{strong,intermediate,weak,antagonist}_filt_arcs.pdf` |

Bars = fraction of variants at that position scoring below the synonymous mean minus one
synonymous SD; rings are 25% increments from 0 to 100. Each panel is a 41.25 mm square —
`coord_polar()` insets the circle, so the square is filled with a negative `plot.margin`
(see the constants at the top of the script).

Drug-name normalisation applied throughout: Oliceridine -> TRV130, Lofentanil ->
Carfentanil, Mitragynine_Pseudoindoxyl -> MP, C6-guano -> C6guano.

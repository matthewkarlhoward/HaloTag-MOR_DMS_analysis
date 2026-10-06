# Figure 8 — Residue networks differ across ligand efficacy

**a** Loss-of-function residue networks for the four efficacy classes, and **b** the
gain-of-function network for antagonists, both rendered on 8EFQ in ChimeraX.

Network definition:
* **node score** = per-position mean disruption count across that class's ligands
  (mutations beyond the synonymous mean ± 2 SD, surface-filtered).
* **nodes** = positions with count >= max(2, Q75 of the non-zero class counts) **and**
  present in the ligand state-contact set (`state_contacts_combined.csv`, 235 positions).
* **edges** = all Cα–Cα pairs < **8.0 Å** among nodes, on PDB **8EFQ**.
* **interactions** = 8EFQ side-chain/side-chain H-bonds, salt bridges and π contacts among
  nodes.
* node radius ∝ √count (area-proportional), range-normalised to 0.55–2.40 Å and scoped
  per direction (LOF anchored to its own max ≈ 17.7, GOF to ≈ 8.0), so equal radii across
  a LOF and a GOF panel do **not** mean equal counts.

| Step | Script | Original path | Reads | Writes |
|---|---|---|---|---|
| 1 | `code/build_merged_class_assets.py` | `network_tools/` | `code/composite_dms_scores_R.csv`, `code/GPCRdb_OPRM1_table.csv`, `data/structures/ca_coords/pdb_8efq_chain_R_ca.csv`, `data/structures/ca_coords/state_contacts_combined.csv`, `data/structures/processed/interactions/`, `data/structures/processed/pi_interactions/` | `chimerax_assets/classes_merged/*.cxc`, `chimerax_assets/classes_merged/*.pb`, `chimerax_assets/attributes/*.defattr` |
| 2 | `code/build_class_variant_assets.py` | `network_tools/` | the assets from step 1 | `chimerax_assets/classes_merged_{nolabel,scaled_labeled}/` display variants |
| 3 | ChimeraX | — | the `.cxc` files | per-class PNG renders |
| 4 | `code/build_class_grid_merged.py` / `build_class_grid_variants.py` | `network_tools/` | the PNGs | stitched class grid |

**Labelling, checked.** Each `class_{lof,gof}_*.cxc` loads the matching
`{lof,gof}_*.defattr` and `.pb`, and the attribute files genuinely hold what their names
say: correlating every deposited `*_count.defattr` against counts recomputed from
`data/dms_scores/composite_dms_scores.csv` gives r = +0.51 to +0.72 for the LOF files
against LOF and +0.83 to +0.93 for the GOF files against GOF, with the cross-terms
negative in all eight cases. An earlier generation of these assets, under
`network_tools/chimerax_scripts/classes/` in the working repository, did carry swapped
labels; the merged set deposited here does not.

**One gotcha.** ChimeraX renders need `windowsize 1600 1200` (matching the save aspect) or the 2D title
   clips and framing shifts run to run. macOS `--nogui` cannot render (no OpenGL) and
   `--offscreen` is Linux-only, so these must be run in the GUI.

**c** ChimeraX zoom on the antagonist GOF cluster (9PXU). No repository code.

**d** TRUPATH Gi1 dose-responses for WT vs A119L across DAMGO, PZM21, nalbuphine and
naloxone. Prism.

**e** Cartoon: wild type to intracellular mutants to A119L to the double mutants, asking
whether the sites are functionally coupled and whether efficacy is recovered. Illustrator.

**f** Heatmap of A119L + mutant activation as % of the wild-type DAMGO Emax, with a
wild-type reference row above a divider.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/heatmap_double_activation.py` | `plots/doubles/` | `data/pharmacology/doubles/doubles_merged_points.csv` | `panels/heatmap_double_activation.pdf` (exactly 40 x 50 mm) |

**d** and **g** Dose-response curves, four ligands per panel, coloured by efficacy class:
wild type and A119L (**d**), K100N and A119L+K100N (**g**). All four panels come from one
script and one file.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/rank_example_curves.py` | `plots/doubles/` | `data/pharmacology/doubles/doubles_merged_points.csv` | `panels/rank_example_curves.pdf` (91 x 97 mm, 2 x 2 grid) |



Both scripts use `code/drc_fit.py` for the F-test-gated sigmoid, and report activation as
the window (-Span, 0 when the F-test fails) as a percentage of the wild-type DAMGO window.
`fig8_bottom_row.py` produced the previous version of this row - a bar chart plus heatmap -
and was dropped when the bar chart left the figure and the heatmap moved to its own script.

Per-variant fitted parameters for every curve are deposited in
`data/pharmacology/doubles/doubles_drc_parameters.csv` and `doubles_metrics_table.csv`,
and in sheet `07_validation_pharmacology` of the supplementary workbook.


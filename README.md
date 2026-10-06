# Distinct activation mechanisms underlie ligand efficacy at a GPCR

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22051140.svg)](https://doi.org/10.5281/zenodo.22051140)

Code and input data for every figure in:

> **Distinct activation mechanisms underlie ligand efficacy at a GPCR**
> Matthew K. Howard, Eve J. Fine, Karthik Srinivasan, Daniel D. Richman, Jerome
> Freudenberg, Zara Weinberg, Jingyou Rao, Ziyue Zou, Christian Macdonald, Patrick
> Rockefeller Grimes, Balazs R. Varga, James S. Fraser, Mark Von Zastrow, Susruta
> Majumdar, Harold Pimentel, Justin English, Ron O. Dror, Aashish Manglik\*, Willow
> Coyote-Maestas\*

A deep mutational pharmacology study of the human μ-opioid receptor: roughly 400,000
measurements of variant effects on signalling and surface expression across 15 opioid
ligands, with full concentration-response scans for morphine, fentanyl and DAMGO.

## Layout

```
data/            all input tables, per-bin variant counts, structures and pharmacology
supplementary_table/  consolidated supplementary data workbook (9 sheets) + CSVs
REPRODUCING.md   what can be regenerated, in what order, and what cannot
MANIFEST.sha256  size + SHA-256 for every file under data/
requirements.txt / r-requirements.txt   Python and R dependencies
CHANGELOG.md     what changed between versions
shared/          upstream pipeline scripts + plotting theme shared by many figures
figure_01 .. 08/ one folder per main figure
figure_S01 .. S11/ one folder per supplemental figure
   code/         the scripts that draw the panels
   panels/       the rendered panel PDFs as they went into the figure
   README.md     panel-by-panel map of script -> input data -> output
```

`data/dms_scores/composite_dms_scores.csv` (10,178 variants x 211 columns) is the master
table that nearly every figure script reads. Its column-by-column provenance is in
`data/dms_scores/COMPOSITE_SCORES_README.md`.

## Figure index

| Figure | Panels | Script | Notes |
|---|---|---|---|
| 1 | a, b, d | — | Illustrator cartoons; flow histogram inset from FlowJo (`data/facs_gating/`). **d** is the library/screening/processing workflow schematic |
| 1 | c | `figure_01/code/compute_scatter_values.py` -> `camp_vs_trupath_plots_black.R` | pooled cAMP vs TRUPATH Gi1 Emax/EC50 |
| 2 | a | — | FlowJo ridgeline export + ChemDraw structures |
| 2 | b | — | Illustrator schematic |
| 2 | c, f | `figure_02/code/dms_dose_response_heatmap.py` | top strip = c, bottom strip = f |
| 2 | d | `figure_02/code/morphine_all_fits_overlay.py` | all ~10k morphine fits + synonymous mean |
| 2 | e | `figure_02/code/position_combo_morphine.py` (`V83`) | example position curves + Emax/EC50 strips |
| 3 | a | `figure_03/code/morphine_variant_ec50_vs_emax.py` | per-variant EC50 x Emax landscape |
| 3 | b, c | `figure_03/code/make_chimerax_param_maps.py` | writes `.defattr` + `.cxc`; rendered in ChimeraX on 8EF6 |
| 3 | d | `figure_03/code/build_overlap_cartoon.py` (structure) + `venn_2SD.py` (Venn) | EC50/Emax overlap; Venn = 29 / 27 / 49 |
| 4 | a | — | cartoon |
| 4 | b | `figure_04/code/variant_morphine_vs_fentanyl_ec50_q126.py` | |
| 4 | c | `figure_04/code/morphine_fentanyl_ec50_reweighting_latest.py` + `chimerax/` | EC50 LOF-bias on 8EF5/8EF6 |
| 4 | d | `figure_04/code/variant_morphine_vs_fentanyl_emax.py` | |
| 4 | e | `figure_04/code/morphine_fentanyl_emax_bias_by_position.py` | least-squares residual track |
| 4 | f | `figure_04/code/build_mor_fent_emax_bias_per_variant.py` | Emax LOF-bias painted in ChimeraX |
| 5 | a | — | cartoon |
| 5 | b | `figure_05/code/dms_stacked_heatmap_magma.R` | 15-ligand heatmap + correlation/ward.D2 dendrogram |
| 5 | c | `figure_05/code/trupath_emax_bar.py` | TRUPATH Gi1 Emax, % DAMGO |
| 5 | d, e, f | `figure_05/code/pca_three_panel_variant_kde.py` | one variant-level PCA fit for all three panels |
| 5 | g | `shared/build_per_position_pc_table.py` | per-position PC2 painted on 8EFQ in ChimeraX; colourbar is now labelled semantically (LOF strong agonists / LOF partial agonists) rather than +/- |
| 6 | a | `figure_06/code/analyze_all_structures.py`, `analyze_chai_predictions.py` | shell definition; buprenorphine cryo-EM model in `data/structures/raw/experimental/buprenorphine.pdb` |
| 6 | b | `figure_06/code/orthosteric_contact_map.R` | first-shell footprint heatmap |
| 6 | c | `figure_06/code/radars_arcs.R` | per-class %LOF radars |
| 7 | a, b | — | ChimeraX views of 8EFQ / 9PXW |
| 7 | c | `figure_07/code/intracellular_gprotein_heatmap_rotated_pctlof_allpos.R` | |
| 7 | d, e, f | — | TRUPATH Gi1 BRET, fitted and drawn in Prism |
| 8 | a, b | `figure_08/code/build_merged_class_assets.py` -> `chimerax_assets/` | LOF/GOF class networks, rendered on 8EFQ |
| 8 | c | — | ChimeraX zoom on 9PXU |
| 8 | d, g | `figure_08/code/rank_example_curves.py` | WT / A119L and K100N / A119L+K100N curve grid (one file, 4 panels) |
| 8 | e | — | cartoon: WT -> intracellular mutants -> A119L -> doubles |
| 8 | f | `figure_08/code/heatmap_double_activation.py` | A119L + mutant activation, WT reference row on top |
| S1 | a | — | ChemDraw |
| S1 | b, c | `figure_S01/code/camp_trupath_composite.py` (or `camp_drc_curves.py` / `trupath_drc_curves.py`) | |
| S1 | d | `figure_S01/code/build_param_table.py` -> `render_param_table.py` | |
| S2 | a | — | DIMPLE schematic |
| S2 | b | `figure_S02/code/per_position_counts.py` | |
| S2 | c | `figure_S02/code/coverage_table.py` | |
| S2 | d, e | — | FlowJo gating (`data/facs_gating/`) |
| S2 | f | `figure_S02/code/04_dose_response_composite.py` | |
| S2 | g | `figure_S02/code/05_multidrug_composite.py` | |
| S2 | h | `figure_S02/code/make_plots.py` | |
| S2 | i | `figure_S02/code/compact_coverage.py` | |
| S3 | a | `figure_S03/code/topdose_vs_surface.py` | |
| S4 | — | `figure_S04/code/per_position_curves_reversed.py` | |
| S5 | a–d | `figure_S05/code/17_concordance_row_figure.py` | Hill vs operational-model fits |
| S6 | a, d | `figure_S06/code/variant_ec50_vs_emax_fen_damgo.py` | |
| S6 | b, c, e, f | `figure_S06/code/make_chimerax_param_maps.py` | |
| S6 | g | `figure_S06/code/distance_vs_delta_ec50_calpha_overlay.py` | |
| S6 | h | `figure_S06/code/venn_altering_positions.py` | |
| S6 | i, j | `figure_S06/code/surface_vs_ec50_emax.py` | |
| S6 | k | `figure_S06/code/surface_residual_structure_maps.py` | |
| S7 | a | `figure_S07/code/emax_vs_effect_scatter.R` | |
| S7 | b, c | `figure_S07/code/build_ligand_overlay.py` + `ligand_overlay.cxc` | |
| S7 | d | `figure_S07/code/orthosteric_contact_map.R` (second-shell output) | |
| S7 | e | `figure_S07/code/pca_biplot_motifs_50mm.py` | |
| S8 | — | — | cryoSPARC processing; no repository code |
| S9 | a–d | `figure_S09/code/build_class_node_edge_plots.py` | |
| S10 | a, b | `figure_S10/code/fig_119_panel.py` | **new** — A119L identification: position waterfall + all 19 substitutions at 119; one script renders the whole figure |
| S11 | a | `figure_S11/code/supp_doubles.py` | |
| S11 | b | — | flow-cytometry surface expression (Prism) |


## What is not here

* **Raw sequencing reads (FASTQ).** Depositing on SRA. The per-bin variant **count** tables derived from them *are*
  included, in `data/variant_counts/` (532 files, 194 MB).

## Citation

If you use this code or data, cite the paper. Author list and metadata are in
`CITATION.cff`.

## License

**Code** (everything under `figure_*/code/`, `shared/`, `supplementary_table/*.py`, and
`data/structures/scripts/`) is MIT — see `LICENSE`.

**Data** (`data/`, `supplementary_table/` outputs, and the rendered panels under
`figure_*/panels/`) is CC BY 4.0 — see `LICENSE-DATA`.

Third-party inputs keep their original terms and are not relicensed: the experimental
structures in `data/structures/raw/experimental/` come from the PDB, and the annotation
tables in `data/annotations/` come from GPCRdb, AlphaMissense, gnomAD and MTR.

Archiving instructions: `ZENODO.md`.

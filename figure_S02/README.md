# Supplemental figure 2 — DMS library generation, screening and QC

**a** DIMPLE workflow schematic. Illustrator.

**b** Positional coverage of the library at ~20x sequencing.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/per_position_counts.py` | `lib_generation/` | `data/library_qc/baseline_1.csv`, `baseline_2.csv` | `panels/per_position_counts.pdf` |

**c** Table of designed vs observed variants.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/coverage_table.py` | `lib_generation/` | `data/library_qc/baseline_{1,2}.csv`, `mor_variants.csv` | `panels/coverage_table.pdf` |

`library_coverage_heatmap.py` produces the position x substitution coverage heatmap used
during library QC.

**d, e** FACS gating for the surface-expression and signalling screens. FlowJo; workspaces,
exported SVGs and three example `.fcs` files are in `data/facs_gating/`.

**f** Score and SE distributions for the three concentration-response screens.
**g** Score and SE distributions for the surface expression and 15-ligand screens.

| Panel | Script | Original path | Reads | Writes |
|---|---|---|---|---|
| f | `code/04_dose_response_composite.py` | `DMS_QC/` | `data/dms_scores/composite_dms_scores.csv` | `panels/04_dose_response_composite.pdf` |
| g | `code/05_multidrug_composite.py` | `DMS_QC/` | same | `panels/05_multidrug_composite.pdf` |

`06_combined_qc_composite.py` stacks the two.

**h** Cells sorted per bin vs gDNA recovered (µg and genome equivalents at 6 pg/cell).
**i** Mean sequencing coverage per condition, stratified by experiment.

| Panel | Script | Original path | Reads | Writes |
|---|---|---|---|---|
| h | `code/make_plots.py` | `figures/variant_counts/plots/` | `data/library_qc/20260522_sort_fastq_gDNA_info.xlsx`, per-bin count tables | `panels/cells_vs_ug_dna.pdf`, `per_bin_coverage_summary.csv` |
| i | `code/compact_coverage.py` | same | `data/library_qc/per_bin_coverage_summary.csv` | `panels/compact_coverage_strip.pdf` |

The per-bin count tables `make_plots.py` walks are included, in `data/variant_counts/`
(532 files), so both panels are reproducible end to end. Its derived output,
`data/library_qc/per_bin_coverage_summary.csv`, is also shipped so panel **i** can be
redrawn without recounting.

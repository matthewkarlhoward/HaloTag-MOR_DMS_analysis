# Changelog

## v1.1.0 — 2026-09-15

Brought in line with the revised manuscript (11 supplemental figures).

**Added**
- `figure_S10/` — Identification of the A119L change-of-efficacy variant. One script
  (`fig_119_panel.py`) renders the whole figure: the position waterfall that nominates
  119, and all 19 substitutions at 119.
- `data/modeling/` — `scores_long.parquet`, `expression.parquet`, `positions.parquet`,
  `ligands.yaml`, the inputs for the above.
- `figure_03` panel d — the EC50/Emax Venn (`venn_2SD.py`; 29 / 27 / 49) and the overlap
  cartoon colouring, plus the ChimeraX assets and the overlap analysis README.
- `figure_08` panels d, f, g — `rank_example_curves.py` (WT / A119L and K100N /
  A119L+K100N) and `heatmap_double_activation.py` (WT reference row above the doubles).
- `data/pharmacology/doubles/doubles_drc_parameters.csv`, `doubles_metrics_table.csv`.
- `requirements.txt`, `r-requirements.txt`, `REPRODUCING.md`, `MANIFEST.sha256`,
  `CHANGELOG.md`.

**Changed**
- `figure_S10` (double-mutant BRET) renumbered to `figure_S11`.
- Figure index updated for Figure 1 panel d and the Figure 5g colourbar, which is now
  labelled semantically rather than as +/-.

**Removed**
- `figure_08/code/fig8_bottom_row.py` and its panel. Superseded: the bar chart is no
  longer in the figure and the heatmap moved to its own script.

**Fixed**
- `07_validation_pharmacology` counted replicates with `rep.nunique()`, but `rep`
  restarts in each of the two double-mutant runs, so pooled WT and A119L were reported
  as n=6 instead of n=10. The per-replicate SEMs were affected the same way.
- `supplementary_table/README.md` repeated the manuscript's ligand-shell definition
  (any heavy atom within 4.5 A). The rule that actually produced the table is side-chain
  atoms within 5.0 A for first shell and 4.5 A for second.

## v1.0.0 — 2026-08-21

Initial export: code and input data for all figures, one folder per main and
supplemental figure, the consolidated supplementary data workbook, and the 532 per-bin
variant count tables.

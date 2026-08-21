# Figure 1 — Deep mutational pharmacology

**a** GPCR activity cartoon and idealised concentration-response curves. Illustrator; no code.

**b** Pooled cAMP inhibition assay schematic. Illustrator. The inset flow-cytometry
histogram (+/- DAMGO eGFP distributions) was exported from FlowJo — workspace and example
`.fcs` in `data/facs_gating/`.

**c** Pooled cAMP vs TRUPATH Gi1 BRET, Emax (% DAMGO) and EC50 (log[ligand]) for the
14-ligand reference panel.

| Step | Script | Original path | Reads | Writes |
|---|---|---|---|---|
| 1 | `code/compute_scatter_values.py` | `trupath_vs_camp_efficacy/` | `data/pharmacology/transcriptional_drc.xlsx`, `data/pharmacology/20260117_mor_wt_trupath_compilation.pzfx` | `camp_trupath_scatter_values.csv` (included in `data/pharmacology/`) |
| 2 | `code/camp_vs_trupath_plots_black.R` | `trupath_vs_camp_efficacy/` | `camp_trupath_scatter_values.csv` | `panels/camp_vs_trupath_plots_black.pdf`, `panels/camp_vs_trupath_legend.pdf` |

`camp_vs_trupath_scatter.R`, `camp_vs_trupath_ec50_scatter.R` and
`camp_vs_trupath_combined.R` are the single-panel and earlier combined variants of the
same plot, kept for reference.

**How the two Emax values are defined.** They are not computed the same way, and this is
deliberate:
* cAMP Emax — three-parameter Hill fit with slope fixed at 1, per ligand, asterisked
  points dropped, reported as `Span` / DAMGO `Span` x 100.
* TRUPATH Emax — peak of the per-dose *mean* activity (1 − normalised BRET), as % of the
  DAMGO peak. Not a fit parameter.

Spearman rho = 0.72 (Emax) and 0.94 (EC50).

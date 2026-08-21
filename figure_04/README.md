# Figure 4 — Ligand-specific determinants of EC50 and Emax

**a** Cartoon. Illustrator.

**b** Variant-level morphine vs fentanyl EC50, Q126 variants highlighted.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/variant_morphine_vs_fentanyl_ec50_q126.py` | `curve_refitting/` | `data/curve_fits/refit_3param_robust_{morphine,fentanyl}.csv` | `panels/variant_morphine_vs_fentanyl_ec50_q126.pdf` |

**c** EC50 LOF-bias mapped on 8EF5 / 8EF6.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/morphine_fentanyl_ec50_reweighting_latest.py` | `dms_scores/` | `data/dms_scores/composite_dms_scores.csv` | per-position EC50 bias table |
| `code/chimerax/recolor_ec50_bias.cxc` + `ec50_bias*.defattr` | `manuscript_revision/chimerax/` | the bias table | rendered side and top-down views |

**d** Variant-level morphine vs fentanyl Emax.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/variant_morphine_vs_fentanyl_emax.py` | `curve_refitting/` | `data/curve_fits/refit_3param_robust_{morphine,fentanyl}.csv` | `panels/variant_morphine_vs_fentanyl_emax.pdf` |

**e** Residuals from the least-squares regression of position-average morphine on fentanyl
Emax; positive (blue) = larger effect on morphine, negative (red) = larger effect on
fentanyl. Black dots mark ligand-contacting positions.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/morphine_fentanyl_emax_bias_by_position.py` | `plots/scatter/ec50_emax/` | `data/dms_scores/composite_dms_scores.csv` | `panels/morphine_fentanyl_emax_bias_by_position.pdf` |

**f** The same Emax bias painted on structure.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/build_mor_fent_emax_bias_per_variant.py` | `dms_scores/` | `data/dms_scores/composite_dms_scores.csv` | `mor_fent_emax_bias_per_variant.csv` -> ChimeraX |

`morphine_fentanyl_emax_reweighting_latest.py` is the per-position aggregation of the
same quantity.

**LOF bias definition.** At each position, measure how strongly mutations reduce one
ligand's potency or efficacy relative to that ligand's own synonymous distribution; the
LOF bias is the difference in that quantity between the two ligands. Positions where
mutations preferentially weaken one ligand carry the largest scores.

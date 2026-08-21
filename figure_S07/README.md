# Supplemental figure 7 — Emax vs saturating-ligand-concentration DMS scores

**a** Position-average Emax from the concentration-response experiment vs the
single-concentration (10 µM) score, for morphine (R² = 0.809), fentanyl (0.764) and
DAMGO (0.679). Validates the saturating-dose screens as an Emax proxy.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/emax_vs_effect_scatter.R` | `plots/scatter/ec50_emax/` | `data/dms_scores/composite_dms_scores.csv` | `panels/emax_vs_effect_scatter.pdf` (165 x 50 mm) |

**b, c** Structural overlays: the Chai-1 predicted poses (butorphanol, nalbuphine,
naltrexone) against experimental structures, and DAMGO/morphine/methadone in the
orthosteric site.

| Script | Original path | Reads |
|---|---|---|
| `code/build_ligand_overlay.py` + `code/ligand_overlay.cxc` | `structures/scripts/chimerax/` | `data/structures/raw/experimental/`, `data/structures/raw/chai_predicted/` |

**d** Positional average missense score heatmap over all **second**-shell contacts across
the assayed ligands — the second output of the figure 6b script.

| Script | Original path | Writes |
|---|---|---|
| `code/orthosteric_contact_map.R` | `plots/interfaces/orthosteric/` | `panels/orthosteric_contact_map_second_shell.pdf` |

**e** PC1/PC2 biplot showing canonical class A activation motifs concentrating at positive
PC1.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/pca_biplot_motifs_50mm.py` | `plots/pca/position/` | `data/pca/pca_position_scores.csv`, `data/annotations/motifs.csv` | `panels/pca_biplot_motifs_50mm.pdf` |

`pca_biplot_motifs_loadings_50mm.py` is the ligand-loading overlay version.

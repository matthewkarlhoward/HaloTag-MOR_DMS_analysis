# Figure 5 — Partial agonists rely on different residues than full agonists

**a** Chemical/efficacy space cartoon. Illustrator + ChemDraw.

**b** Heatmap of variant effects for the untreated and 15 ligand conditions, with the
ligand dendrogram on the right.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/dms_stacked_heatmap_magma.R` | `plots/heatmaps/` | `data/dms_scores/composite_dms_scores.csv` | `panels/dms_heatmap_stacked_magma.pdf` (+ `_sym`, `_lof` scalings) |
| `code/dms_heatmap_stacked_magma_colorbar.py` | `plots/heatmaps/` | same | `panels/dms_heatmap_stacked_magma_colorbar_h.pdf` |

Clustering: **1 − Pearson** distance with **ward.D2** linkage on each ligand's
per-position mean missense effect, positions 65–355. Correlation distance clusters by
profile *shape* rather than effect magnitude; plain Euclidean distance produced a
magnitude artefact that stranded DAMGO with TRV130/carfentanil. The tree is rotated
toward an efficacy ladder and DAMGO pinned to the bottom row. Antagonist, weak partial
and efficacious-agonist groups resolve cleanly; the intermediate-vs-strong agonist split
is not a distinct fingerprint.

**c** TRUPATH Gi1 Emax per ligand, normalised to DAMGO, bars coloured by the cluster
membership from **b**.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/trupath_emax_bar.py` | `trupath_vs_camp_efficacy/` | `data/pharmacology/20260117_mor_wt_trupath_compilation.pzfx` | `panels/trupath_emax_bar.pdf` |

**d, e, f** All three panels come from a single **variant-level** PCA fit (PC1 49.1%,
PC2 10.2%), so the percentages are consistent across the row.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/pca_three_panel_variant_kde.py` | `plots/pca/variant/` | `data/dms_scores/composite_dms_scores.csv`, `data/pca/pca_variant_scores.csv`, `data/pca/pca_position_scores.csv`, `data/pca/pca_loadings_table.csv`, `data/pharmacology/WT_cAMP_vs_TRUPATH_All_Drugs.xlsx` | `panels/pca_three_panel_variant_contour.pdf` (+ `_density` variant) |

* **d** PC1 x PC2 contour of the ~5.5k complete-case missense variants.
* **e** Variant PC1 vs per-variant surface expression score, Spearman rho = −0.51.
* **f** Per-ligand PC2 **loading** vs TRUPATH Gi1 Emax, Spearman rho = −0.96.

`pca_three_panel_variant.py` is the scatter (non-KDE) version;
`pca_pc2_vs_trupath_efficacy.py` is the position-level equivalent of panel f.
Terminology: variants get **scores**, ligands get **loadings**.

Sign convention: variant PCs are flipped, if needed, to match the position-level PCA
(per-position mean variant PC1 vs position PC1, r = +1.00).

**g** Per-position PC2 painted on 8EFQ in ChimeraX (green = positive, purple = negative).

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `shared/build_per_position_pc_table.py` | `plots/pca/` | `data/pca/pca_variant_scores.csv`, `data/pca/pca_position_scores.csv` | `data/pca/per_position_pc_scores.csv` |

## Upstream

`shared/pca_analysis.py` (position level) and `shared/pca_variant_analysis.py` (variant
level) build the score/loading tables in `data/pca/`. Features are the 15 ligand effect
columns, z-scored per ligand, then SVD.

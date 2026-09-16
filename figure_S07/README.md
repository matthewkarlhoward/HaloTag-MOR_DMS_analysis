# Supplemental figure 7 — Emax vs saturating-ligand-concentration DMS scores

**a** Position-average Emax from the concentration-response experiment vs the
single-concentration (10 µM) score, for morphine (R² = 0.809), fentanyl (0.764) and
DAMGO (0.679). Validates the saturating-dose screens as an Emax proxy.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/emax_vs_effect_scatter.R` | `plots/scatter/ec50_emax/` | `data/dms_scores/composite_dms_scores.csv` | `panels/emax_vs_effect_scatter.pdf` (165 x 50 mm) |

**Morphinan-core RMSD to naloxone.** `code/morphinan_core_rmsd.py` superposes each Chai-1
model's receptor onto a naloxone reference (Kabsch on shared Cα, positions 65–355), applies
that transform to the ligand, and takes the RMSD over the maximum common substructure with
naloxone — so agreement is measured in the pocket frame rather than by independently
best-fitting the two ligands. Results in `code/morphinan_core_rmsd.csv`.

Against **9PXY** (naloxone, active, Gi1-bound) the predicted cores agree to
**0.62 / 0.75 / 0.80 Å** for naltrexone / nalbuphine / butorphanol. Against **9PXU**
(inactive, Nb6-bound, the structure used for contact assignment) the same poses give
1.59 / 1.63 / 1.97 Å — the difference is dominated by the active-versus-inactive receptor
(Cα RMSD 3.0 Å vs 1.3 Å), not by the ligand pose. Methadone is not a morphinan and is
excluded from the statement.

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

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

**Two gotchas, both preserved deliberately.**
1. The LOF/GOF labels in the `.cxc` filenames are **swapped** relative to
   `shared/network_analysis_shared.R`: `class_lof_*.cxc` carries the shared script's
   `count_gof_2sd_surf` data and vice versa. Verified against the class-sharing summary.
2. ChimeraX renders need `windowsize 1600 1200` (matching the save aspect) or the 2D title
   clips and framing shifts run to run. macOS `--nogui` cannot render (no OpenGL) and
   `--offscreen` is Linux-only, so these must be run in the GUI.

**c** ChimeraX zoom on the antagonist GOF cluster (9PXU). No repository code.

**d** TRUPATH Gi1 dose-responses for WT vs A119L across DAMGO, PZM21, nalbuphine and
naloxone. Prism.

**e, f** Double-mutant summary: %-activation bars for the K100N series and the A119L
rescue heatmap across all partners.

| Step | Script | Original path | Reads | Writes |
|---|---|---|---|---|
| 1 | `code/extract_prism_points.py` | `plots/doubles/` | the two source `.prism` projects (not distributed) | `data/pharmacology/doubles/doubles{1,2}_points.csv` |
| 2 | `code/drc_fit.py` | `plots/doubles/` | — | `fit_drc()`: unconstrained three-parameter sigmoid plus an extra-sum-of-squares F-test against a flat line. Fitted Span kept (any sign) if p < 0.05, otherwise "no response", Span = 0, drawn flat. |
| 3 | `code/fig8_bottom_row.py` | `plots/doubles/` | `data/pharmacology/doubles/doubles_merged_points.csv` | `panels/fig8_bottom_row.pdf` (exact 165 x 40 mm) |

The heatmap carries a **wild-type reference row on top**, separated from the double-mutant
rows by a rule, so each A119L + partner value can be read against WT in the same units.

Signal-down assay: baseline Bottom ≈ 1, activation drives the BRET ratio down, so the
activation window is `1 − Top = −Span`. Percentages are window / WT-DAMGO window x 100.
The two double-mutant runs are co-scaled (shared-anchor r = 0.98) and pooled for WT and
A119L (n = 4 + 6 = 10); they are **not** re-normalised to DAMGO, because WT is the
noisiest anchor.

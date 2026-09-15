# data/

All inputs consumed by the figure scripts. Original repository paths are given so that
each file can be traced back to the analysis tree it came from.

## variant_counts/ — primary read counts

**532 per-bin count files** across the five screens (surface, DAMGO/morphine/fentanyl
concentration-response, and the 15-ligand panel), plus sequencing sample sheets. These
are the Lilace inputs that every score in `dms_scores/` is derived from. One file per
sorted bin; `hgvs` is the join key. Full layout, file format and naming conventions in
`variant_counts/README.md`.

## dms_scores/ — variant scores

| File | Original path | Contents |
|---|---|---|
| `composite_dms_scores.csv` | `dms_scores/composite_dms_scores.csv` | **Master table**, 10,178 variants x 211 columns. Per-ligand saturating-concentration effects, per-dose rescaled scores for the three DRC ligands, surface expression, Hill parameters, and merged annotation. Read by almost every figure script. |
| `COMPOSITE_SCORES_README.md` | same | Column-by-column provenance of the master table. |
| `drc/Morphine_wide_hill_ordinal_df.csv` | `dms_scores/drc/` | Lilace scores for the 8-point morphine dose series (wide). |
| `drc/Fentanyl_wide_hill_ordinal_df.csv` | `dms_scores/drc/` | Lilace scores for the 8-point fentanyl dose series (wide). |
| `drc/DAMGO_drc.csv` | `dms_scores/drc/` | Lilace scores for the 6-point DAMGO dose series. |
| `surface/surface_expression_scores.tsv` | `dms_scores/surface/` | HaloTag/JF635i surface expression screen scores. |
| `multidrug/all_drug_scores_mapped_normalized.tsv.gz` | `dms_scores/multidrug/` | 15-ligand saturating-concentration Lilace scores (gzipped; 83 MB uncompressed). |
| `lof_gof_scores.csv`, `lof_gof_group_scores.csv` | `plots/lof_gof/` | Per-position fraction of LOF/GOF missense variants (>2 SD from the synonymous mean), per ligand and per efficacy class. Produced by `shared/compute_lof_gof.R`. |

## curve_fits/ — three-parameter Hill fits

`refit_3param_robust_{morphine,fentanyl,damgo}.csv`
(original: `curve_refitting/`). Per-variant `fitted_emin`, `fitted_ec50_logM`,
`fitted_emax`, fit diagnostics and `new_curve_type`. Produced by
`shared/refit_3param_robust.py` and `shared/refit_3param_robust_damgo.py` from
`composite_dms_scores.csv`. Emax is on an activity scale with synonymous
Emin = 0 / Emax = 1; EC50 is log10 M.

## operational_model/ — Black–Leff fits

`{morphine,fentanyl,damgo}_operational_fits.csv` (original: `operational_model/`).
Per-variant tau (operational efficacy), rho (expression-corrected per-receptor
coupling), KA, baseline and fit diagnostics. Produced by
`shared/01_fit_operational.py`; see `shared/README_operational_model.md` for the model
and its anchors.

## pca/

| File | Original path | Contents |
|---|---|---|
| `pca_position_scores.csv` | `plots/pca/position/` | Position-level PCA scores (PC1 66.8%, PC2 15.6%). |
| `pca_loadings_table.csv` | `plots/pca/position/` | Per-ligand PC loadings + efficacy-class assignment. |
| `pca_variant_scores.csv` | `plots/pca/variant/` | Variant-level PCA scores (PC1 49.1%, PC2 10.2%) — the fit used in figure 5d–f. |
| `per_position_pc_scores.csv` | `plots/pca/` | Position-averaged variant PCs; the table painted on the structure in figure 5g. |

Terminology: positions/variants get **scores**, ligands get **loadings**.

## modeling/ — inputs for the A119L identification figure

`scores_long.parquet` is the long-format per-variant x per-ligand score table used by the
modeling analyses, and `ligands.yaml` carries each ligand's TRUPATH Emax and holdout flag.
Together they drive supplemental figure 10 (`figure_S10/code/fig_119_panel.py`): the
waterfall that nominates position 119, and the per-substitution efficacy slopes.
Original paths: `mor_efficacy/data/processed/` and `mor_efficacy/configs/`.

## annotations/

GPCRdb OPRM1 numbering + secondary structure (`GPCRdb_OPRM1_table.csv`), class A motif
membership (`motifs.csv`), human↔mouse position map (`oprm_human_mouse_pos.csv`),
AlphaMissense, MTR and gnomAD tables. Original: `annotations/`.

## pharmacology/ — low-throughput assays on wild-type and single variants

| File | Contents |
|---|---|
| `transcriptional_drc.xlsx` | Raw pooled-cAMP transcriptional reporter dose-responses, wild-type μOR, all ligands. |
| `20260117_mor_wt_trupath_compilation.pzfx` | Prism project with the wild-type TRUPATH Gi1 BRET dose-responses. |
| `WT_cAMP_vs_TRUPATH_All_Drugs.xlsx` | Per-ligand Emax/EC50 summary for both assays. **Note:** in the `cAMP vs TRUPATH (WT)` sheet the two Emax headers are swapped — the column labelled `cAMP Response` holds the TRUPATH Gi1 Emax. |
| `camp_trupath_scatter_values.csv` | Derived table behind figure 1c (`compute_scatter_values.py`). |
| `param_table.csv` | Derived table behind supplemental figure 1d. |
| `doubles/doubles{1,2}.xlsx` | Prism fit tables for the two double-mutant BRET runs. |
| `doubles/doubles{1,2}_points.csv` | Replicate-level points extracted from the source `.prism` projects by `extract_prism_points.py`. |
| `doubles/doubles_merged_points.csv` | The two runs merged (WT and A119L pooled, n = 10); input to figures 8e–f and S11a. |

Signal-down convention: baseline ≈ 1, activation drives BRET ratio down, so the
activation window is `1 − Top = −Span`.

## library_qc/

`baseline_1.csv`, `baseline_2.csv`, `mor_variants.csv`, `dimple_paper_oligo.csv`
(original: `lib_generation/`) — designed vs observed library composition and per-position
counts. `per_bin_coverage_summary.csv` (original: `figures/variant_counts/plots/`) —
mean/median coverage per sorted bin, precomputed from the count tables in
`variant_counts/`; its `variant_counts_filename` column is the join key back to them, and
it doubles as the sample manifest. `20260522_sort_fastq_gDNA_info.xlsx` — cells sorted and
gDNA recovered per sample. `20260522_MOR_oligos_primers.xlsx` — oligo and primer tables.

## structures/

Copied from `structures/` in the analysis repo; see `structures/README.md` for the full
convention (all numbering harmonised to **human** μOR, mouse PDBs offset +2, key column
`residue_number_human`).

* `raw/experimental/` — 42 mmCIF files plus `buprenorphine.pdb`, the μOR–buprenorphine–miniGi
  model determined in this study.
* `raw/chai_predicted/` — Chai-1 predictions for butorphanol, methadone, nalbuphine and naltrexone.
* `processed/ligand_distances/` — per-residue ligand distances; `experimental_and_chai_combined.csv`
  is the canonical first/second-shell source for figure 6.
* `processed/gprotein_distances/`, `processed/intermediate_gprotein_distances/` — receptor–Gα
  contacts in active and latent/intermediate states (figure 7c).
* `processed/interactions/`, `processed/pi_interactions/` — H-bond, salt-bridge and π contacts
  (network edge annotation, figure 8).
* `metadata/pdb_state_table.csv` — canonical active/intermediate/inactive classification.
* `ca_coords/pdb_8ef6_chain_R_ca.csv`, `pdb_8efq_chain_R_ca.csv` — Cα coordinates used for the
  8 Å contact graphs; `state_contacts_combined.csv` — per-position ligand-shell flags used to
  restrict network nodes.
* `scripts/` — the producers for everything under `processed/`.

## facs_gating/

FlowJo workspaces (`.acs`), exported gating SVGs and three example `.fcs` files behind
figure 1b and supplemental figures 2d–e.

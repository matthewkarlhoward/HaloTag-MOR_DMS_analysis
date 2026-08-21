# shared/

Upstream pipeline scripts and the plotting theme used by more than one figure. Originals
live at the repository paths given below.

| File | Original path | Role |
|---|---|---|
| `create_composite_scores.py` | repo root | Builds `data/dms_scores/composite_dms_scores.csv` from the Lilace outputs (three dose-response screens, the 15-ligand saturating panel, surface expression) plus the annotation tables. Filter: `mapped == 'raw'` and `lilace_run == drug`; sign flip retained; per-ligand columns named `{Drug}_effect`. Column provenance in `data/dms_scores/COMPOSITE_SCORES_README.md`. |
| `refit_3param_robust.py` | `curve_refitting/` | Three-parameter Hill fits for morphine and fentanyl -> `data/curve_fits/`. Slope fixed at 1, bounded least squares, one robust refit pass (points beyond 2x RMS residual dropped, at most two per variant). Emin/Emax bounded to [−2, 2], logEC50 to the assayed range. Only fits reaching both asymptotes within the tested range are kept. Emax is placed on an activity scale with synonymous Emin = 0 and Emax = 1; EC50 is log10 M. |
| `refit_3param_robust_damgo.py` | `curve_refitting/` | Same for the six-point DAMGO series (ligand-free condition included as a baseline anchor at −12). |
| `01_fit_operational.py` | `operational_model/` | Black–Leff operational fits -> `data/operational_model/`. See `README_operational_model.md`. |
| `README_operational_model.md` | `operational_model/README.md` | Model definition, fixed anchors, and the expression decomposition. |
| `pca_analysis.py` | `plots/pca/position/` | Position-level PCA -> `data/pca/pca_position_scores.csv`, `pca_loadings_table.csv`. PC1 66.8%, PC2 15.6%. |
| `pca_variant_analysis.py` | `plots/pca/variant/` | Variant-level PCA -> `data/pca/pca_variant_scores.csv`. PC1 49.1%, PC2 10.2%. Sign-aligned to the position PCA. |
| `build_per_position_pc_table.py` | `plots/pca/` | Position-averaged variant PCs -> `data/pca/per_position_pc_scores.csv` (figure 5g). |
| `compute_lof_gof.R` | `plots/lof_gof/` | Per-position fraction of missense variants beyond ±1 SD / ±2 SD of the synonymous mean, per ligand and per efficacy class -> `data/dms_scores/lof_gof_scores.csv`, `lof_gof_group_scores.csv`. |
| `network_analysis_shared.R` | `network_tools/` | Shared data prep for the network figures: surface filtering, disruption counts, GPCRdb SSE join, class colours and ordering. Sourced by the network plot scripts. |
| `build_structural_master.py` | `structures/scripts/` | Joins every per-residue structural table into `data/structures/processed/structural_master.csv`. |
| `howard_theme.R` | repo root | ggplot2 house theme. |
| `fig_template.R` | repo root | Minimal starting template for a new R panel. |

## Figure style

Helvetica 6 pt, black text, 0.5 pt rules throughout; panel letters bold 8 pt. PDFs are
written with `pdf.fonttype = 42` (TrueType, not Type 3) and `transparent = True` so they
stay editable in Illustrator.

Two recurring mechanical gotchas, both worked around in the scripts:
* matplotlib on the `macosx` backend snaps `figsize` to two decimal places; panels that
  must be an exact millimetre size force `mpl.use('Agg')`, set the size explicitly and
  never use `bbox_inches='tight'`.
* R's `pdf()` device floors the MediaBox to integer points, so a 165 mm request comes out
  at 164.75 mm. Use SVG where the exact width matters.

## Ligand naming

Normalised everywhere to: Oliceridine -> **TRV130**, Lofentanil -> **Carfentanil**,
Mitragynine Pseudoindoxyl -> **MP**, C6-guano -> **C6guano**.

Efficacy classes (from the figure 5b clustering, confirmed against TRUPATH Emax):
* **antagonist** — naloxone, naltrexone
* **weak** — buprenorphine, nalbuphine, butorphanol, MP
* **intermediate** — PZM21, TRV130, methadone
* **strong** — C6-guano, morphine, fentanyl, carfentanil, SR-17018, DAMGO

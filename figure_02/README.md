# Figure 2 — Concentration-response deep mutational scan of the μOR

**a** Ridgeline of library eGFP distributions across the morphine dose series; exported
from FlowJo (`data/facs_gating/`). Ligand structures drawn in ChemDraw.

**b** Surface expression workflow schematic. Illustrator.

**c** and **f** Stacked variant heatmaps. Both strips come from one script.

| Panel | Script | Original path | Reads | Writes |
|---|---|---|---|---|
| c | `code/dms_dose_response_heatmap.py` | `plots/heatmaps/` | `data/dms_scores/composite_dms_scores.csv` | `dms_dose_response_heatmap_top.pdf` (morphine / fentanyl / DAMGO dose series + No Ligand + surface expression) |
| f | same script | | same | `dms_dose_response_heatmap_bottom.pdf` (per-ligand Emin / EC50 / Emax) |
| c, f | `code/dms_heatmap_stacked_magma_colorbar.py` | `plots/heatmaps/` | same | standalone colourbars, vector-safe (256 rects) |

Colour: matplotlib `magma` resampled at `t**0.55` ("magma_fastyellow"). Scores are
signal-down, so dark = high signalling and bright yellow = loss of signalling. White tiles
are missing data or parameters that could not be determined.

**d** All ~10,000 morphine concentration-response fits overlaid in black with the
synonymous mean in red.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/morphine_all_fits_overlay.py` | `curve_refitting/` | `data/curve_fits/refit_3param_robust_morphine.csv` | `panels/morphine_all_fits_overlay.pdf` |

**e** Example position (V83, 1x45): per-variant curves plus raw points, with the derived
EC50 and Emax distributions to the right.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/position_combo_morphine.py` | `curve_refitting/` | `data/dms_scores/composite_dms_scores.csv`, `data/curve_fits/refit_3param_robust_morphine.csv` | `panels/position_combo_morphine_V83.pdf` |

`position_beeswarm_morphine.py` and `position_showcase_morphine.py` are the two halves of
that panel drawn separately.

## Upstream

The Hill parameters plotted in **e** and **f** come from `shared/refit_3param_robust.py`
(morphine, fentanyl) and `shared/refit_3param_robust_damgo.py`. Three-parameter Hill,
slope fixed at 1, bounded least squares; each curve is fit once, points deviating by more
than twice the RMS residual are removed (at most two per variant) and the curve refit.
Only fits reaching both asymptotes inside the tested range are retained.

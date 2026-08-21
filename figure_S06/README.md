# Supplemental figure 6 — Additional EC50 and Emax data across morphine, fentanyl and DAMGO

| Panel | Content | Script | Original path | Reads |
|---|---|---|---|---|
| a, d | Per-variant EC50 x Emax landscapes for fentanyl and DAMGO | `code/variant_ec50_vs_emax_fen_damgo.py` | `curve_refitting/` | `data/curve_fits/refit_3param_robust_{fentanyl,damgo}.csv` |
| b, c | Fentanyl EC50 / Emax per-position averages on 8EF5 | `code/make_chimerax_param_maps.py` | `curve_refitting/` | `data/curve_fits/refit_3param_robust_fentanyl.csv` |
| e, f | DAMGO EC50 / Emax per-position averages on 8EFQ | same | same | `..._damgo.csv` |
| g | Cα distance to ligand vs ΔEC50, per ligand, with exponential fits | `code/distance_vs_delta_ec50_calpha_overlay.py` | `plots/scatter/distance_shell/` | `data/dms_scores/composite_dms_scores.csv`, `data/structures/processed/ligand_distances/` |
| h | Venn diagrams of positions altering EC50 and Emax across the three ligands | `code/venn_altering_positions.py` | `curve_refitting/` | `data/curve_fits/refit_3param_robust_*.csv` |
| i, j | Variant-level surface expression vs morphine EC50 (rho = −0.30, n = 6,254) and Emax (rho = 0.63, n = 7,104), loess in red | `code/surface_vs_ec50_emax.py` | `plots/scatter/surface/` | `data/dms_scores/composite_dms_scores.csv`, `data/curve_fits/refit_3param_robust_morphine.csv` |
| k | Expression-adjusted morphine EC50 and Emax on 8EF6 | `code/surface_residual_structure_maps.py` | `plots/scatter/surface/` | same |

Panels b, c, e, f share the ChimeraX generator with figure 3 — running
`make_chimerax_param_maps.py` once writes the `.defattr`/`.cxc` for every ligand and
parameter. Panel k writes its own attribute files into `code/chimerax/`;
`render_residual_maps.cxc` renders both views.

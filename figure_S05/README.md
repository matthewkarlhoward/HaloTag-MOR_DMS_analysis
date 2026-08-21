# Supplemental figure 5 — Hill vs operational fits

**a–c** Variant-level concordance between the three-parameter Hill fit and the Black–Leff
operational fit for potency, observed efficacy and baseline (n = 7,146).
**d** Operational Emax vs log10 rho (per-receptor coupling), coloured by surface expression:
low-Emax/rho≈0 variants are expression-limited, low-Emax/rho<0 variants have lost coupling.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/17_concordance_row_figure.py` | `operational_model/` | `data/operational_model/morphine_operational_fits.csv`, `data/curve_fits/refit_3param_robust_morphine.csv`, surface-expression bin counts | `panels/fig_operational_concordance_row.pdf` |

`06_hill_concordance.py` and `11_rho_vs_emax.py` are the standalone versions of panels
a–c and d; `14_supp_composite.py` is an earlier multi-row layout.

## The model

`resp(X) = Em · τ·A / (K_A + (1+τ)·A)`, `score(X) = B − resp(X)` (signal-down), with the
transducer slope fixed at 1 and `Em` fixed as a global anchor (q0.995 of panel activation
= 0.695). Efficacy is identifiable without affinity because
`Emax_obs = Em·τ/(1+τ)`. The expression split is a straight subtraction on the log scale,
`log10 ρ = Δlog10 τ − expression_effect`, with `τ_WT` = median τ over well-fit synonymous
variants = 1.53. Affinity is **not** a per-variant deliverable — it is not identifiable
from dose-response data alone. Full description in `shared/README_operational_model.md`.

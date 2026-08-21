# Operational-model efficacy decomposition — morphine DMS (prototype)

Fits the Black–Leff **operational model of agonism** to every morphine variant's
8-point dose-response curve, then uses the **surface-expression** score to split the
efficacy change into an *expression* part and a *per-receptor coupling* part.

We have **no affinity/binding data**, so affinity is deliberately **not** a per-variant
output. The identifiable deliverables are efficacy (τ) and expression-corrected
efficacy (ρ).

## Model (n = 1, one ligand = morphine)

```
A = 10^X                                   X = log10 [morphine] (M)
resp(X) = Em · τ·A / (KA + (1+τ)·A)        operational activation magnitude
score(X) = B − resp(X)                     signal-DOWN (Gi/cAMP) orientation
```

| symbol | meaning | how treated |
|---|---|---|
| `Em`  | system-max activation | **FIXED** anchor = q0.995 of panel activation = **0.695** |
| `n`   | transducer slope | fixed = 1 |
| `B`   | per-variant baseline | free (nuisance) |
| `τ`   | operational efficacy (receptor reserve) | free per variant — **deliverable** |
| `KA`  | agonist–receptor Kd | free but **not identifiable when saturated**; flag channel only |

Closed forms (why efficacy is identifiable without affinity):
```
Emax_obs = Em·τ/(1+τ)      ⟹  τ = Emax_obs/(Em − Emax_obs)     efficacy from the Emax drop
EC50     = KA/(1+τ)                                            EC50 shift; KA is a by-product
```

## Expression decomposition

`τ_v = τ_WT · S_v · ρ_v`, with `S_v` = surface expression relative to WT (linear in
receptor number). Because the expression `effect` is on a log-fold scale (WT ≈ 0):

```
log10 ρ_v = log10 τ_v − log10 τ_WT − log10 S_v
          = Δlog10 τ_v − expr_effect_v         (a straight subtraction)
```

- `τ_WT` = median τ over well-fit **synonymous** variants = **1.53** (log10 = 0.185).
- `ρ_v` = per-receptor coupling efficiency (=1 for WT). **This is the key output:** it
  separates "signals poorly because it doesn't reach the surface" (S↓, ρ≈1) from
  "reaches the surface but couples worse" (S≈1, ρ<1).

## Files

| file | contents |
|---|---|
| `01_fit_operational.py` | fits all variants, merges expression, writes results |
| `02_plots.py` | all summary figures + position table |
| `morphine_operational_fits.csv` | per-variant: B, τ, logKA, logEC50, Emax_obs, SEs, log10_tau, dlog10_tau, log10_S, **log10_rho**, flags (`saturating`, `ka_unidentified`, `clean`) |
| `morphine_position_summary.csv` | position-averaged Δlog10 τ, log10 ρ, expression |
| `plots/fig1_heatmaps` | DMS maps: log10 τ, Δlog10 τ, expression, **log10 ρ** |
| `plots/fig2_decomposition` | efficacy-vs-expression scatter (colored by ρ) + mechanistic classification |
| `plots/fig3_Em_sensitivity` | how WT reserve & ρ dynamic range depend on the Em anchor |
| `plots/fig4_example_fits` | representative operational fits (WT, coupling-loss, expression-loss, super-responder) |
| `plots/fig5_diagnostics` | KA-identifiability collapse on saturation; operational τ vs raw Emax |
| `03_compare.py` / `plots/fig6_compare_params` | operational τ/ρ vs existing Hill EC50/Emax |

## How the new params relate to the old EC50/Emax (fig6)

Spearman, clean missense (n≈6.6k):

| comparison | ρ_s | reading |
|---|---:|---|
| operational EC50 vs Hill EC50 | **+0.86** | potency is the same info, just reparameterized — no new potency axis |
| log10 τ vs Hill Emax | **+0.84** | τ tracks Emax but **nonlinearly** — same ΔEmax means more efficacy change near the WT ceiling (reserve compression) |
| log10 ρ vs Hill Emax | +0.42 | partly related... |
| log10 ρ vs Hill EC50 | **+0.13** | ...but **largely orthogonal** to the existing params — ρ is the genuinely new axis |
| log10 ρ vs expression | −0.48 | by construction (ρ removes expression) |

**Reclassification:** of 533 variants that raw Emax would call LoF, **60% are actually
expression-limited (coupling intact, ρ≈0)** and only **23% are genuine coupling loss**.
That split is invisible in EC50/Emax and is the main thing the operational+expression
model buys you — subject to the survivorship caveat below.

Quality gate for deliverable maps (`clean` column): sigmoid curve_type, RMSE < 0.12,
τ not railed, not saturating → 6589 / 7565 missense.

## Caveats (read before interpreting)

1. **Em is the linchpin.** Absolute τ (and WT reserve) scale strongly with the Em
   anchor (fig3). We chose q0.995 of panel activation, which sits in the *stable*
   regime; DAMGO (a fuller agonist) would set Em properly and is the right long-term fix.
2. **No affinity.** Any EC50 shift not explained by τ shows up as an unstable `logKA`
   (see `ka_unidentified` / fig5-left). Treat it as a *flag* for candidate
   orthosteric/affinity effects, never as a measured Kd.
3. **Survivorship bias in ρ.** At buried positions most substitutions kill expression
   and go flat → filtered out. The surviving fits are biased toward variants that
   retained function, inflating the positive-ρ ("coupling-gain") signal in the TM core
   of fig1-D. Interpret the core ρ enhancement cautiously; a censored-data / hierarchical
   fit would handle the dead variants properly instead of dropping them.
4. **Expression scale base** assumed log10. If the FACS pipeline used log2/ln, rescale
   `EXPR_LOG_BASE` in `01_fit_operational.py`; this linearly rescales the ρ magnitude
   (not the qualitative pattern).

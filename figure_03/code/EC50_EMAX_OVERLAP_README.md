# Are EC50- and Emax-altering positions separable?

Short answer: **no, not as sets of positions** — and the figure-3 panels should not
be described that way. What *is* separable is a potency-specific *component* of
the effect, which is localised to the orthosteric pocket.

All numbers are morphine, per-position means of missense variants (n >= 5),
using the same curve-type filters as figure 3B/C
(`curve_refitting/make_chimerax_param_maps.py`): EC50 from `sigmoid`/`no_baseline`
(a flat curve has no identifiable EC50), Emax from `sigmoid`/`flat`/`no_baseline`
(a dead variant is a real efficacy effect). N = 394 positions.

## What the data say

**1. The two sets overlap far more than chance.** At the paper's own LoF
definition (effect > 1 SD of the synonymous distribution): 150 EC50, 143 Emax,
**96 in both** vs 54 expected by chance (P = 3e-19). 64% of EC50-altering
positions also alter Emax. The overlap is at or above chance at *every*
stringency tested (2-6x from top-20 through top-150); it is never below chance.

**2. Position effects lie on one shared axis.** r = -0.54 between EC50 shift and
Emax (-0.61 corrected for split-half reliability: 0.82 EC50, 0.96 Emax). A
noise-corrected (Deming) slope of **-0.98 log EC50 per unit of Emax** is what the
operational model predicts: dropping coupling efficiency lowers Emax and
right-shifts EC50 through the same parameter. The shared axis accounts for 29% of
the position-level EC50 variance.

**3. The shared class is the conserved activation machinery.** All of the
class A motifs land in "both": Na+ site (D116, D149), DRY (D166, Y168), PIF
(F291), CWxP (V293/C294/W295/T296), NPxxY (N334/P335/L337/Y338), plus the pocket
positions H299, Y301, Y328, W320.

**4. "Selective" positions are the weak, unstable ones.** Median effect in the
single-parameter classes is barely over threshold (46% of EC50-only and 43% of
Emax-only positions are within 1.5 SD, vs 8% of "both"). Bootstrapping the
variants at each position, EC50-only positions keep their class in a median 76%
of resamples and Emax-only in 67%, vs 89% for "both"; 27% of selective positions
land in "both" in more than a quarter of resamples.

**5. At the 1 SD threshold there is no spatial separation either.** Median
sidechain distance to morphine: EC50-only 17.0 A, Emax-only 16.5 A (P = 0.76).
Spatial separation appears only in the extreme tail: at top-15/20 the EC50 set
sits at 8.4-8.9 A vs 19.9 A for Emax (P = 0.002-0.003, ~50% vs ~5% within 8 A),
and it is gone by top-50 (P = 0.12).

**6. What IS separable: the potency-specific residual.** Regressing the EC50
shift on Emax (Deming, error ratio from the reliabilities) and taking the
residual isolates the EC50 change that the efficacy loss does not explain. That
residual is pocket-local: median +0.31 log within 8 A vs +0.02 beyond
(P = 1.3e-5; rho = -0.17 with distance, P = 0.005). Top potency-specific
positions are H299, V302, F239, Y328, W320, W295, L221, Y301 (all 4-8 A from
morphine) plus the G-protein face R279, Y254, P311. Top efficacy-specific are
R167 (DRY), A170, S164, N88, M257, L261 — distal, in the core.

## Suggested wording

Not: "positions altering EC50 and Emax are separable / occupy distinct sets".

Instead: *"Most positions alter potency and efficacy together, along a single
coupling axis of slope ~1 log unit per unit Emax, as expected if both read out
the efficiency of receptor activation; the conserved activation motifs behave
this way. Superimposed on that shared axis is a potency-specific component —
EC50 shifts larger than the accompanying Emax loss predicts — that is
concentrated at the orthosteric pocket and, secondarily, at the G-protein
coupling site."*

That claim is supported by 6 above, survives the reviewer's overlap question, and
is more informative than a dissociation claim the data do not support.

## Files

| file | what |
|---|---|
| `overlap_core.py` | position table, classes, Deming residual, bootstrap |
| `fig_overlap_2d.py` -> `fig_overlap_2d.pdf/.png` | 8-panel quantitative figure |
| `build_chimerax_overlap.py` | writes the ChimeraX assets below |
| `chimerax/overlap_classes_{side,top}.cxc` | 3-class spheres: blue EC50-only (54), red Emax-only (47), purple both (96) |
| `chimerax/overlap_shared_only_{side,top}.cxc` | the 96 shared positions alone |
| `chimerax/potency_selective_{side,top}.cxc` | continuous map of the potency-specific residual |
| `chimerax/render_all.cxc` | renders all six views to PNG |
| `snakeplot_overlap.py` -> `snakeplot_overlap.svg/.png` | 2D snakeplot recoloured by class; no ChimeraX needed |
| `build_panel3d_variants.py` | v1 published / v2 refiltered / v3 overlap versions of panel 3D |
| `fig_panel3d_compare.py` | assembles the three rendered variants side by side |
| `ec50_emax_overlap_positions.csv` | per-position values, class, residual, distance |
| `stringency_sweep.csv` | overlap and distance vs set size |

Cameras are the captured figure-3 views, so these panels drop into the existing
figure without re-orienting. Render with ChimeraX open (headless rendering is not
available on macOS):

    open plots/ec50_emax_overlap/chimerax/render_all.cxc       # the three overlap maps
    open plots/ec50_emax_overlap/chimerax/render_panel3d.cxc   # the three panel-3D variants

Both scripts begin with `close session`, so run them in a ChimeraX window you do
not mind clearing. Then assemble with `fig_panel3d_compare.py`.

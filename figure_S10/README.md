# Supplemental figure 10 — Identification of the A119L change-of-efficacy variant

Both panels come from **one script**, which renders the complete 165 × 175 mm figure
rather than panels to be assembled elsewhere.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/fig_119_panel.py` | `mor-dms-modeling/src/` | `data/modeling/scores_long.parquet`, `data/modeling/ligands.yaml` | `panels/figF_position119_all_substitutions.pdf` (+ `.png`) |

Supporting modules, copied alongside: `code/common.py` (paths, the 65–355 position window,
zone definitions) from `mor-dms-modeling/src/`, and `code/plots.py` (the `figure()` /
`save()` millimetre canvas helpers) from `mor_efficacy/src/`.
`code/fig_waterfall_pos.py` is the standalone version of panel **a**, kept for reference —
`fig_119_panel.py` redraws the waterfall internally rather than importing it.

**a** Position-level waterfall ranking receptor positions by the *antagonist-to-agonist
switch* axis: the mean naloxone/naltrexone score minus that variant's own basal
(forskolin-only) score. Position 119 ranks **#1 of 246**. Only positions with at least
**5** scored substitutions enter the ranking, and naloxone and naltrexone must agree to
within **0.05** — drawing the two antagonists separately rather than pooling them makes
that agreement visible as a consistency check.

**b** All 19 substitutions at position 119, each panel plotting the DMS score against
wild-type ligand efficacy (TRUPATH Emax, % DAMGO) with a fitted slope. Panels are ordered
by that slope, most negative first: a steep negative slope is the "every ligand becomes a
fuller agonist" phenotype, strongest at **A119I (−0.130)**, **A119P (−0.125)** and
**A119L (−0.096)**. The basal (no ligand, FSK) point sits on its own broken segment of
the x axis rather than being carried across as a level line, because it is not a ligand
efficacy. Per-substitution surface expression is annotated in grey at the bottom left of
each panel.

## Position window

Analysis is restricted to positions **65–355** — TM1 through helix 8. The N-terminus
(1–64) and C-terminal tail (356–400) are disordered in every μOR structure, carry no
Ballesteros–Weinstein number and have no reliable structural covariates. They are
filtered *before* analysis rather than after, because they change the permutation/FDR
calibration, the conditioning quantiles, the expression LOWESS fit and the zone
denominators. This is the same window used for the ligand-profile clustering in
figure 5b.

## Data

`data/modeling/scores_long.parquet` is the long-format per-variant × per-ligand score
table used by the modeling analyses, and `data/modeling/ligands.yaml` carries each
ligand's TRUPATH Emax and its holdout flag — the x-axis values in panel **b** and the
ligand filtering in panel **a** both come from it.

## Note on paths

`fig_119_panel.py` prepends two absolute `sys.path` entries pointing into the original
analysis repository. Repoint those at `figure_S10/code/` before re-running, since
`common.py` and `plots.py` are copied here.

# SR-17018 TruPath Gi1

SR-17018 is the one ligand in the 15-ligand panel without a full concentration–response
in the TruPath Gi1 assay: its solubility was poor, so it is reported in **Figure 5c** as a
single measurement at **100 µM**.

## The published value

**98.5 ± 19.3 % of DAMGO (n = 6).**

From the **20260502** plate, rows **E / J / O** at 1e-4 M, each normalised to those rows'
own vehicle wells and divided by that plate's DAMGO peak.

| File | Role |
|---|---|
| `sr17018_trupath_wells.csv` | Well-level readings. The six E/J/O values behind the published number |
| `extract_sr17018_trupath.py` | Reads the raw Prism/CLARIOstar files and writes the CSVs here |
| `sr17018_trupath_prism_values.csv`, `..._prism_fits.csv`, `..._dose_summary.csv` | Per-dose values and Prism fits for both runs |
| `emax_pct_damgo.py`, `sr17018_emax_pct_damgo.csv` | Percent-of-DAMGO under several baseline choices |

The value is computed by `figure_05/code/trupath_emax_bar_points.py` (function
`sr_points()`) and carried as a constant into
`figure_05/code/trupath_emax_bar_40mm.py`, which renders the panel.

## Caveats that belong with the number

- One concentration only. There is no saturation, so Emax is not constrained by a fit.
- The 95% CI is wide and is not distinguishable from DAMGO's maximum.
- SR-17018 is therefore placed with the strong agonists in Figure 5c, but on a single
  point rather than a curve. The figure legend states the 100 µM caveat.

## A second run exists and gives a different answer

The earlier **20260117** plates — the same plate set as the published 14-ligand panel —
give **42.4 ± 6.8 %** when the 100 µM point is measured against DAMGO's baseline. That
analysis, and the arguments either way (buprenorphine's control value, DAMGO's assay
window, the donor-channel behaviour in SR-17018's top-dose wells), live in the working
analysis repository under `sr17018_trupath/` and are not reproduced here, since Figure 5c
uses the 20260502 value.

Recorded so the choice is explicit rather than implicit: **the paper reports 98.5%.**

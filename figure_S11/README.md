# Supplemental figure 11 — Double mutation BRET experiments

**a** TRUPATH Gi1 dose-responses testing the addition of A119L to K100D/N, V175N, R278D
and I280A/K, across DAMGO, PZM21, nalbuphine and naloxone. Four traces per panel: WT,
A119L, the partner single, and the A119L + partner double.

| Step | Script | Original path | Reads | Writes |
|---|---|---|---|---|
| 1 | `code/extract_prism_points.py` | `plots/doubles/` | the two source `.prism` projects (not distributed) | `data/pharmacology/doubles/doubles{1,2}_points.csv` |
| 2 | `code/drc_fit.py` | `plots/doubles/` | — | the shared `fit_drc()` sigmoid + F-test |
| 3 | `code/supp_doubles.py` | `plots/doubles/` | `data/pharmacology/doubles/doubles_merged_points.csv` | `panels/supp_doubles_rescue_all.pdf` (also `_series1`, `_series2`) |

`drc_doubles_grid_supp.py` is the earlier per-series grid layout.

Assay conventions are the same as figure 8e–f: signal-down, activation window
`1 − Top = −Span`, the two runs co-scaled and pooled for WT and A119L (n = 10), no
re-normalisation to DAMGO. Fitting keeps the fitted Span (any sign) when the
extra-sum-of-squares F-test against a flat line gives p < 0.05, otherwise the variant is
called "no response" and drawn flat.

**b** Relative surface expression of every variant tested, by flow cytometry (mean
SureLight APC fluorescence, CMV-Flag-μOR). All variants are within two-fold of wild type
(red dotted lines). Drawn in Prism; no repository code.

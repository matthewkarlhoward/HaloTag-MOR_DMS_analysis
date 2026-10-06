# K100 / V175 / R278 TruPath curves, normalised

Panels d-f of the G-protein figure plot **raw** 515/410 BRET. This folder finds
the replicates behind them and puts them on the same normalised scale as the
WT / A119L / double-mutant curves.

## Where the data is

Raw replicates are in two Prism projects on the Desktop, not in the repo:

| source in the CSVs | file (`~/Desktop/mor_trupath/`) | sheets | what it holds |
|---|---|---|---|
| `doubles_raw` (run 20260427) | `20260125_mor_ICL_mutants.prism` | `20260427 data: * R278` (three byte-identical copies) | verbatim raw copy of the 20260427 doubles plate: WT / 119L / 100D / 100N / 278D / 119L+X x DAMGO, PZM21, Nalbuphine, Naloxone, 6 reps. **This is what panels d and f plot**, and `figures/doubles1_points.csv` is its first four replicates after normalisation. |
| `ICL` (run 20260125) | `20260125_mor_ICL_mutants.prism` | `<Ligand> read 2 comp` | WT / K100D / K100N / V175E / V175N x DAMGO, Morphine, Butorphanol, Fentanyl, PZM21, Nalbuphine, 6 reps. **This is what panel e plots** (read 2, not read 1), and an independent K100 replicate set. |
| `R278_run1` (run 20260121) | `20260121_mkh_mkh_mutant_trupath.prism` | `<Ligand> read 2 comp` | WT / R278D / I280D / I280K / I280P x DAMGO, Buprenorphine, Morphine, MP, Fentanyl, Nalbuphine, 6 reps. Earlier independent R278D run, no PZM21. |

V175E exists only on the ICL plate; the doubles plates carry 175N (20260421)
but never 175E.

## The normalisation

Reverse-engineered from the 20260427 plate, whose raw and normalised forms are
both in hand: **each curve (ligand x variant) is divided by its own fitted
no-drug plateau**, i.e. fit

    Y = Bottom + (Top - Bottom) / (1 + 10^(LogEC50 - X))      (Hill = 1)

to that curve's own raw replicates and divide every point by `Bottom` (the
dilute-end plateau; this is a signal-down assay, so it is the larger plateau).
Baseline becomes 1.0 and the agonist response reads downward.

`normalize_singles.py --validate` reproduces all 32 normalised 20260427 curves
from the raw plate: 31/32 to <1e-4, the flat `Nalbuphine_100N` non-responder to
2.7% (its plateau is not identifiable, so Prism's optimiser landed elsewhere).

This is the per-curve scaling for the dose-response panels. It is **not** the
WT-DAMGO anchoring used for the efficacy bars and heatmaps, which is applied
later on top of these normalised spans.

Two deliberate differences from `figures/doubles*_points.csv`:
- all **6** replicates are used (the doubles CSV keeps the first 4), matching
  what panels d/f plot, so the divisors differ from the doubles ones by ~5%;
- the fit spans every real dose (the vehicle well, `logM == 0`, is dropped),
  not a fixed -13..-4 window.

## Pipeline

```bash
python3 plots/singles_trupath/extract_singles_prism.py   # Prism -> figures/singles_raw_points.csv
python3 plots/singles_trupath/normalize_singles.py       # -> figures/singles_norm_points.csv + _params.csv
python3 plots/singles_trupath/normalize_singles.py --validate
python3 plots/singles_trupath/drc_singles_grid.py        # main 3x3 panel
python3 plots/singles_trupath/drc_singles_supp.py        # every ligand each plate ran
```

Curves are re-fit on the normalised points with the doubles F-test rule
(`plots/doubles/drc_fit.py`): a sigmoid is drawn only where it beats a flat line
at p < 0.05, otherwise the trace is drawn flat at the data mean. Fits whose
LogEC50 pins to a dose bound are flagged `pinned` in
`figures/singles_norm_params.csv`.

## Outputs

- `singles_norm_grid.pdf` / `.png` (118 x 100 mm) - the direct normalised
  replacement for panels d-f: rows K100 / V175 / R278, columns DAMGO / PZM21 /
  Nalbuphine. The V175 row carries the ICL plate's own WT, which is a different
  WT from the 20260427 one in rows 1 and 3.
- `singles_supp_icl.pdf`, `singles_supp_r278.pdf` - all ligands each plate ran.
- `singles_norm_grid_params.csv`, `singles_supp_params.csv` - fitted LogEC50,
  Emax, activation window, R^2 and F-test p for every drawn curve.

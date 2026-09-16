# Reproducing the figures

What can be regenerated from this repository, in what order, and what cannot.

## Paths

Every script is copied **verbatim** from the analysis repository, so each one references
its inputs either by an absolute path (`/Users/mkh/GitHub/mor_dms_analysis/...`) or by
walking up from its own location (`Path(__file__).resolve().parents[N]`). Neither
resolves here. Repoint them at `data/` before running. Each figure README lists the
original path of every script together with the inputs it reads, so the mapping is
unambiguous.

This was a deliberate choice: rewriting ~70 scripts' paths without being able to re-run
them all would have risked silent breakage in files that are otherwise exactly what
produced the published panels.

## Order

Almost everything reads `data/dms_scores/composite_dms_scores.csv`. If you are starting
from the counts rather than the deposited tables:

```
data/variant_counts/            532 per-bin count files (the primary data)
  -> Lilace (external)          per-condition variant scores
  -> shared/create_composite_scores.py
                                data/dms_scores/composite_dms_scores.csv
```

Then the four derived tables, each independent of the others:

```
shared/refit_3param_robust.py        -> data/curve_fits/          (morphine, fentanyl)
shared/refit_3param_robust_damgo.py  -> data/curve_fits/          (DAMGO)
shared/01_fit_operational.py         -> data/operational_model/
shared/pca_analysis.py               -> data/pca/                 (position level)
shared/pca_variant_analysis.py       -> data/pca/                 (variant level)
shared/build_per_position_pc_table.py-> data/pca/per_position_pc_scores.csv
shared/compute_lof_gof.R             -> data/dms_scores/lof_gof_scores.csv
```

All of these are already deposited, so figure scripts can be run directly without
recomputing them.

Figure scripts are then independent of each other and can be run in any order, with two
exceptions:

* `figure_08/code/build_merged_class_assets.py` must run before
  `build_class_variant_assets.py` and before the ChimeraX render.
* `figure_S02/code/make_plots.py` writes `per_bin_coverage_summary.csv`, which
  `compact_coverage.py` then reads. The summary is deposited, so panel i can be drawn
  without recounting.

## Structure panels

Figures 3b–d, 4c, 4f, 5g, 6a, 7a–b, 8a–c and S6b–f, S7b–c are ChimeraX renders. The
scripts write `.defattr` attribute files and `.cxc` command scripts; ChimeraX 1.7+ then
has to be driven to produce the image. Rendering must be done in the **GUI** — macOS
`--nogui` cannot render (no OpenGL) and `--offscreen` is Linux-only — and the window must
be `1600 x 1200` to match the saved camera aspect, or framing shifts between runs.

## What cannot be regenerated from this repository

| Panels | Why |
|---|---|
| 1a, 1b, 1d, 2b, 4a, 5a, 8e | Illustrator cartoons |
| 2a, S2d, S2e | FlowJo exports; workspaces and example `.fcs` are in `data/facs_gating/` |
| 1b inset, 2a structures, S1a | ChemDraw |
| 7d, 7e, 7f | Single-mutant TRUPATH BRET, fitted and drawn in Prism. **The underlying points are not in this repository.** |
| S11b | Surface-expression flow cytometry, drawn in Prism. Same. |
| S8 | cryoSPARC processing. The refined model is deposited as `data/structures/raw/experimental/buprenorphine.pdb`. |
| 5g | The PC2 colouring was painted by hand in ChimeraX. `data/pca/per_position_pc_scores.csv` holds the values, but there is no script that writes the attribute file. |

The double-mutant BRET (8d, 8f, 8g, S11a) **is** reproducible — the replicate-level points
were extracted from the Prism projects into
`data/pharmacology/doubles/doubles_merged_points.csv`.

## Verifying the data

`MANIFEST.sha256` lists every file under `data/` with its size and SHA-256.

```bash
python3 - <<'EOF'
import csv, hashlib
bad = 0
for r in csv.DictReader(open("MANIFEST.sha256")):
    h = hashlib.sha256(open(r["path"], "rb").read()).hexdigest()
    if h != r["sha256"]:
        print("MISMATCH", r["path"]); bad += 1
print("ok" if not bad else f"{bad} mismatched")
EOF
```

## Environment

`requirements.txt` (Python) and `r-requirements.txt` (R). Note that **`pyarrow` and
`pyyaml` are required for supplemental figure 10** and were not installed on the machine
that assembled this export, so that script is the one most likely to fail first on a
fresh environment.

The figures were produced over roughly a year and no version is pinned in the code, so
treat the recorded versions as a reference point rather than a guarantee.

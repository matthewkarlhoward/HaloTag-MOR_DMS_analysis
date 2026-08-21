# Supplemental figure 1 — Ligands used in the study and wild-type signalling data

**a** Chemical structures of all 15 opioid ligands. ChemDraw; no code.

**b** Pooled cAMP transcriptional reporter concentration-responses on wild-type μOR.
**c** TRUPATH Gi1 BRET concentration-responses on wild-type μOR.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/camp_trupath_composite.py` | `trupath_vs_camp_efficacy/` | `data/pharmacology/transcriptional_drc.xlsx`, `data/pharmacology/20260117_mor_wt_trupath_compilation.pzfx` | `panels/camp_trupath_composite.pdf` (both panels, shared legend) |
| `code/camp_drc_curves.py` | same | `transcriptional_drc.xlsx` | panel b alone |
| `code/trupath_drc_curves.py` | same | the `.pzfx` | panel c alone |

**d** Table of pharmacologic parameters (logEC50, Emax % DAMGO, n) for both assays.

| Step | Script | Original path | Reads | Writes |
|---|---|---|---|---|
| 1 | `code/build_param_table.py` | `trupath_vs_camp_efficacy/` | `transcriptional_drc.xlsx`, the `.pzfx` | `data/pharmacology/param_table.csv` |
| 2 | `code/render_param_table.py` | same | `param_table.csv` | `panels/param_table.pdf` |

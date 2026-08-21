# Figure 3 — Variant effects on morphine EC50 and Emax

**a** Left: schematic. Right: per-variant EC50 x Emax landscape for morphine.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/morphine_variant_ec50_vs_emax.py` | `curve_refitting/` | `data/curve_fits/refit_3param_robust_morphine.csv` | `panels/morphine_variant_ec50_vs_emax.pdf` |

**b, c, d** Per-position average effects painted on the μOR–morphine structure (PDB 8EF6):
Emax (red), EC50 (blue), and the two overlaid.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/make_chimerax_param_maps.py` | `curve_refitting/` | `data/curve_fits/refit_3param_robust_{morphine,fentanyl,damgo}.csv` | `code/chimerax_maps/attr_*.defattr` + per-view `.cxc` |

`code/chimerax_maps/` holds the generated attribute files and the ChimeraX command
scripts as run (`render_all.cxc` renders every view; `workspace.cxc` sets up the session).
Open the structure from `data/structures/raw/experimental/8ef6.cif`, `open` the
`.defattr`, then run the matching `.cxc`.

Panel **d** is the overlay view: EC50-sensitive positions (blue spheres) and
Emax-sensitive positions (red spheres) on the same model.

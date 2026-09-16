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

Panel **d** has two parts. The structure is the overlay view — EC50-sensitive positions
(blue) and Emax-sensitive positions (red) on the same model, coloured by
`code/build_overlap_cartoon.py` (with `code/build_chimerax_overlap.py` and
`chimerax_overlap/` for the ChimeraX assets). The Venn beside it is:

| Script | Original path | Writes |
|---|---|---|
| `code/venn_2SD.py` | `plots/ec50_emax_overlap/` | `panels/ec50_emax_venn_2SD.pdf` — EC50-only 29, both 27, Emax-only 49 |

Both read `code/overlap_core.py`, which loads
`data/curve_fits/refit_3param_robust_morphine.csv` and
`data/structures/processed/ligand_distances/per_pdb_legacy/8ef6_distances.csv`.

> **Threshold mismatch to check.** The figure legend says positions are coloured when
> their effect "exceeds one standard deviation of the synonymous population", but the
> Venn numbers (29/27/49) come from `venn_2SD.py`, which is the **2 SD** cutoff. At 1 SD
> the counts are 150 EC50 / 143 Emax / 96 both (see
> `code/EC50_EMAX_OVERLAP_README.md`). Either the structure and the Venn use different
> thresholds, or the legend needs correcting.

`code/EC50_EMAX_OVERLAP_README.md` documents the fuller analysis behind this panel, and
is worth reading before describing the two sets as separable: they overlap 2-6x above
chance at every stringency tested.

# Figure 7 — Mutations to intracellular residues have efficacy-dependent effects

**a, b** ChimeraX views of the canonical active state (8EFQ) and the latent/intermediate
state with the Gαi C-terminus bound extrahelically (9PXW). No repository code; both CIFs
are in `data/structures/raw/experimental/`.

**c** Heatmap of G protein-contacting and non-contacting intracellular residues, showing
the fraction of loss-of-function variants (below two SD from the synonymous mean) for each
ligand.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/intracellular_gprotein_heatmap_rotated_pctlof_allpos.R` | `plots/interfaces/intracellular/` | `data/dms_scores/composite_dms_scores.csv`; `data/structures/processed/gprotein_distances/*_gprotein_distances.csv` (9PXW excluded); `data/structures/processed/intermediate_gprotein_distances/intermediate_contacts_combined.csv` | `panels/intracellular_gprotein_heatmap_rotated_pctlof_allpos.pdf` |

The `_allpos` version keeps every intracellular position (ICL1–3 and H8 as well as the
contact set); `intracellular_gprotein_heatmap_rotated_pctlof.R` is the contacts-only
version. The G-protein contact class shown along the y-axis is Active / Intermediate /
Both / NA, derived from which state's first shell a position appears in.

| Script | Original path | Writes |
|---|---|---|
| `code/analyze_gpcr_gprotein_distances.py` | `structures/scripts/gprotein/` | `data/structures/processed/gprotein_distances/` |

**d, e, f** TRUPATH Gi1 BRET dose-responses for K100D/K100N, V175E/V175N and R278D
against DAMGO, PZM21 and nalbuphine. Fit and drawn in Prism; no repository code and the
Prism project is not distributed here.

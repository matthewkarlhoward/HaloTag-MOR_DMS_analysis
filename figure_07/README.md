# Figure 7 — Mutations to intracellular residues have efficacy-dependent effects

**a, b** ChimeraX views of the canonical active state (8EFQ) and the latent/intermediate
state with the Gαi C-terminus bound extrahelically (9PXW). No repository code; both CIFs
are in `data/structures/raw/experimental/`.

**c** Heatmap of G protein-contacting and non-contacting intracellular residues, showing
the fraction of loss-of-function variants (below two SD from the synonymous mean) for each
ligand.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/intracellular_gprotein_heatmap_rotated_pctlof.R` | `plots/interfaces/intracellular/` | `data/dms_scores/composite_dms_scores.csv`; `data/structures/processed/gprotein_distances/*_gprotein_distances.csv` (9PXW excluded); `data/structures/processed/intermediate_gprotein_distances/intermediate_contacts_combined.csv` | `panels/intracellular_gprotein_heatmap_rotated_pctlof.pdf` — **this is the published panel**, 60 x 110 mm, 30 positions |

`..._allpos.R` is a companion that keeps **every** candidate intracellular position with
no LoF filter (60 x 186 mm). It is not in the paper; it is retained because the published
panel's position set is defined by subtracting a filter from it.

**Position selection.** Candidates are the union of: G protein first-shell positions in
the active state, first-shell positions in the latent/intermediate state, and all
remaining ICL1/ICL2/ICL3/H8 positions from the GPCRdb secondary-structure annotation.
The published panel then keeps only positions whose mean missense effect across the 15
ligands falls more than `LOF_THRESHOLD = -0.05` below the synonymous mean, leaving 30.

| Step | Positions |
|---|---|
| G protein first shell, 29 active structures | 35 |
| G protein first shell, latent state (9ODL, 9PXW), Gα only | 18 |
| — active only / both / intermediate only | 27 / 8 / 10 |
| Remaining ICL1–3 and H8 with no G protein contact | 14 |
| **Candidate set** | **59** |
| Retained by the mean-effect filter | **30** |

**Active set (29):** 6DDE, 6DDF, 7SBF, 7SCG, 7T2G, 7T2H, 7U2K, 7U2L, 8Y72, 8Y73, 8EF5,
8EF6, 8EFB, 8EFL, 8EFO, 8EFQ, 8F7Q, 8F7R, 8K9K, 8K9L, 9PXV, 9PXX, 9PXY, 9PY2, 9PY3, 9PY4,
9WST, 9WSW, 9BQJ — Gi1 throughout except 8K9K/8K9L (Gi3) and 9WST/9WSW (Gz).
**Latent set (2):** 9ODL, 9PXW.
Structures with no resolved Gα are excluded: 4DKL, 5C1M, 7UL4, 8E0G, 8QOT, 9MQI, 9PXU,
9BJK, 9WSV, 9WSX.

Recomputed from the deposited data; reproduces the panel's 30 rows exactly. Note the two
distinct thresholds: `LOF_THRESHOLD = -0.05` on the mean effect selects *which positions
appear*, while `LOF_SD = 2` (per ligand) defines the *cell values*. No surface-expression
filter is applied here, unlike the figure 8 network analysis. The G-protein contact class shown along the y-axis is Active / Intermediate /
Both / NA, derived from which state's first shell a position appears in.

| Script | Original path | Writes |
|---|---|---|
| `code/analyze_gpcr_gprotein_distances.py` | `structures/scripts/gprotein/` | `data/structures/processed/gprotein_distances/` |

**d, e, f** TRUPATH Gi1 BRET dose-responses for K100D/K100N, V175E/V175N and R278D
against DAMGO, PZM21 and nalbuphine. Fit and drawn in Prism; no repository code and the
Prism project is not distributed here.

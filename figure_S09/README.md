# Supplemental figure 9 — Functional residue network mapping

**a** Node counts and **b** edge counts in the LOF and GOF networks for each ligand
efficacy class. **c** Top-down and **d** side views of the LOF and GOF nodes, with
marginal kernel densities along each axis; grey dots are positions that are not a
functional node in any class.

| Script | Original path | Reads | Writes |
|---|---|---|---|
| `code/build_class_node_edge_plots.py` | `network_tools/` | `data/structures/ca_coords/pdb_8efq_chain_R_ca.csv`, `code/GPCRdb_OPRM1_table.csv` (in `figure_08/code/`), the node/edge sets from `build_merged_class_assets.py` | `panels/class_node_counts_bar.pdf`, `class_edge_counts_bar.pdf`, `class_node_density_overlay_{lof,gof}.pdf` |

Node and edge definitions are identical to figure 8 — see `figure_08/README.md` for the
Q75 + state-contact node rule, the 8 Å Cα edge rule, and the LOF/GOF filename swap.
`build_merged_class_assets.py` is copied here as well so this figure can be regenerated
standalone.

Counts as published: nodes LOF/GOF = 26/24 (antagonist), 44/4 (weak), 57/4
(intermediate), 56/0 (strong); edges = 27/30, 78/1, 108/1, 104/0.

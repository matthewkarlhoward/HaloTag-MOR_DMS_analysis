## Stacked per-condition DMS effect heatmaps (No Ligand + all ligands),
## rows ordered by hierarchical clustering of each ligand's per-position mean
## missense effect, with the clustering dendrogram drawn on the right.
## Outputs (one PDF per fill scaling, everything else identical):
##   dms_heatmap_stacked_magma_sym.pdf      symmetric about the synonymous midpoint
##                                          (same limits as the RdBu original)
##   dms_heatmap_stacked_magma.pdf          empirical 2-98% range of the effects
##   dms_heatmap_stacked_magma_lof.pdf      LOF-contrast: q05 -> synonymous midpoint
##
## MAGMA VARIANT of dms_stacked_heatmap.R: identical clustering / layout, but the
## fill uses the same "magma_fastyellow" ramp as the 3-ligand dose-response figure
## (dms_dose_response_heatmap.py). There the plotted score is signal-DOWN, so dark
## = high signaling and bright yellow = loss of signaling. `*_effect` here runs the
## other way (low = LOF), so the ramp is REVERSED to keep the visual convention:
## LOF -> bright yellow, neutral -> mid red/purple, GOF -> near black.
##
## A sequential ramp on a neutral-heavy, LOF-skewed distribution is sensitive to
## where the limits sit, so three scalings are rendered (see SCALINGS below).
##
## Clustering uses CORRELATION distance (1 - Pearson) + WARD.D2 linkage: this
## clusters ligands by profile *shape* (which positions are mutation-sensitive),
## not overall effect magnitude. Raw Euclidean instead produced a magnitude
## artifact that stranded DAMGO with {TRV130, Carfentanil}; under a shape distance
## DAMGO's nearest profiles are all strong agonists. The tree is rotated toward an
## efficacy-ladder order, then DAMGO is pinned to the bottom row.
## Note: antagonist / partial / efficacious-agonist resolve cleanly, but the
## intermediate-vs-strong agonist split is NOT a distinct DMS fingerprint.
##
## Layout: the heatmap is ONE faceted plot (facets touching, panel.spacing = 0)
## so a single panel.border supplies both the line between rows and the outer
## box. The dendrogram is a second plot in a plain 2-column patchwork, which
## aligns its panel to the facet grid exactly (a multi-row `design` span does
## NOT align reliably). Dendro y-window [0,16] => leaf j sits at y = 15.5 - j,
## the vertical centre of heatmap row j+1 (row 1 = No Ligand, no leaf).

library(tidyverse)
library(patchwork)
library(dendextend)

OUT  <- "/Users/mkh/GitHub/mor_dms_analysis/plots/heatmaps"
FONT <- "Helvetica"
TXT  <- 6             # all text (pt)
LW05 <- 0.5 / 2.13    # ggplot linewidth that renders ~0.5 pt (border, ticks, dendro)
SZ   <- 1.3           # pharmacology-group circle at each dendrogram leaf tip

## Pharmacological group + colour for the leaf-tip circles (match the legend).
grp     <- c(Naltrexone="antag", Naloxone="antag",
             Nalbuphine="partial", Buprenorphine="partial", Butorphanol="partial", MP="partial",
             TRV130="intermed", PZM21="intermed", Methadone="intermed",
             C6guano="strong", Morphine="strong", SR17018="strong",
             Fentanyl="strong", Carfentanil="strong", DAMGO="strong")
grp_col <- c(antag="#7b3fa0", partial="#3779b9", intermed="#9b1c1d", strong="black")

## "magma_fastyellow" (matplotlib magma resampled at t^0.55, 33 stops), listed
## dark -> bright exactly as in dms_dose_response_heatmap.py.
MAGMA_FASTYELLOW <- c(
  "#000004","#251255","#420f75","#59157e","#6b1d81","#7c2382","#8c2981",
  "#9b2e7f","#a8327d","#b5367a","#c23b75","#cd4071","#d8456c","#e24d66",
  "#e95462","#f05f5e","#f4695c","#f8745c","#fa7d5e","#fc8961","#fd9467",
  "#fe9d6c","#fea772","#feb078","#febb81","#fec287","#fecc8f","#fed597",
  "#fddea0","#fde7a9","#fceeb0","#fcf7b9","#fcfdbf")

## Clustering knobs (see header). DIST: "pearson" | "spearman" | "euclidean" |
## "euclidean_z" (row z-scored).  LINKAGE: any hclust method.
DIST    <- "pearson"
LINKAGE <- "ward.D2"
POS_LO  <- 65         # position window, used for BOTH clustering and the heatmap
POS_HI  <- 355
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

dms_scores <- read_csv(
  "/Users/mkh/GitHub/mor_dms_analysis/dms_scores/composite_dms_scores.csv",
  show_col_types = FALSE)

aa_order <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

## ---- Conditions -----------------------------------------------------------
## FSK is the "No Ligand" baseline (pinned to the top, excluded from clustering).
no_ligand <- "FSK"
ligands   <- c("Naltrexone","Naloxone","Nalbuphine","Buprenorphine","Butorphanol",
               "MP","TRV130","PZM21","Methadone","C6guano","Fentanyl","Carfentanil",
               "Morphine","SR17018","DAMGO")               # Emax efficacy order (rotation target)
conditions <- c(no_ligand, ligands)
eff_cols   <- paste0(conditions, "_effect")

## display labels (left of each row)
disp <- c(FSK = "No Ligand", C6guano = "C6-Guano", SR17018 = "SR-17018")
label_of <- function(x) ifelse(x %in% names(disp), disp[x], x)

## ---- Hierarchical clustering of ligands (missense only) -------------------
## Feature per ligand = vector of per-position MEAN MISSENSE effect.
mis_mat <- dms_scores %>%
  filter(type == "missense", position >= POS_LO, position <= POS_HI) %>%
  select(position, all_of(paste0(ligands, "_effect"))) %>%
  group_by(position) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

mat <- t(as.matrix(mis_mat[, -1]))
rownames(mat) <- ligands
mat <- mat[, colSums(is.na(mat)) == 0]              # positions present for all ligands

dmat <- switch(DIST,
  euclidean   = dist(mat),
  euclidean_z = dist(t(scale(t(mat)))),             # row z-score, then Euclidean
  pearson     = as.dist(1 - cor(t(mat))),           # shape, magnitude-invariant
  spearman    = as.dist(1 - cor(t(mat), method = "spearman")),
  stop("unknown DIST"))
hc    <- hclust(dmat, method = LINKAGE)
dend0 <- as.dendrogram(hc)                           # rotate the FRESH tree each time
## Rotate branches toward the Emax efficacy order (clean tree, no crossings).
## ward.D2 locks a few within-group orders, so this is the CLOSEST rotation, not
## an exact efficacy match (Buprenorphine/Nalbuphine, Methadone, Morphine stay put).
clust_order <- labels(rotate(dend0, ligands))
if (tail(clust_order, 1) != "DAMGO")                 # pin DAMGO to the very bottom row
  clust_order <- labels(rotate(dend0, c(setdiff(clust_order, "DAMGO"), "DAMGO")))
row_order   <- c(no_ligand, clust_order)             # final order; No Ligand on top
message("Row order: ", paste(label_of(row_order), collapse = " | "))

## Single stacked panel: each condition is a BH-row block, separated by GAP_Y
## empty rows so the white gap is real space (no tiles covered). One panel =>
## panel.border is the OUTER box only.
disp_levels  <- label_of(row_order)                 # 16 labels, top -> bottom
NBLOCK       <- length(disp_levels)                 # 16 conditions
BH           <- length(aa_order)                    # 20 substitutions / block
GAP_Y        <- 1.8                                 # empty space between blocks (~0.5pt)
PITCH        <- BH + GAP_Y                          # block-to-block stride
block_center <- (NBLOCK - seq_len(NBLOCK)) * PITCH + (BH + 1) / 2
Y_LIM        <- c(0.5, NBLOCK * BH + (NBLOCK - 1) * GAP_Y + 0.5)

## ---- Long data for the heatmap fill --------------------------------------
dms_long <- dms_scores %>%
  select(position, type, mutation, all_of(eff_cols)) %>%
  pivot_longer(all_of(eff_cols), names_to = "condition", values_to = "effect") %>%
  mutate(condition = str_remove(condition, "_effect")) %>%
  filter(mutation %in% aa_order, position >= POS_LO, position <= POS_HI) %>%
  mutate(
    mutation = factor(mutation, levels = rev(aa_order)),
    cond_i   = match(label_of(condition), disp_levels),                    # 1 = top
    y        = (NBLOCK - cond_i) * PITCH + as.integer(mutation)            # blocks + gaps
  )

## Divergent fill centred on the mean synonymous effect, shared across panels.
syn_baseline <- dms_scores %>%
  filter(type == "synonymous") %>%
  summarise(across(all_of(eff_cols), ~ mean(.x, na.rm = TRUE))) %>%
  unlist()
midpoint_val <- mean(syn_baseline, na.rm = TRUE)
max_abs      <- quantile(abs(dms_long$effect - midpoint_val), 0.98, na.rm = TRUE)
q <- function(p) unname(quantile(dms_long$effect, p, na.rm = TRUE))

## Fill scalings, all with oob = squish.
##  sym   - the original RdBu limits. Faithful to the divergent version, but the
##          GOF tail is far shorter than the LOF tail, so the dark half of the
##          ramp goes unused and the bulk piles into the flat salmon middle.
##  range - empirical 2-98%: the full ramp is spent on the observed effects.
##  lof   - LOF contrast. The window is pushed entirely onto the LOF tail
##          (2nd -> 20th percentile of effects): everything from mildly-LOF
##          upwards - neutral, functional, and all GOF - saturates to black, and
##          only the worst tiles have any colour at all. Dead positions (219,
##          167, ...) then read as bright columns on black, and rows with little
##          LOF (the antagonists) go essentially dark. Cutting closer to the
##          baseline (hi = median or 0) instead speckles the whole map yellow,
##          because per-tile noise is the same size as the LOF signal; a gamma
##          warp of a wider window washes out to flat purple rather than fixing
##          that. GOF is discarded here by construction.
SCALINGS <- list(
  sym   = list(lo = midpoint_val - max_abs, hi = midpoint_val + max_abs,
               file = "dms_heatmap_stacked_magma_sym.pdf"),
  range = list(lo = q(0.02),                hi = q(0.98),
               file = "dms_heatmap_stacked_magma.pdf"),
  lof   = list(lo = q(0.02),                hi = q(0.20),
               file = "dms_heatmap_stacked_magma_lof.pdf")
)

X_LIM    <- c(min(dms_long$position) - 0.5, max(dms_long$position) + 0.5)
X_BREAKS <- c(100, 200, 300)

## ---- Heatmap: ONE stacked panel (no facets => no lines between rows) -----
## Condition labels sit at each block centre; panel.border = outer box only.
make_heat <- function(lo, hi, gamma = 1) ggplot(dms_long, aes(position, y, fill = effect)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = rev(MAGMA_FASTYELLOW),                 # LOF = yellow, GOF = black
    ## gamma > 1 pulls the bright end into the bottom of the window
    values = seq(0, 1, length.out = length(MAGMA_FASTYELLOW))^gamma,
    limits = c(lo, hi),
    oob = scales::squish, na.value = "white") +
  scale_x_continuous(limits = X_LIM, breaks = X_BREAKS, expand = c(0, 0)) +
  scale_y_continuous(breaks = block_center, labels = disp_levels,
                     limits = Y_LIM, expand = c(0, 0)) +
  labs(x = "Position") +
  theme_classic(base_size = TXT, base_family = FONT) +
  theme(
    legend.position     = "none",
    panel.border        = element_rect(fill = NA, colour = "black", linewidth = 2 * LW05),
    axis.title.y        = element_blank(),
    axis.text.y         = element_text(size = TXT, colour = "black", hjust = 1,
                                        family = FONT, margin = margin(r = 1.5)),
    axis.ticks.y        = element_blank(),
    axis.line           = element_blank(),
    axis.ticks.x        = element_line(linewidth = LW05, colour = "black"),
    axis.ticks.length.x = unit(1.5, "pt"),
    axis.text.x         = element_text(size = TXT, colour = "black", margin = margin(t = 1)),
    axis.title.x        = element_text(size = TXT, colour = "black", margin = margin(t = 1)),
    plot.margin         = margin(0, 0, 0, 0)
  )

## ---- Dendrogram (right), leaves pinned to block centres -----------------
## ligand j (clust order) is condition i = j+1, whose block centre is leaf_y.
leaf_y <- setNames(block_center[seq_along(clust_order) + 1], clust_order)

node_y <- numeric(nrow(hc$merge))
pos_of <- function(k) if (k < 0) c(0, leaf_y[[hc$labels[-k]]]) else c(hc$height[k], node_y[k])
segs <- list()
for (m in seq_len(nrow(hc$merge))) {
  a <- pos_of(hc$merge[m, 1]); b <- pos_of(hc$merge[m, 2]); xm <- hc$height[m]
  node_y[m] <- (a[2] + b[2]) / 2
  segs[[length(segs) + 1]] <- data.frame(x = a[1], xend = xm, y = a[2], yend = a[2])
  segs[[length(segs) + 1]] <- data.frame(x = b[1], xend = xm, y = b[2], yend = b[2])
  segs[[length(segs) + 1]] <- data.frame(x = xm,   xend = xm, y = a[2], yend = b[2])
}
seg_df <- bind_rows(segs)

## push the whole tree right so its leaf branches clear the leaf circles at x = 0
DEND_SHIFT  <- max(hc$height) * 0.12
seg_df$x    <- seg_df$x + DEND_SHIFT
seg_df$xend <- seg_df$xend + DEND_SHIFT

## group-coloured circle at each leaf tip (x = 0), aligned to its ligand row
leaf_pts <- tibble(x = 0, y = leaf_y[clust_order],
                   col = unname(grp_col[grp[clust_order]]))

dendro <- ggplot(seg_df) +
  geom_segment(aes(x, y, xend = xend, yend = yend), linewidth = LW05, colour = "black",
               lineend = "butt") +
  geom_point(data = leaf_pts, aes(x, y, colour = col), size = SZ) +
  scale_colour_identity() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +   # leaves flush at x=0
  scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +        # match heatmap panel
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 1, unit = "mm"))         # 1 mm gap from heatmap

## ---- Assemble: heatmap | dendrogram (plain 2-column => reliable align) ----
## Overall 1pt margin only (the patchwork default is 5.5pt); just enough that the
## 0.5pt box border isn't half-clipped at the canvas edge.
for (nm in names(SCALINGS)) {
  sc <- SCALINGS[[nm]]
  final_plot <- make_heat(sc$lo, sc$hi, sc$gamma %||% 1) + dendro +
    plot_layout(widths = c(8, 2)) +
    plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1)))
  ggsave(file.path(OUT, sc$file), final_plot,
         width = 72, height = 40, units = "mm", limitsize = FALSE)
  message(sprintf("Saved: %-38s fill limits [%.3f, %.3f]  gamma %.1f",
                  sc$file, sc$lo, sc$hi, sc$gamma %||% 1))
}

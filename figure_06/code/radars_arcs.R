#!/usr/bin/env Rscript
# =============================================================================
# radars_arcs.R  —  per-class LOF circular barplots ("radars") with structural arcs
# -----------------------------------------------------------------------------
# One panel per efficacy class: circular barplot of LOF burden at the 26 filtered
# orthosteric key positions, with
#   - light->dark gradient bars (class colour)
#   - radial position-number labels + per-position mini tick marks
#   - structural-region arcs (TM-helix / ECL-loop bands) with end caps + labels
# The figure fills a 41.25 mm square: a negative plot margin defeats coord_polar's
# built-in dead border (it always insets the circle to radius 0.4 of the panel).
#
# Outputs (NEW files; the original radar scripts/outputs are untouched):
#   circular_lof_strong_filt_arcs.pdf
#   circular_lof_intermediate_filt_arcs.pdf
#   circular_lof_weak_filt_arcs.pdf
#   circular_lof_antagonist_filt_arcs.pdf
#
# LOF = % missense variants with effect < (per-drug synonymous mean - 1 SD),
# averaged across the drugs in each class. "Carfentanil" spelled correctly, so
# the strong class includes it (unlike radars.Rmd's lowercase "carfentanil").
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(colorspace)
  library(scales)
})

setwd("/Users/mkh/GitHub/mor_dms_analysis/plots/radars")
dms_scores <- read_csv(
  "/Users/mkh/GitHub/mor_dms_analysis/dms_scores/composite_dms_scores.csv",
  show_col_types = FALSE)

# ── Effect columns (standard: drop SE + surface) ──────────────────────────────
alldrug_cols <- colnames(dms_scores)[grepl("effect", colnames(dms_scores))]
alldrug_cols <- alldrug_cols[!grepl("_se",                alldrug_cols)]
alldrug_cols <- alldrug_cols[!grepl("Surface_effect",     alldrug_cols)]
alldrug_cols <- alldrug_cols[!grepl("Surface_pos_effect", alldrug_cols)]

# ── Efficacy classes: ligands + colour (matches the class-average radars) ──────
# arc_col = base colour of the structural arcs + end caps + region labels (black on
# the black strong panel; grey to match the axis on the coloured panels).
# hl = region(s) drawn in black instead, to highlight them on that panel.
classes <- list(
  strong       = list(ligands = c("DAMGO", "Carfentanil", "Fentanyl",
                                  "Morphine", "C6guano"),       col = "#231f20", arc_col = "black", hl = character(0)),
  intermediate = list(ligands = c("Methadone", "TRV130", "PZM21"),
                                                                col = "#9b1c1d", arc_col = "grey",  hl = character(0)),
  weak         = list(ligands = c("MP", "Buprenorphine",
                                  "Nalbuphine", "Butorphanol"), col = "#3779b9", arc_col = "grey",  hl = "TM3"),
  antagonist   = list(ligands = c("Naloxone", "Naltrexone"),    col = "#6f08a3", arc_col = "grey",  hl = "TM3")
)
all_ligands <- unlist(map(classes, "ligands"), use.names = FALSE)

# ── Key positions, filtered set (drops 129/130/220/235) ───────────────────────
key_positions  <- c(77, 116, 119, 122, 123, 126, 129, 130, 135, 145, 146, 149,
                    150, 152, 153, 156, 220, 235, 238, 295, 298, 299, 302, 320,
                    321, 324, 327, 328, 330, 331)
filt_positions <- setdiff(key_positions, c(129, 130, 220, 235))

# ── Structural region per residue (human OPRM1 TM / loop boundaries) ──────────
region_of <- function(p) dplyr::case_when(
  p %in% c(77)                          ~ "TM1",
  p %in% c(116, 119, 122, 123, 126, 129, 130) ~ "TM2",
  p %in% c(135)                         ~ "ECL1",
  p %in% c(145, 146, 149, 150, 152, 153, 156) ~ "TM3",
  p %in% c(220, 235, 238)               ~ "TM5",
  p %in% c(295, 298, 299, 302)          ~ "TM6",
  p %in% c(320, 321, 324, 327, 328, 330, 331) ~ "TM7",
  TRUE ~ NA_character_)

# ── House style: black, Helvetica, 6pt text, 0.5pt lines ─────────────────────
.PT  <- 72.27 / 25.4            # ggplot's points-per-mm constant (.pt)
gg   <- function(pt) pt / .PT   # point size -> ggplot mm units
TXT  <- 6                       # text point size
LWD  <- 0.5                     # line weight (pt)
FONT <- "Helvetica"

# ── Radii (data-y space; bars run 0..100) ─────────────────────────────────────
Y_HOLE <- -8      # inner hollow (lower data-y limit -> centre)
R_NUM  <- 111     # position-number labels
R_ARC  <- 156     # structural-region arcs (clear of the radial number text)
R_LAB  <- 183     # region labels (gap so side TM labels clear the arcs)
Y_MAX  <- 208     # outer data-y limit (label text reaches the canvas edge)
PAD    <- 0.42    # arc half-extension past end bars (~half a bar width)
MARG   <- -14     # plot margin (pt); negative defeats coord_polar's dead border

# ── % LOF per (drug, position) over the filtered set ──────────────────────────
dms_long <- dms_scores %>%
  filter(type %in% c("missense", "synonymous")) %>%
  pivot_longer(all_of(alldrug_cols), names_to = "drug", values_to = "effect") %>%
  mutate(drug = str_remove(drug, "_effect$"))

syn_stats <- dms_long %>%
  filter(type == "synonymous") %>%
  group_by(drug) %>%
  summarise(syn_mean = mean(effect, na.rm = TRUE),
            syn_sd   = sd(effect,   na.rm = TRUE), .groups = "drop")

pct_by_drug <- dms_long %>%
  filter(type == "missense", position %in% filt_positions, drug %in% all_ligands) %>%
  left_join(syn_stats, by = "drug") %>%
  mutate(is_lof = effect < (syn_mean - syn_sd)) %>%
  group_by(drug, position) %>%
  summarise(pct_lof = mean(is_lof, na.rm = TRUE) * 100, .groups = "drop")

# ── Position skeleton (shared by every panel: positions are identical) ─────────
n      <- length(filt_positions)
step   <- 360 / n
pos_df <- tibble(position = sort(filt_positions)) %>%
  mutate(idx    = row_number(),
         region = region_of(position),
         theta  = (idx - 0.5) * step,
         flip   = theta > 90 & theta < 270,
         angle  = if_else(flip, theta + 180, theta),
         hjust  = if_else(flip, 1, 0))

region_df <- pos_df %>%
  group_by(region) %>%
  summarise(i0 = min(idx), i1 = max(idx),
            mid = (min(idx) + max(idx)) / 2, .groups = "drop")

# Smooth arcs: many points at constant radius across each region's span
arc_df <- region_df %>%
  mutate(pts = pmap(list(region, i0, i1), function(rg, a, b)
    tibble(region = rg, x = seq(a - PAD, b + PAD, length.out = 80), y = R_ARC))) %>%
  pull(pts) %>% bind_rows()

# Arc end-caps (short radial segments at each arc end)
cap_df <- bind_rows(
  region_df %>% transmute(region, x = i0 - PAD),
  region_df %>% transmute(region, x = i1 + PAD)
) %>% mutate(y = R_ARC - 4, yend = R_ARC + 4)

rings <- data.frame(y = c(25, 50, 75, 100))

# ── Build + save one class panel ──────────────────────────────────────────────
make_radar <- function(name, ligands, col, arc_col, hl) {

  dat <- pct_by_drug %>%
    filter(drug %in% ligands) %>%
    group_by(position) %>%
    summarise(value = mean(pct_lof, na.rm = TRUE), .groups = "drop") %>%
    left_join(select(pos_df, position, idx), by = "position") %>%
    arrange(idx)

  # per-region colour: highlight (hl) regions in black, the rest in arc_col
  region_col <- function(rg) if_else(rg %in% hl, "black", arc_col)
  arc_d <- mutate(arc_df,    rcol = region_col(region))
  cap_d <- mutate(cap_df,    rcol = region_col(region))
  reg_d <- mutate(region_df, rcol = region_col(region))

  p <- ggplot(dat) +
    geom_hline(aes(yintercept = y), data = rings, color = "grey", linewidth = gg(LWD)) +
    geom_col(aes(x = idx, y = value, fill = value), width = 0.85) +
    scale_fill_gradientn(
      colours = c(lighten(col, 0.80), lighten(col, 0.50),
                  lighten(col, 0.20), col),
      values  = rescale(c(0, 20, 60, 100)),
      limits  = c(0, 100), guide = "none") +
    # per-position mini tick marks (restored), at the bar-field edge — grey to
    # match the radial reference rings (the radar "y-axis")
    geom_segment(data = pos_df, aes(x = idx, xend = idx, y = 100, yend = 108),
                 linewidth = gg(LWD), color = "grey") +
    # position numbers
    geom_text(data = pos_df,
              aes(x = idx, y = R_NUM, label = position, angle = angle, hjust = hjust),
              size = gg(TXT), family = FONT, color = "black") +
    # structural-region arcs + end caps + labels — coloured per region (rcol):
    # arc_col base, with `hl` regions (e.g. TM3) in black to highlight them
    geom_path(data = arc_d, aes(x = x, y = y, group = region, color = rcol),
              linewidth = gg(LWD), lineend = "round") +
    geom_segment(data = cap_d, aes(x = x, xend = x, y = y, yend = yend, color = rcol),
                 linewidth = gg(LWD)) +
    geom_text(data = reg_d, aes(x = mid, y = R_LAB, label = region, color = rcol),
              size = gg(TXT), family = FONT) +
    scale_color_identity(guide = "none") +
    scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(Y_HOLE, Y_MAX), expand = c(0, 0)) +
    coord_polar(theta = "x", start = -pi / 2, direction = -1, clip = "off") +
    theme_void(base_family = FONT, base_size = TXT) +
    theme(plot.margin = margin(MARG, MARG, MARG, MARG))

  # 41.25 mm square -> 117 pt (= 1.625 in = 41.275 mm); pdf() floors MediaBox to
  # integer pt, so 117 pt is the nearest to nominal. (SVG for a bit-exact canvas.)
  ggsave(sprintf("circular_lof_%s_filt_arcs.pdf", name), p,
         width = 1.625, height = 1.625, units = "in", bg = "white")
  message("Saved: circular_lof_", name, "_filt_arcs.pdf")
}

walk(names(classes), function(nm) {
  cl <- classes[[nm]]
  make_radar(nm, cl$ligands, cl$col, cl$arc_col, cl$hl)
})

message("Done.")

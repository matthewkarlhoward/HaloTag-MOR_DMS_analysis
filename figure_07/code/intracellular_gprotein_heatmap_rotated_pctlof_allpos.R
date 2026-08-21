## intracellular_gprotein_heatmap_rotated_pctlof_allpos.R  (SUPPLEMENTAL — all positions)
##
## Same as intracellular_gprotein_heatmap_rotated_pctlof.R (rotated, % LOF at 2 SD) but
## WITHOUT the LoF filter: shows ALL ~59 candidate positions (every Gα contact + other
## intracellular residue), not just the 30 in the main-text figure. Taller to fit the
## extra rows; width and cell size unchanged.
##
## % LOF variant of the rotated intracellular Gα heatmap. Same layout, positions,
## axes, labels, and 60×110mm geometry as intracellular_gprotein_heatmap_rotated.R
## — the ONLY difference is the heatmap fill metric:
##   value = % of missense variants at each position (per ligand) that are
##           loss-of-function, i.e. effect more than LOF_SD (=2) SD below THAT
##           ligand's synonymous mean (per-ligand mean − LOF_SD·SD threshold).
## Coloured with a sequential white→red scale (0% → max%). The position set/order
## still come from the mean-effect LoF filter, so it matches the companion figures.
## Keep the mean-effect versions (intracellular_gprotein_heatmap*.R) alongside this.

source("/Users/mkh/GitHub/mor_dms_analysis/howard_theme.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
grDevices::pdf.options(encoding = "ISOLatin1")

OUT_DIR <- "/Users/mkh/GitHub/mor_dms_analysis/plots/interfaces/intracellular"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
save_pdf <- function(plot, filename, width_mm, height_mm) {
  ggsave(file.path(OUT_DIR, filename), plot, width = width_mm, height = height_mm,
         units = "mm", device = "pdf", bg = "white")
  message("Saved: ", file.path(OUT_DIR, filename))
}

DIST_DIR <- "/Users/mkh/GitHub/mor_dms_analysis/structures/processed/gprotein_distances"
INT_FILE <- file.path("/Users/mkh/GitHub/mor_dms_analysis/structures/processed",
                      "intermediate_gprotein_distances/intermediate_contacts_combined.csv")

# ── Contacts + categories (identical to the mean-effect versions) ─────────────
gp_files <- list.files(DIST_DIR, pattern = "_gprotein_distances\\.csv$", full.names = TRUE)
gp_files <- gp_files[!grepl("9PXW", gp_files)]
active_pos <- map_dfr(gp_files, ~ read_csv(.x, show_col_types = FALSE)) %>%
  filter(shell == "first_shell", !is.na(residue_number_human)) %>%
  pull(residue_number_human) %>% as.integer() %>% unique()
inter_pos <- read_csv(INT_FILE, show_col_types = FALSE) %>%
  filter(shell == "first_shell", partner_type == "galpha", !is.na(residue_number_human)) %>%
  pull(residue_number_human) %>% as.integer() %>% unique()
dms <- read_csv("/Users/mkh/GitHub/mor_dms_analysis/dms_scores/composite_dms_scores.csv",
                show_col_types = FALSE)
pos_meta <- dms %>% distinct(position, wildtype, SSE, GPCRdb)
icl_h8 <- pos_meta %>% filter(SSE %in% c("ICL1", "ICL2", "ICL3", "H8")) %>%
  pull(position) %>% as.integer()

both_pos    <- intersect(active_pos, inter_pos)
active_only <- setdiff(active_pos, inter_pos)
inter_only  <- setdiff(inter_pos, active_pos)
other_ic    <- setdiff(icl_h8, union(active_pos, inter_pos))
cat_levels  <- c("Active", "Intermediate", "Both", "NA")
col_meta <- bind_rows(
  tibble(position = active_only, category = "Active"),
  tibble(position = both_pos,    category = "Both"),
  tibble(position = inter_only,  category = "Intermediate"),
  tibble(position = other_ic,    category = "NA")
) %>%
  left_join(pos_meta, by = "position") %>%
  mutate(category = factor(category, levels = cat_levels)) %>%
  arrange(category, position) %>%
  mutate(rec_label = paste0(toupper(wildtype), position))

no_ligand    <- "FSK"
ligand_order <- c("Naltrexone", "Naloxone",
                  "Nalbuphine", "Buprenorphine", "Butorphanol", "MP",
                  "TRV130", "PZM21", "Methadone",
                  "C6guano", "Morphine", "SR17018", "Fentanyl", "Carfentanil", "DAMGO")
ligand_order  <- ligand_order[paste0(ligand_order, "_effect") %in% names(dms)]
display_order <- c(no_ligand, ligand_order)            # left -> right

lig_grp   <- c(Naltrexone = "antag", Naloxone = "antag",
               Nalbuphine = "partial", Buprenorphine = "partial",
               Butorphanol = "partial", MP = "partial",
               TRV130 = "intermed", PZM21 = "intermed", Methadone = "intermed",
               C6guano = "strong", Morphine = "strong", SR17018 = "strong",
               Fentanyl = "strong", Carfentanil = "strong", DAMGO = "strong")
grp_col   <- c(antag = "#7b3fa0", partial = "#3779b9", intermed = "#9b1c1d", strong = "black")
grp_label <- c(antag = "Antagonist", partial = "Weak",
               intermed = "Intermediate", strong = "Strong")
eff_levels <- c("Antagonist", "Weak", "Intermediate", "Strong")

disp        <- c(FSK = "No Ligand", C6guano = "C6-Guano", SR17018 = "SR-17018")
label_of    <- function(x) ifelse(x %in% names(disp), unname(disp[x]), x)
disp_levels <- label_of(display_order)

# ── Position selection: same mean-effect LoF filter as the companion figures ──
mis <- dms %>% filter(type == "missense", position %in% col_meta$position)
mean_eff <- map_dfr(display_order, function(d) {
  col <- paste0(d, "_effect")
  mis %>% group_by(position) %>%
    summarise(value = mean(.data[[col]], na.rm = TRUE), .groups = "drop") %>%
    mutate(ligand = d)
})
syn_mid <- dms %>% filter(type == "synonymous") %>%
  summarise(m = mean(unlist(across(all_of(paste0(ligand_order, "_effect")),
                                   ~ mean(.x, na.rm = TRUE))))) %>% pull(m)
# Supplemental: ALL candidate positions — no LoF filter (cf. main figure's 30).
keep_pos <- col_meta$position
message("All positions (supplemental): ", length(keep_pos), " shown.")

pos_order  <- sort(keep_pos)
lab_levels <- col_meta$rec_label[match(pos_order, col_meta$position)]   # K100 .. F349
y_levels   <- rev(lab_levels)                                          # K100 at TOP
gpcrdb_raw  <- col_meta$GPCRdb[match(pos_order, col_meta$position)]
gpcrdb_list <- lapply(gpcrdb_raw, function(g) {
  if (is.na(g) || !grepl("x", g)) return("")        # loop residues lack a GPCRdb number
  p <- strsplit(g, "x")[[1]]
  bquote(.(as.integer(p[1])) %*% .(as.integer(p[2])))
})
gpcrdb_y <- as.expression(rev(gpcrdb_list))

# ── Metric: % LOF missense variants per position per ligand ───────────────────
# LOF = effect > LOF_SD SD below that ligand's synonymous mean.
LOF_SD  <- 2
syn_df  <- dms %>% filter(type == "synonymous")
lof_thr <- sapply(display_order, function(d) {
  v <- syn_df[[paste0(d, "_effect")]]
  mean(v, na.rm = TRUE) - LOF_SD * sd(v, na.rm = TRUE)
})
mis_keep <- mis %>% filter(position %in% keep_pos)
pct_lof <- map_dfr(display_order, function(d) {
  col <- paste0(d, "_effect"); thr <- lof_thr[[d]]
  mis_keep %>% group_by(position) %>%
    summarise(value = 100 * sum(.data[[col]] < thr, na.rm = TRUE) /
                            sum(!is.na(.data[[col]])),
              .groups = "drop") %>%
    mutate(ligand = d)
})
message(sprintf("%% LOF range: %.1f to %.1f", min(pct_lof$value, na.rm = TRUE),
                max(pct_lof$value, na.rm = TRUE)))

plot_df <- pct_lof %>%
  left_join(col_meta, by = "position") %>%
  mutate(lig_disp  = factor(label_of(ligand), levels = disp_levels),  # x: No Ligand left
         rec_label = factor(rec_label, levels = y_levels))            # y: K100 top
strip_df <- col_meta %>% mutate(rec_label = factor(rec_label, levels = y_levels))

class_df <- tibble(cond = display_order) %>%
  mutate(lig_disp = factor(label_of(cond), levels = disp_levels),
         grp      = ifelse(cond %in% names(lig_grp), unname(lig_grp[cond]), "none"),
         klass    = factor(ifelse(grp == "none", NA, unname(grp_label[grp])),
                           levels = eff_levels))
class_pts <- dplyr::filter(class_df, !is.na(klass))
class_nl  <- dplyr::filter(class_df,  is.na(klass))
max_pct   <- max(plot_df$value, na.rm = TRUE)

# ── Plots (transposed) ────────────────────────────────────────────────────────
lw  <- 0.5 / 2.835
txt <- element_text(family = "Helvetica", size = 6, color = "black")
cat_colors <- c("Active" = "#1b7837", "Intermediate" = "#762a83",
                "Both" = "#e08214", "NA" = "#999999")

# Top row: Efficacy circles (one per ligand column) + ligand names on TOP
p_class <- ggplot(class_df, aes(x = lig_disp, y = 1)) +
  geom_point(data = class_pts, aes(color = klass), size = 1.4) +
  geom_point(data = class_nl, shape = 4, stroke = 0.5, size = 1.2, color = "black") +
  scale_color_manual(values = setNames(unname(grp_col), unname(grp_label)),
                     name = "Efficacy", breaks = eff_levels, drop = FALSE,
                     guide = guide_legend(title.position = "top", ncol = 1)) +
  scale_x_discrete(position = "top", limits = disp_levels,
                   expand = expansion(add = 0.5)) +
  scale_y_continuous(breaks = NULL, limits = c(0.5, 1.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(
    text            = txt,
    axis.text.x.top = element_text(family = "Helvetica", size = 6, color = "black",
                                   angle = 45, hjust = 0, vjust = 0),
    axis.text.y     = element_blank(),
    axis.ticks      = element_blank(),
    axis.line       = element_blank(),
    plot.margin     = margin(0, 0, 0, 0),
    legend.position = "bottom")

# Left column: G protein contact (one tile per position row) + residue labels LEFT
p_strip <- ggplot(strip_df, aes(x = "Contact", y = rec_label, fill = category)) +
  geom_tile(color = "white", linewidth = 0.12) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(limits = y_levels, expand = c(0, 0)) +
  scale_fill_manual(values = cat_colors, name = "G protein contact",
                    guide = guide_legend(title.position = "top", ncol = 1)) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(
    text         = txt,
    axis.text.y  = element_text(family = "Helvetica", size = 6, color = "black"),  # residue
    axis.text.x  = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    plot.margin  = margin(0, 0, 0, 1),
    legend.position = "bottom")

# Heatmap: positions (y, GPCRdb on RIGHT) × ligands (x); fill = % LOF variants
p_heat <- ggplot(plot_df, aes(x = lig_disp, y = rec_label, fill = value)) +
  geom_tile(color = "white", linewidth = 0.12) +
  scale_x_discrete(limits = disp_levels, expand = c(0, 0)) +
  scale_y_discrete(position = "right", limits = y_levels, expand = c(0, 0), labels = gpcrdb_y) +
  scale_fill_gradient(
    low = "white", high = "#9b1c1d",
    limits = c(0, max_pct), breaks = c(0, 50, 100), na.value = "#dddddd", name = "% LOF",
    guide = guide_colorbar(title.position = "top", barwidth = unit(0.30, "cm"),
                           barheight = unit(0.85, "cm"))
  ) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(
    text              = txt,
    axis.text.x       = element_blank(),
    axis.text.y.right = element_text(family = "Helvetica", size = 6, color = "black"),  # GPCRdb
    axis.ticks.y.right = element_line(color = "black", linewidth = lw),   # 0.5pt
    axis.ticks.x      = element_blank(),
    axis.line         = element_blank(),
    panel.border      = element_rect(fill = NA, color = "black", linewidth = lw),  # 0.5pt box
    plot.margin       = margin(0, 0, 0, 0),
    legend.position   = "bottom")

# ── Assemble: circles (top-right) | strip (bottom-left) | heat (bottom-right) ──
design <- "
#A
BC
"
combined <- p_class + p_strip + p_heat +
  plot_layout(design = design, widths = c(1, 16), heights = c(1, length(pos_order)),
              guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5))) &
  theme(legend.position    = "bottom",
        legend.box         = "horizontal",
        legend.direction   = "vertical",
        legend.box.just    = "top",
        legend.title       = element_text(family = "Helvetica", size = 6, color = "black"),
        legend.text        = txt,
        legend.key.size    = unit(0.26, "cm"),
        legend.box.spacing = unit(2, "pt"),
        legend.spacing.x   = unit(0.10, "cm"),
        legend.margin      = margin(0, 0, 0, 0))
save_pdf(combined, "intracellular_gprotein_heatmap_rotated_pctlof_allpos.pdf", width_mm = 60, height_mm = 186)

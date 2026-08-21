# ============================================================================
# network_analysis_shared.R
# ----------------------------------------------------------------------------
# Shared data-prep + helpers for every plot_network_*.R script in this folder.
# Sourced by them automatically — you normally don't run this directly.
#
# Inputs (must sit in the working directory):
#   - composite_dms_scores_R.csv    (override via env var MOR_COMPOSITE_CSV)
#   - GPCRdb_OPRM1_table.csv        (BW numbers, TM/ICL/ECL assignments)
#
# Provides to downstream scripts:
#   - composite / dms_long / dms_surface       (tidy DMS tables)
#   - count_lof_2sd_surf, count_gof_2sd_surf   (per-position disruption matrices)
#   - pos_to_sse                               (position -> SSE label)
#   - class_order, class_colors                (plotting aesthetics)
#   - out_dir / network_dir                    (output roots, auto-created)
#
# Cutoffs defined here: resi_min=66, resi_max=352. See README §4 for the rest.
# Run with working directory == this folder (network_tools/).
# ============================================================================

library(dplyr)
library(tidyr)
library(readr)

resi_min <- 66L
resi_max <- 352L
out_dir <- "dms_plots_output"
network_dir <- file.path(out_dir, "network_analysis")
if (!dir.exists(out_dir)) dir.create(out_dir)
if (!dir.exists(network_dir)) dir.create(network_dir)

# --- Ligand class definitions (from tuning-residue clustering; order = full -> strong -> intermediate -> weak -> antagonist -> no_ligand) ---
class_to_cols <- list(
  full         = c("Fentanyl_effect", "carfentanil_effect", "SR17018_effect", "DAMGO_effect"),
  strong       = c("Morphine_effect", "C6guano_effect"),
  intermediate = c("TRV130_effect", "PZM21_effect", "Methadone_effect"),
  weak         = c("Buprenorphine_effect", "Nalbuphine_effect", "MP_effect", "Butorphanol_effect"),
  antagonist   = c("Naloxone_effect", "Naltrexone_effect"),
  no_ligand    = c("FSK_effect")
)

all_ligand_cols <- unlist(class_to_cols, use.names = FALSE)

col_to_drug <- c(
  FSK_effect = "FSK",
  Naloxone_effect = "Naloxone", Naltrexone_effect = "Naltrexone",
  Buprenorphine_effect = "Buprenorphine", MP_effect = "MP",
  Butorphanol_effect = "Butorphanol", Nalbuphine_effect = "Nalbuphine",
  Morphine_effect = "Morphine", C6guano_effect = "C6guano",
  PZM21_effect = "PZM21", TRV130_effect = "TRV130",
  SR17018_effect = "SR17018", Fentanyl_effect = "Fentanyl",
  Methadone_effect = "Methadone", DAMGO_effect = "DAMGO",
  carfentanil_effect = "carfentanil"
)

col_to_class <- setNames(
  rep(names(class_to_cols), lengths(class_to_cols)),
  all_ligand_cols
)

class_colors <- c(
  full         = "#a6ea66",
  strong       = "#e63946",
  intermediate = "#ff9f0b",
  weak         = "#fde032",
  antagonist   = "#85fcff",
  no_ligand    = "#888888",
  shared       = "#2d508b"
)

# Plot/panel order: full, strong, intermediate, weak, antagonist, no_ligand
class_order <- c("full", "strong", "intermediate", "weak", "antagonist", "no_ligand")

# --- Load composite data ---
composite_csv <- Sys.getenv("MOR_COMPOSITE_CSV", unset = "composite_dms_scores_R.csv")
cat("Loading composite DMS data from ", composite_csv, "...\n", sep = "")
composite_raw <- read_csv(composite_csv, show_col_types = FALSE) %>%
  mutate(position = as.integer(position))
# Harmonize column name with class_to_cols (some exports use capital C)
if ("Carfentanil_effect" %in% names(composite_raw) && !"carfentanil_effect" %in% names(composite_raw)) {
  composite_raw <- composite_raw %>% rename(carfentanil_effect = Carfentanil_effect)
}

avail_cols <- intersect(all_ligand_cols, names(composite_raw))

# --- Synonymous variant statistics (computed on full unfiltered data) ---
synonymous <- composite_raw %>% filter(wildtype == mutation)

# Surface expression synonymous stats (data-driven cutoff)
surface_syn_mean <- mean(synonymous$Surface_effect, na.rm = TRUE)
surface_syn_sd   <- sd(synonymous$Surface_effect, na.rm = TRUE)
if (!is.finite(surface_syn_sd) || surface_syn_sd <= 0) surface_syn_sd <- 1e-10

surface_cutoff_1sd <- surface_syn_mean - surface_syn_sd
surface_cutoff_2sd <- surface_syn_mean - 2 * surface_syn_sd

cat("  Surface synonymous: mean =", round(surface_syn_mean, 4),
    ", SD =", round(surface_syn_sd, 4), "\n")
cat("  Surface cutoffs: 1SD =", round(surface_cutoff_1sd, 4),
    ", 2SD =", round(surface_cutoff_2sd, 4), "\n")

# `composite` = base dataset (all mutations with valid Surface_effect)
composite <- composite_raw %>% filter(!is.na(Surface_effect))

# `composite_surf` = surface-filtered: only mutations with WT-like
# surface expression (above syn_mean - 2SD), removing variants whose
# ligand effects are confounded by folding/trafficking disruption
composite_surf <- composite %>% filter(Surface_effect >= surface_cutoff_2sd)

cat("  Total mutations:", nrow(composite_raw), "\n")
cat("  After NA filter:", nrow(composite), "\n")
cat("  After surface 2SD filter:", nrow(composite_surf),
    sprintf("(removed %d surface-disrupted)\n",
            nrow(composite) - nrow(composite_surf)))
syn_stats <- tibble(
  col = avail_cols,
  syn_mean = sapply(avail_cols, function(c) mean(synonymous[[c]], na.rm = TRUE)),
  syn_sd   = sapply(avail_cols, function(c) sd(synonymous[[c]], na.rm = TRUE))
)
syn_stats$syn_sd[!is.finite(syn_stats$syn_sd) | syn_stats$syn_sd <= 0] <- 1e-10
syn_stats$cutoff_lof <- syn_stats$syn_mean + syn_stats$syn_sd
syn_stats$cutoff_gof <- syn_stats$syn_mean - syn_stats$syn_sd

# --- Position-level mean effect matrix ---
cat("Computing position-level mean effect matrix...\n")
pos_mean <- composite %>%
  filter(position >= resi_min, position <= resi_max) %>%
  group_by(position) %>%
  summarise(across(any_of(avail_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

pos_mean_mat <- as.matrix(pos_mean[, avail_cols])
rownames(pos_mean_mat) <- pos_mean$position

# --- Class-level mean effect (average across ligands within each class per position) ---
pos_class_mean <- pos_mean %>%
  mutate(
    no_ligand_mean    = rowMeans(across(any_of(class_to_cols$no_ligand)), na.rm = TRUE),
    antagonist_mean   = rowMeans(across(any_of(class_to_cols$antagonist)), na.rm = TRUE),
    weak_mean         = rowMeans(across(any_of(class_to_cols$weak)), na.rm = TRUE),
    intermediate_mean = rowMeans(across(any_of(class_to_cols$intermediate)), na.rm = TRUE),
    strong_mean       = rowMeans(across(any_of(class_to_cols$strong)), na.rm = TRUE),
    full_mean         = rowMeans(across(any_of(class_to_cols$full)), na.rm = TRUE)
  ) %>%
  select(position, ends_with("_mean"))

class_mean_cols <- paste0(class_order, "_mean")

# --- Selectivity: variance across the 5 class means per position ---
pos_selectivity <- pos_class_mean %>%
  rowwise() %>%
  mutate(
    selectivity_var = var(c_across(all_of(class_mean_cols))),
    selectivity_sd  = sd(c_across(all_of(class_mean_cols))),
    overall_mean    = mean(c_across(all_of(class_mean_cols))),
    dominant_class  = class_order[which.min(c_across(all_of(class_mean_cols)))]
  ) %>%
  ungroup()

# --- GPCR topology ---
cat("Loading GPCR topology...\n")
gpcrdb <- read_csv("GPCRdb_OPRM1_table.csv", show_col_types = FALSE) %>%
  filter(!is.na(pos), pos >= resi_min, pos <= resi_max) %>%
  mutate(pos = as.integer(pos))

sse_order <- c("N-tail", "TM1", "ICL1", "TM2", "ECL1", "TM3", "ICL2", "TM4",
               "ECL2", "TM5", "ICL3", "TM6", "ECL3", "TM7", "H8", "C-tail")

pos_to_sse <- gpcrdb %>% select(position = pos, SSE)

# GPCR profile layout (reused from connectome script)
sse_to_x <- c(
  "N-tail" = 0.5, TM1 = 1, ICL1 = 1.5, TM2 = 2, ECL1 = 2.5, TM3 = 3, ICL2 = 3.5,
  TM4 = 4, ECL2 = 4.5, TM5 = 5, ICL3 = 5.5, TM6 = 6, ECL3 = 6.5, TM7 = 7, H8 = 7.3, "C-tail" = 7.5
)
sse_to_y_range <- list(
  "N-tail" = c(1, 1), TM1 = c(0.9, 0.1), ICL1 = c(0.1, 0.1), TM2 = c(0.1, 0.9),
  ECL1 = c(0.9, 0.9), TM3 = c(0.9, 0.1), ICL2 = c(0.1, 0.1), TM4 = c(0.1, 0.9),
  ECL2 = c(0.9, 0.9), TM5 = c(0.9, 0.1), ICL3 = c(0.1, 0.1), TM6 = c(0.1, 0.9),
  ECL3 = c(0.9, 0.9), TM7 = c(0.9, 0.1), H8 = c(0.1, 0.1), "C-tail" = c(0.1, 0.1)
)

layout_df <- gpcrdb %>%
  mutate(SSE = factor(SSE, levels = sse_order)) %>%
  group_by(SSE) %>%
  mutate(n_in_sse = n(), idx = row_number(), t = (idx - 1) / pmax(n_in_sse - 1, 1)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    x_sse = if (SSE %in% names(sse_to_x)) sse_to_x[as.character(SSE)] else 4,
    y_range = list(if (as.character(SSE) %in% names(sse_to_y_range)) sse_to_y_range[[as.character(SSE)]] else c(0.5, 0.5)),
    y = y_range[[1]] + t * (y_range[[2]] - y_range[[1]])
  ) %>%
  ungroup() %>%
  select(position = pos, SSE, x = x_sse, y)

outline_segments <- tibble(
  seg = paste0("TM", 1:7), x = 1:7,
  xmin = (1:7) - 0.35, xmax = (1:7) + 0.35
)

# --- LoF/GoF position sets per ligand and per class ---
get_lof_positions <- function(col, min_var = 3L) {
  cutoff <- syn_stats$cutoff_lof[syn_stats$col == col]
  hits <- composite %>%
    filter(position >= resi_min, position <= resi_max, !is.na(.data[[col]]), .data[[col]] > cutoff) %>%
    group_by(position) %>% summarise(n = n(), .groups = "drop") %>%
    filter(n >= min_var) %>% pull(position)
  hits
}

get_gof_positions <- function(col, min_var = 3L) {
  cutoff <- syn_stats$cutoff_gof[syn_stats$col == col]
  hits <- composite %>%
    filter(position >= resi_min, position <= resi_max, !is.na(.data[[col]]), .data[[col]] < cutoff) %>%
    group_by(position) %>% summarise(n = n(), .groups = "drop") %>%
    filter(n >= min_var) %>% pull(position)
  hits
}

get_class_lof_positions <- function(cls, min_var = 3L) {
  cols <- class_to_cols[[cls]]
  pos_lists <- lapply(cols, get_lof_positions, min_var = min_var)
  Reduce(union, pos_lists)
}

get_class_gof_positions <- function(cls, min_var = 3L) {
  cols <- class_to_cols[[cls]]
  pos_lists <- lapply(cols, get_gof_positions, min_var = min_var)
  Reduce(union, pos_lists)
}

# =====================================================================
# DISRUPTION COUNT MATRICES (v2 core metric)
# For each (position, ligand): count how many mutations exceed threshold
# =====================================================================
cat("Computing disruption count matrices (LoF/GoF x 1SD/2SD)...\n")

syn_stats$cutoff_lof_2sd <- syn_stats$syn_mean + 2 * syn_stats$syn_sd
syn_stats$cutoff_gof_2sd <- syn_stats$syn_mean - 2 * syn_stats$syn_sd

compute_disruption_counts <- function(direction, n_sd) {
  positions <- sort(unique(composite$position[composite$position >= resi_min & composite$position <= resi_max]))
  mat <- matrix(0L, nrow = length(positions), ncol = length(avail_cols),
                dimnames = list(positions, avail_cols))
  for (col in avail_cols) {
    if (direction == "lof") {
      cutoff <- syn_stats$syn_mean[syn_stats$col == col] + n_sd * syn_stats$syn_sd[syn_stats$col == col]
      hits <- composite %>%
        filter(position >= resi_min, position <= resi_max, !is.na(.data[[col]]), .data[[col]] > cutoff)
    } else {
      cutoff <- syn_stats$syn_mean[syn_stats$col == col] - n_sd * syn_stats$syn_sd[syn_stats$col == col]
      hits <- composite %>%
        filter(position >= resi_min, position <= resi_max, !is.na(.data[[col]]), .data[[col]] < cutoff)
    }
    counts <- hits %>% group_by(position) %>% summarise(n = n(), .groups = "drop")
    mat[as.character(counts$position), col] <- counts$n
  }
  mat
}

count_lof_1sd <- compute_disruption_counts("lof", 1)
count_lof_2sd <- compute_disruption_counts("lof", 2)
count_gof_1sd <- compute_disruption_counts("gof", 1)
count_gof_2sd <- compute_disruption_counts("gof", 2)

# Total mutations tested per position (for normalization context)
n_mutations_per_pos <- composite %>%
  filter(position >= resi_min, position <= resi_max) %>%
  group_by(position) %>% summarise(n_total = n(), .groups = "drop")

# --- Class-level disruption counts (mean count across ligands within class) ---
compute_class_counts <- function(count_mat) {
  positions <- as.integer(rownames(count_mat))
  df <- as.data.frame(count_mat) %>% mutate(position = positions)
  for (cls in class_order) {
    cls_cols <- intersect(class_to_cols[[cls]], colnames(count_mat))
    df[[paste0(cls, "_count")]] <- rowMeans(df[, cls_cols, drop = FALSE], na.rm = TRUE)
  }
  df %>% select(position, all_of(paste0(class_order, "_count")))
}

class_count_lof_1sd <- compute_class_counts(count_lof_1sd)
class_count_lof_2sd <- compute_class_counts(count_lof_2sd)
class_count_gof_1sd <- compute_class_counts(count_gof_1sd)
class_count_gof_2sd <- compute_class_counts(count_gof_2sd)

class_count_cols <- paste0(class_order, "_count")

# --- Count-based selectivity (variance of class-level disruption count) ---
compute_count_selectivity <- function(class_count_df, label) {
  class_count_df %>%
    rowwise() %>%
    mutate(
      selectivity_var = var(c_across(all_of(class_count_cols))),
      selectivity_sd  = sd(c_across(all_of(class_count_cols))),
      overall_count   = mean(c_across(all_of(class_count_cols))),
      dominant_class  = class_order[which.max(c_across(all_of(class_count_cols)))]
    ) %>%
    ungroup()
}

sel_lof_1sd <- compute_count_selectivity(class_count_lof_1sd, "lof_1sd")
sel_lof_2sd <- compute_count_selectivity(class_count_lof_2sd, "lof_2sd")
sel_gof_1sd <- compute_count_selectivity(class_count_gof_1sd, "gof_1sd")
sel_gof_2sd <- compute_count_selectivity(class_count_gof_2sd, "gof_2sd")

# =====================================================================
# SURFACE-FILTERED DISRUPTION COUNTS
# Only count mutations with WT-like surface expression (above
# syn_mean - 2SD), removing variants whose ligand effects are
# confounded by folding/trafficking disruption.
# =====================================================================
cat("Computing surface-filtered disruption count matrices...\n")

compute_disruption_counts_filtered <- function(data_subset, direction, n_sd) {
  positions <- sort(unique(data_subset$position[data_subset$position >= resi_min & data_subset$position <= resi_max]))
  mat <- matrix(0L, nrow = length(positions), ncol = length(avail_cols),
                dimnames = list(positions, avail_cols))
  for (col in avail_cols) {
    if (direction == "lof") {
      cutoff <- syn_stats$syn_mean[syn_stats$col == col] + n_sd * syn_stats$syn_sd[syn_stats$col == col]
      hits <- data_subset %>%
        filter(position >= resi_min, position <= resi_max, !is.na(.data[[col]]), .data[[col]] > cutoff)
    } else {
      cutoff <- syn_stats$syn_mean[syn_stats$col == col] - n_sd * syn_stats$syn_sd[syn_stats$col == col]
      hits <- data_subset %>%
        filter(position >= resi_min, position <= resi_max, !is.na(.data[[col]]), .data[[col]] < cutoff)
    }
    counts <- hits %>% group_by(position) %>% summarise(n = n(), .groups = "drop")
    matched_pos <- as.character(counts$position[as.character(counts$position) %in% rownames(mat)])
    mat[matched_pos, col] <- counts$n[as.character(counts$position) %in% rownames(mat)]
  }
  mat
}

count_lof_1sd_surf <- compute_disruption_counts_filtered(composite_surf, "lof", 1)
count_lof_2sd_surf <- compute_disruption_counts_filtered(composite_surf, "lof", 2)
count_gof_1sd_surf <- compute_disruption_counts_filtered(composite_surf, "gof", 1)
count_gof_2sd_surf <- compute_disruption_counts_filtered(composite_surf, "gof", 2)

class_count_lof_1sd_surf <- compute_class_counts(count_lof_1sd_surf)
class_count_lof_2sd_surf <- compute_class_counts(count_lof_2sd_surf)
class_count_gof_1sd_surf <- compute_class_counts(count_gof_1sd_surf)
class_count_gof_2sd_surf <- compute_class_counts(count_gof_2sd_surf)

sel_lof_1sd_surf <- compute_count_selectivity(class_count_lof_1sd_surf, "lof_1sd_surf")
sel_lof_2sd_surf <- compute_count_selectivity(class_count_lof_2sd_surf, "lof_2sd_surf")
sel_gof_1sd_surf <- compute_count_selectivity(class_count_gof_1sd_surf, "gof_1sd_surf")
sel_gof_2sd_surf <- compute_count_selectivity(class_count_gof_2sd_surf, "gof_2sd_surf")

cat("  Surface-filtered matrices: ", nrow(count_lof_1sd_surf), " positions x ",
    ncol(count_lof_1sd_surf), " ligands\n", sep = "")

# --- Surface-filtered LoF/GoF position getters ---
get_lof_positions_surf <- function(col, min_var = 3L) {
  cutoff <- syn_stats$cutoff_lof[syn_stats$col == col]
  hits <- composite_surf %>%
    filter(position >= resi_min, position <= resi_max, !is.na(.data[[col]]), .data[[col]] > cutoff) %>%
    group_by(position) %>% summarise(n = n(), .groups = "drop") %>%
    filter(n >= min_var) %>% pull(position)
  hits
}

get_gof_positions_surf <- function(col, min_var = 3L) {
  cutoff <- syn_stats$cutoff_gof[syn_stats$col == col]
  hits <- composite_surf %>%
    filter(position >= resi_min, position <= resi_max, !is.na(.data[[col]]), .data[[col]] < cutoff) %>%
    group_by(position) %>% summarise(n = n(), .groups = "drop") %>%
    filter(n >= min_var) %>% pull(position)
  hits
}

get_class_lof_positions_surf <- function(cls, min_var = 3L) {
  cols <- class_to_cols[[cls]]
  Reduce(union, lapply(cols, get_lof_positions_surf, min_var = min_var))
}

get_class_gof_positions_surf <- function(cls, min_var = 3L) {
  cols <- class_to_cols[[cls]]
  Reduce(union, lapply(cols, get_gof_positions_surf, min_var = min_var))
}

cat("  Disruption count matrices: ", nrow(count_lof_1sd), " positions x ", ncol(count_lof_1sd), " ligands\n", sep = "")

# --- Curated positions of biological interest ---
# 100: partial agonists far more sensitive than full agonists
# 119: validated -- mutations convert antagonists/partial agonists to full agonists
# 149: known functional site
curated_positions_base <- c(100L, 119L, 149L)

cat("Shared data preparation complete.\n")
cat("  Positions:", nrow(pos_mean), "\n")
cat("  Ligands:", length(avail_cols), "\n")
cat("  Classes:", paste(names(class_to_cols), collapse = ", "), "\n")

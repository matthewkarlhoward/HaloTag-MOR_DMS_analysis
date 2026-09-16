## howard_theme.R — Shared ggplot2 theme + helpers for all plots in this repo.
##
## Source at the top of every plot script:
##   source("/Users/mkh/GitHub/mor_dms_analysis/howard_theme.R")
## then add `+ base_theme` to your ggplot() calls. Per-script tweaks layer on
## top via `+ theme(axis.text.x = element_text(angle = 45, hjust = 1))` etc.
##
## Provides:
##   lw          line width (0.5pt converted to mm — matches existing plots)
##   base_theme  the canonical ggplot theme used across all panels
##   PT_SIZE     default point size for scatters
##   PAL_*       color palettes (diverging / expression / discrete)
##   FIG_*       figure widths (mm) for Cell-family journals
##   save_panel()       publication PDF export helper
##   new_run_dir()      timestamped iterative output dir
##   summarize_missing(), write_run_summary()
##
## Style: 6pt text, 0.5pt axis lines, white background, no grid, classic axis.
## Matches the inline theme used by qc_report.R and most existing plots, so
## figure outputs are unchanged.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
})

# ── Core theme ───────────────────────────────────────────────────────────────

# 0.5pt in ggplot linewidth units (1pt = 0.353 mm; 0.5pt = 0.176 mm)
lw      <- 0.5 / 2.835
PT_SIZE <- 0.8

base_theme <- theme_classic(base_size = 6) +
  theme(
    text             = element_text(color = "black", size = 6),
    axis.text        = element_text(color = "black", size = 6),
    axis.title       = element_text(color = "black", size = 6),
    axis.line        = element_line(color = "black", linewidth = lw),
    axis.ticks       = element_line(color = "black", linewidth = lw),
    strip.background = element_blank(),
    strip.text       = element_text(size = 6),
    legend.text      = element_text(size = 6),
    legend.title     = element_text(size = 6),
    legend.key.size  = unit(0.25, "cm"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Set as default so plots without an explicit `+ base_theme` still inherit it.
theme_set(base_theme)

# ── Color palettes ───────────────────────────────────────────────────────────

# Diverging: for sign-flip / pharmacology comparisons
PAL_DIVERGING <- c(
  "#2166AC",  # blue  — negative / morphine-like
  "#F7F7F7",  # white — neutral
  "#D6604D"   # red   — positive / fentanyl-like
)

# Sequential: for expression / surface scores
PAL_EXPRESSION <- c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#084594")

# Discrete: ligand classes or receptor systems
PAL_DISCRETE <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
  "#FF7F00", "#A65628", "#F781BF", "#999999"
)

# ── Figure dimensions (Cell-family journals, mm) ─────────────────────────────

FIG_SINGLE      <- 85    # single column
FIG_ONEHALF     <- 114   # 1.5 column
FIG_DOUBLE      <- 180   # full width
FIG_HEIGHT_UNIT <- 50    # typical panel height

# ── Export helper ────────────────────────────────────────────────────────────

# Save a PDF at correct dimensions for publication.
# width and height in mm (Cell single column = 85mm, double = 180mm)
save_panel <- function(plot, filename, width_mm, height_mm,
                       output_dir = "figures/main/") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  ggsave(
    filename = file.path(output_dir, filename),
    plot     = plot,
    width    = width_mm,
    height   = height_mm,
    units    = "mm",
    device   = cairo_pdf
  )
  message("Saved: ", file.path(output_dir, filename))
}

# ── Run-tracking helpers ─────────────────────────────────────────────────────

# Create a timestamped output directory for iterative plot runs
new_run_dir <- function(root = "plots/iterative") {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_dir <- file.path(root, ts)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir
}

# Return a tibble of NA counts and fractions for selected columns in df
summarize_missing <- function(df, cols) {
  cols <- cols[cols %in% names(df)]
  if (length(cols) == 0) {
    return(data.frame(column = character(), missing_n = integer(),
                      missing_frac = numeric()))
  }
  tibble::tibble(column = cols) %>%
    rowwise() %>%
    mutate(missing_n    = sum(is.na(df[[column]])),
           missing_frac = missing_n / nrow(df)) %>%
    ungroup()
}

# Write a plain-text run summary (timestamp, row/col counts, missingness, notes)
write_run_summary <- function(df, out_dir, used_cols, notes = NULL) {
  miss_tbl <- summarize_missing(df, used_cols)
  summary_path <- file.path(out_dir, "run_summary.txt")

  lines <- c(
    paste0("timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("rows: ", nrow(df)),
    paste0("columns: ", ncol(df)),
    "",
    "columns_used:",
    paste0("- ", used_cols)
  )

  if (nrow(miss_tbl) > 0) {
    lines <- c(lines, "", "missingness:")
    miss_lines <- sprintf("- %s: %d (%.2f%%)",
                          miss_tbl$column, miss_tbl$missing_n,
                          miss_tbl$missing_frac * 100)
    lines <- c(lines, miss_lines)
  }

  if (!is.null(notes) && length(notes) > 0) {
    lines <- c(lines, "", "notes:", paste0("- ", notes))
  }

  writeLines(lines, summary_path)
  message("Saved: ", summary_path)
  summary_path
}

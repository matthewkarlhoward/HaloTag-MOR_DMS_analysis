## fig_template.R — starting point for a new plot script.
## Copy into plots/<subdir>/<name>.R and fill in.

# ── Setup ────────────────────────────────────────────────────────────────────
source("/Users/mkh/GitHub/mor_dms_analysis/howard_theme.R")
library(tidyverse)

# ── Load data ────────────────────────────────────────────────────────────────
dms <- read_csv(
  "/Users/mkh/GitHub/mor_dms_analysis/dms_scores/composite_dms_scores.csv",
  show_col_types = FALSE
)

# ── Filter / transform ───────────────────────────────────────────────────────


# ── Plot ─────────────────────────────────────────────────────────────────────
p <- ggplot(dms, aes(x = , y = )) +
  geom_point(size = PT_SIZE, stroke = 0) +
  labs(x = "", y = "", title = "") +
  base_theme

# ── Export ───────────────────────────────────────────────────────────────────
# save_panel() writes a publication PDF (Cell-style mm dimensions, cairo_pdf)
save_panel(p, "figNx.pdf",
           width_mm  = FIG_SINGLE,
           height_mm = FIG_HEIGHT_UNIT,
           output_dir = "plots/<subdir>")

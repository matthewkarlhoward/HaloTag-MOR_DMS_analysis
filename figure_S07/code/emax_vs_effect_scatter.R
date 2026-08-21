library(tidyverse)
library(patchwork)

setwd("/Users/mkh/GitHub/mor_dms_analysis/dms_scores")
dms <- read_csv("composite_dms_scores.csv", show_col_types = FALSE) %>%
  mutate(
    Fentanyl_ec50_logM = Fentanyl_ec50 - 13,
    Morphine_ec50_logM = Morphine_ec50 - 12.5,
    DAMGO_ec50_logM    = DAMGO_ec50    - 10
  )

missense <- dms %>% filter(type == "missense")

# ── Per-position means ────────────────────────────────────────────────────────
drugs <- list(
  list(
    name        = "Morphine",
    effect_col  = "Morphine_effect",
    emax_col    = "Morphine_emax",
    curve_col   = "Morphine_curve_type"
  ),
  list(
    name        = "Fentanyl",
    effect_col  = "Fentanyl_effect",
    emax_col    = "Fentanyl_emax",
    curve_col   = "Fentanyl_curve_type"
  ),
  list(
    name        = "DAMGO",
    effect_col  = "DAMGO_effect",
    emax_col    = "DAMGO_emax",
    curve_col   = "DAMGO_curve_type"
  )
)

# Synonymous references
syn <- dms %>% filter(type == "synonymous")

make_pos_df <- function(drug) {
  sig_syn    <- syn %>% filter(.data[[drug$curve_col]] == "sigmoid")
  syn_effect <- mean(sig_syn %>% pull(all_of(drug$effect_col)), na.rm = TRUE)
  syn_emax   <- mean(sig_syn %>% pull(all_of(drug$emax_col)),   na.rm = TRUE)

  # Both metrics averaged over sigmoid missense variants only
  sig_missense <- missense %>% filter(.data[[drug$curve_col]] == "sigmoid")

  sig_missense %>%
    group_by(position) %>%
    summarise(
      mean_effect = mean(.data[[drug$effect_col]], na.rm = TRUE) - syn_effect,
      mean_emax   = mean(.data[[drug$emax_col]],   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(drug = drug$name, syn_emax = syn_emax)
}

plot_data <- map_dfr(drugs, make_pos_df)

# ── Base theme ────────────────────────────────────────────────────────────────
lw <- 0.5 / 2.835

base_theme <- theme_classic(base_size = 6) +
  theme(
    text             = element_text(color = "black", size = 6, face = "plain"),
    axis.text        = element_text(color = "black", size = 6, face = "plain"),
    axis.title       = element_text(color = "black", size = 6, face = "plain"),
    axis.line        = element_line(color = "black", linewidth = lw),
    axis.ticks       = element_line(color = "black", linewidth = lw),
    plot.title       = element_text(color = "black", size = 6, face = "plain", hjust = 0.5),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# ── Plots ─────────────────────────────────────────────────────────────────────
x_lab <- expression("Position mean signaling score from 10 " * mu * "M ligand screen")
y_lab <- expression(atop("Position mean E"[max] * " from",
                         "concentration-response experiment"))

plots <- map(drugs, function(drug) {
  drug_name <- drug$name

  df <- plot_data %>% filter(drug == drug_name)

  r2 <- cor(df$mean_effect, df$mean_emax, use = "complete.obs")^2
  r2_label <- paste0("R\u00b2 = ", round(r2, 3))

  ggplot(df, aes(x = mean_effect, y = mean_emax)) +
    geom_hline(yintercept = df$syn_emax[1], color = "black",
               linewidth = lw, linetype = "dashed") +
    geom_vline(xintercept = 0, color = "black",
               linewidth = lw, linetype = "dashed") +
    geom_point(color = "black", size = 0.6, alpha = 0.8) +
    labs(title    = drug_name,
         subtitle = r2_label,
         x        = x_lab,
         y        = y_lab) +
    base_theme +
    theme(plot.subtitle = element_text(color = "black", size = 6, face = "plain",
                                       hjust = 0.5, margin = margin(t = 0, b = 2)))
})

combined <- (plots[[1]] | plots[[2]] | plots[[3]]) +
  plot_layout(axis_titles = "collect")

out_stem <- "/Users/mkh/GitHub/mor_dms_analysis/plots/scatter/ec50_emax/emax_vs_effect_scatter"
ggsave(paste0(out_stem, ".pdf"), combined, width = 165, height = 50, units = "mm")
ggsave(paste0(out_stem, ".png"), combined, width = 165, height = 50, units = "mm",
       dpi = 600, bg = "white")
message("Saved: emax_vs_effect_scatter.pdf / .png")

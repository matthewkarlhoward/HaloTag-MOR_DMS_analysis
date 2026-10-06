library(tidyverse)
library(readxl)
library(patchwork)

setwd("/Users/mkh/GitHub/mor_dms_analysis/trupath_vs_camp_efficacy")

# ── Efficacy data ──────────────────────────────────────────────────────────────
# Source spreadsheet headers were flipped:
#   TRUPATH_Gi1_Emax_pct_DAMGO is already normalised to DAMGO
#   cAMP_response_raw is the raw cAMP BRET deltaEmax and needs normalising
eff <- read_excel("WT_cAMP_vs_TRUPATH_All_Drugs.xlsx",
                  sheet = "cAMP vs TRUPATH (WT)") %>%
  rename(trupath_pct = `TRUPATH_Gi1_Emax_pct_DAMGO`,
         camp_raw    = `cAMP_response_raw`) %>%
  filter(!is.na(trupath_pct), !is.na(camp_raw))

damgo_camp <- eff$camp_raw[eff$Ligand == "DAMGO"]
eff <- eff %>% mutate(camp_pct = camp_raw / damgo_camp * 100)

# ── EC50 data ──────────────────────────────────────────────────────────────────
ec50 <- read_excel("cAMP_trupath_ec50_copmp.xlsx", sheet = "Sheet1") %>%
  rename(Ligand       = ligand,
         trupath_ec50 = `TRUPATH EC50`,
         camp_ec50    = `cAMP EC50`) %>%
  filter(!is.na(trupath_ec50), !is.na(camp_ec50))

# ── Shared styling ─────────────────────────────────────────────────────────────
lw <- 0.5 / 2.835

txt <- element_text(color = "black", size = 6, face = "plain", family = "Helvetica")

base_theme <- theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(
    text              = txt,
    axis.text         = txt,
    axis.text.x       = txt,
    axis.text.y       = txt,
    axis.title        = txt,
    axis.title.x      = txt,
    axis.title.y      = txt,
    axis.line         = element_line(color = "black", linewidth = lw),
    axis.ticks        = element_line(color = "black", linewidth = lw),
    plot.title        = element_text(color = "black", size = 6, face = "plain", family = "Helvetica", hjust = 0.5),
    plot.background   = element_rect(fill = NA, color = NA),
    panel.background  = element_rect(fill = NA, color = NA),
    plot.margin       = margin(1, 1, 1, 1),
    strip.text        = txt,
    legend.position   = "bottom",
    legend.title      = element_blank(),
    legend.text       = txt,
    legend.key.size   = unit(2, "mm"),
    legend.key        = element_rect(fill = NA, color = NA),
    legend.background = element_rect(fill = NA, color = NA),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1, "mm"),
    legend.spacing.x  = unit(1, "mm"),
    legend.spacing.y  = unit(0.4, "mm"),
    legend.box.margin = margin(0, 0, 0, 0),
    aspect.ratio      = 1
  )

ligand_palette <- c(
  "DAMGO"         = "#1F77B4",
  "Fentanyl"      = "#FF7F0E",
  "Morphine"      = "#2CA02C",
  "Naloxone"      = "#D62728",
  "Buprenorphine" = "#9467BD",
  "Methadone"     = "#8C564B",
  "PZM21"         = "#E377C2",
  "Carfentanil"   = "#17BECF",
  "Nalbuphine"    = "#BCBD22",
  "Butorphanol"   = "#7F7F7F",
  "TRV130"        = "#AEC7E8",
  "MP"            = "#FFBB78"
)

# Legend order driven by efficacy panel: descending TRUPATH % DAMGO
legend_order <- eff %>% arrange(desc(trupath_pct)) %>% pull(Ligand)

# ── Efficacy panel ─────────────────────────────────────────────────────────────
rho_eff <- cor(eff$camp_pct, eff$trupath_pct, method = "spearman")
lab_eff <- sprintf("rho == %.2f", rho_eff)

p_eff <- ggplot(eff, aes(x = trupath_pct, y = camp_pct, color = Ligand)) +
  geom_abline(slope = 1, intercept = 0,
              color = "grey70", linewidth = lw, linetype = "dashed") +
  geom_point(shape = 16, size = 1.0, alpha = 0.95) +
  scale_color_manual(values = ligand_palette, breaks = legend_order) +
  annotate("text", x = Inf, y = -Inf, label = lab_eff, parse = TRUE,
           family = "Helvetica",
           hjust = 1.1, vjust = -0.8, size = 6 / .pt, color = "black") +
  scale_x_continuous(limits = c(0, 105), expand = c(0, 0),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) ifelse(x %in% c(25, 75), "", as.character(x))) +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) ifelse(x %in% c(25, 75), "", as.character(x))) +
  coord_cartesian(clip = "off") +
  labs(title = expression(E[max]*", % DAMGO"), x = NULL, y = "Pooled cAMP") +
  base_theme +
  theme(axis.title.x = element_blank())

# ── EC50 panel ─────────────────────────────────────────────────────────────────
rho_ec50 <- cor(ec50$camp_ec50, ec50$trupath_ec50, method = "spearman")
lab_ec50 <- sprintf("rho == %.2f", rho_ec50)

p_ec50 <- ggplot(ec50, aes(x = trupath_ec50, y = camp_ec50, color = Ligand)) +
  geom_abline(slope = 1, intercept = 0,
              color = "grey70", linewidth = lw, linetype = "dashed") +
  geom_point(shape = 16, size = 1.0, alpha = 0.95) +
  scale_color_manual(values = ligand_palette, breaks = legend_order) +
  annotate("text", x = Inf, y = -Inf, label = lab_ec50, parse = TRUE,
           family = "Helvetica",
           hjust = 1.1, vjust = -0.8, size = 6 / .pt, color = "black") +
  scale_x_continuous(limits = c(-13, -5), expand = c(0, 0),
                     breaks = seq(-12, -6, by = 1),
                     labels = function(x) ifelse(x %% 2 == 0, as.character(x), "")) +
  scale_y_continuous(limits = c(-13, -5), expand = c(0, 0),
                     breaks = seq(-12, -6, by = 1),
                     labels = function(x) ifelse(x %% 2 == 0, as.character(x), "")) +
  coord_cartesian(clip = "off") +
  labs(title = expression(EC[50]*", log[ligand]"), x = NULL, y = NULL) +
  base_theme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

# ── File 1: two plots side-by-side, NO legend, shared x-axis caption ─────────
plots_only <- (p_eff + p_ec50) +
  plot_annotation(
    caption = "TRUPATH Gi1",
    theme = theme(
      plot.caption = element_text(color = "black", size = 6, family = "Helvetica",
                                  face = "plain", hjust = 0.5,
                                  margin = margin(t = 1, b = 0))
    )
  ) &
  theme(legend.position = "none")

ggsave("camp_vs_trupath_plots.pdf", plots_only,
       width = 63, height = 37, units = "mm",
       bg = "transparent")
message("Saved: ", file.path(getwd(), "camp_vs_trupath_plots.pdf"))

# ── File 2: just the legend ──────────────────────────────────────────────────
# Build a throw-away plot styled identically, then extract its guide-box grob.
p_for_legend <- p_eff +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 2, byrow = FALSE,
                              keyheight = unit(2, "mm"),
                              keywidth  = unit(2, "mm")))

g <- ggplotGrob(p_for_legend)
legend_idx <- grep("guide-box", sapply(g$grobs, function(x) x$name))
# Prefer a non-zero-sized box (newer ggplot stores empty placeholders too)
legend_grob <- NULL
for (i in legend_idx) {
  cand <- g$grobs[[i]]
  if (!inherits(cand, "zeroGrob")) { legend_grob <- cand; break }
}
if (is.null(legend_grob)) stop("Could not locate legend grob")

# Size the page to the legend's natural footprint (+ a 1 mm margin on each side).
leg_w_mm <- grid::convertWidth(sum(legend_grob$widths),  "mm", valueOnly = TRUE)
leg_h_mm <- grid::convertHeight(sum(legend_grob$heights), "mm", valueOnly = TRUE)
pad_mm   <- 1
grDevices::pdf("camp_vs_trupath_legend.pdf",
               width  = (leg_w_mm + 2 * pad_mm) / 25.4,
               height = (leg_h_mm + 2 * pad_mm) / 25.4,
               bg = "transparent")
grid::grid.newpage()
grid::pushViewport(grid::viewport(
  x = grid::unit(pad_mm, "mm"), y = grid::unit(pad_mm, "mm"),
  width  = grid::unit(leg_w_mm, "mm"),
  height = grid::unit(leg_h_mm, "mm"),
  just = c("left", "bottom")
))
grid::grid.draw(legend_grob)
grid::popViewport()
dev.off()
message(sprintf("Saved: %s  (%.1f x %.1f mm)",
                file.path(getwd(), "camp_vs_trupath_legend.pdf"),
                leg_w_mm + 2 * pad_mm, leg_h_mm + 2 * pad_mm))

library(tidyverse)
library(readxl)

setwd("/Users/mkh/GitHub/mor_dms_analysis/trupath_vs_camp_efficacy")

dat <- read_excel("cAMP_trupath_ec50_copmp.xlsx", sheet = "Sheet1") %>%
  rename(Ligand      = ligand,
         trupath_ec50 = `TRUPATH EC50`,
         camp_ec50    = `cAMP EC50`) %>%
  filter(!is.na(trupath_ec50), !is.na(camp_ec50))

# 0.5 pt linewidth in ggplot2 units
lw <- 0.5 / 2.835

base_theme <- theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(
    text             = element_text(color = "black", size = 6, face = "plain", family = "Helvetica"),
    axis.text        = element_text(color = "black", size = 6, face = "plain", family = "Helvetica"),
    axis.title       = element_text(color = "black", size = 6, face = "plain", family = "Helvetica"),
    axis.line        = element_line(color = "black", linewidth = lw),
    axis.ticks       = element_line(color = "black", linewidth = lw),
    plot.title       = element_text(color = "black", size = 6, face = "plain", family = "Helvetica", hjust = 0.5),
    plot.background  = element_rect(fill = NA, color = NA),
    panel.background = element_rect(fill = NA, color = NA),
    legend.position  = "right",
    legend.title     = element_blank(),
    legend.text      = element_text(color = "black", size = 6, family = "Helvetica"),
    legend.key.size  = unit(2, "mm"),
    legend.key       = element_rect(fill = NA, color = NA),
    legend.background = element_rect(fill = NA, color = NA),
    legend.margin    = margin(0, 0, 0, 0),
    legend.box.spacing = unit(4, "mm"),
    aspect.ratio     = 1
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

# Spearman correlation across all WT ligands
rho_val <- cor(dat$camp_ec50, dat$trupath_ec50, method = "spearman")
corr_label <- sprintf("rho == %.2f", rho_val)

# Legend ordered by TRUPATH potency (most potent = most negative log[M] first)
legend_order <- dat %>% arrange(trupath_ec50) %>% pull(Ligand)

p <- ggplot(dat, aes(x = camp_ec50, y = trupath_ec50, color = Ligand)) +
  geom_abline(slope = 1, intercept = 0,
              color = "grey70", linewidth = lw, linetype = "dashed") +
  geom_point(shape = 16, size = 1.0, alpha = 0.95) +
  scale_color_manual(values = ligand_palette,
                     breaks = legend_order) +
  annotate("text", x = Inf, y = -Inf, label = corr_label, parse = TRUE,
           family = "Helvetica",
           hjust = 1.1, vjust = -0.8, size = 6 / .pt, color = "black") +
  scale_x_continuous(limits = c(-13, -5), expand = c(0, 0),
                     breaks = seq(-12, -6, by = 1),
                     labels = function(x) ifelse(x %% 2 == 0, as.character(x), "")) +
  scale_y_continuous(limits = c(-13, -5), expand = c(0, 0),
                     breaks = seq(-12, -6, by = 1),
                     labels = function(x) ifelse(x %% 2 == 0, as.character(x), "")) +
  coord_cartesian(clip = "off") +
  labs(
    x = "cAMP EC50 (log[M])",
    y = "TRUPATH Gi1 EC50 (log[M])"
  ) +
  base_theme

ggsave("camp_vs_trupath_ec50_scatter.pdf", p,
       width = 68, height = 40, units = "mm",
       bg = "transparent")
message("Saved: ", file.path(getwd(), "camp_vs_trupath_ec50_scatter.pdf"))

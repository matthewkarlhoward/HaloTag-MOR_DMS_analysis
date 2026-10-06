library(tidyverse)
library(readxl)

setwd("/Users/mkh/GitHub/mor_dms_analysis/trupath_vs_camp_efficacy")

# Headers in the source spreadsheet are flipped:
#   TRUPATH_Gi1_Emax_pct_DAMGO is already normalised to DAMGO.
#   cAMP_response_raw is the raw cAMP BRET deltaEmax and needs normalising.
dat <- read_excel("WT_cAMP_vs_TRUPATH_All_Drugs.xlsx",
                  sheet = "cAMP vs TRUPATH (WT)") %>%
  rename(trupath_pct = `TRUPATH_Gi1_Emax_pct_DAMGO`,
         camp_raw    = `cAMP_response_raw`) %>%
  filter(!is.na(trupath_pct), !is.na(camp_raw))

# Normalise cAMP BRET ΔEmax to % DAMGO (DAMGO = 100)
damgo_camp <- dat$camp_raw[dat$Ligand == "DAMGO"]
dat <- dat %>%
  mutate(camp_pct = camp_raw / damgo_camp * 100)

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

# Pearson correlation across all WT ligands
rho_val <- cor(dat$camp_pct, dat$trupath_pct, method = "spearman")
corr_label <- sprintf("rho == %.2f", rho_val)

legend_order <- dat %>% arrange(desc(trupath_pct)) %>% pull(Ligand)

p <- ggplot(dat, aes(x = camp_pct, y = trupath_pct, color = Ligand)) +
  geom_abline(slope = 1, intercept = 0,
              color = "grey70", linewidth = lw, linetype = "dashed") +
  geom_point(shape = 16, size = 1.0, alpha = 0.95) +
  scale_color_manual(values = ligand_palette,
                     breaks = legend_order) +
  annotate("text", x = Inf, y = -Inf, label = corr_label, parse = TRUE,
           family = "Helvetica",
           hjust = 1.1, vjust = -0.8, size = 6 / .pt, color = "black") +
  scale_x_continuous(limits = c(0, 105), expand = c(0, 0),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) ifelse(x %in% c(25, 75), "", as.character(x))) +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) ifelse(x %in% c(25, 75), "", as.character(x))) +
  coord_cartesian(clip = "off") +
  labs(
    x = "cAMP response (% DAMGO)",
    y = "TRUPATH Gi1 (% DAMGO)"
  ) +
  base_theme

ggsave("camp_vs_trupath_scatter.pdf", p,
       width = 68, height = 40, units = "mm",
       bg = "transparent")
message("Saved: ", file.path(getwd(), "camp_vs_trupath_scatter.pdf"))

## Orthosteric (ligand-binding pocket) mean-effect heatmaps.
##
## Mirrors the style of the other interface heatmaps:
##   x = receptor position
##   y = ligand
##   fill = mean missense effect for that ligand at that position
##   top axis  = GPCRdb numbering
##   bottom axis = MOR residue (single-letter + position)
##
## Two outputs:
##   orthosteric_contact_map_first_shell.pdf
##   orthosteric_contact_map_second_shell.pdf
##
## Shell membership is read directly from the per-structure ligand-distance
## CSVs (experimental + Chai-1 predicted poses for compounds with no co-crystal),
## plus a separate buprenorphine.pdb analysis (no CIF available). All distance
## CSVs already encode positions in human-MOR numbering via residue_number_human
## (mouse structures get +2 offset based on D147→D149 detection in
## analyze_all_structures.py / analyze_buprenorphine_distances.py).
##
## Effect values come from composite_dms_scores.csv — the same dataset used
## by the gprotein / arrestin / interface plots — keyed on `position` (also
## human-MOR numbering).

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

aa3to1 <- c(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E",
            GLY="G", HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F",
            PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")

OUT_DIR <- "/Users/mkh/GitHub/mor_dms_analysis/plots/interfaces/orthosteric"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ── Load shell membership (per-ligand, human numbering) ─────────────────────
DIST_DIR <- "/Users/mkh/GitHub/mor_dms_analysis/structures/processed/ligand_distances"

drug_rename <- c(
  "Oliceridine"               = "TRV130",
  "Lofentanil"                = "Carfentanil",
  "Mitragynine_Pseudoindoxyl" = "MP",
  "C6-guano"                  = "C6guano"
)

dist_main <- read_csv(file.path(DIST_DIR, "experimental_and_chai_combined.csv"),
                      show_col_types = FALSE)
dist_bup  <- read_csv(file.path(DIST_DIR, "per_pdb_legacy/buprenorphine_distances.csv"),
                      show_col_types = FALSE) %>%
  mutate(drug_name = "Buprenorphine")

dist_all <- bind_rows(dist_main, dist_bup) %>%
  mutate(drug_name = recode(drug_name, !!!drug_rename))

## Per (drug, position) take the *closest* shell label across all structures
## of that ligand: first_shell > second_shell > none.
shell_rank <- c(first_shell = 1L, second_shell = 2L, none = 3L)

shell_long <- dist_all %>%
  filter(!is.na(shell)) %>%
  mutate(shell_rk = shell_rank[shell]) %>%
  group_by(drug = drug_name, position = residue_number_human) %>%
  summarise(shell = names(shell_rank)[min(shell_rk)], .groups = "drop") %>%
  filter(shell %in% c("first_shell", "second_shell"))

## ── Composite DMS effect values ──────────────────────────────────────────────
dms <- read_csv("/Users/mkh/GitHub/mor_dms_analysis/dms_scores/composite_dms_scores.csv",
                show_col_types = FALSE)

## Ligand order (mirrors interface_mean_effect_bar.R, FSK removed since it
## is not an orthosteric ligand)
ligand_order <- c("FSK", "Naloxone", "Naltrexone", "Nalbuphine", "Buprenorphine",
                  "Butorphanol", "MP", "TRV130", "PZM21", "Methadone",
                  "C6guano", "Morphine", "SR17018", "Fentanyl",
                  "Carfentanil", "DAMGO")

## Per-position, per-drug mean missense effect
effect_col_for <- function(drug) {
  col <- paste0(drug, "_effect")
  if (col %in% names(dms)) return(col)
  alt <- paste0(tolower(substr(drug,1,1)), substr(drug,2,nchar(drug)), "_effect")
  if (alt %in% names(dms)) return(alt)
  NA_character_
}

mean_effect <- map_dfr(ligand_order, function(drug) {
  col <- effect_col_for(drug)
  if (is.na(col)) return(NULL)
  dms %>%
    filter(type == "missense") %>%
    group_by(position) %>%
    summarise(mean_effect = mean(.data[[col]], na.rm = TRUE), .groups = "drop") %>%
    mutate(drug = drug)
})

## Synonymous baseline per drug → midpoint of fill scale
syn_baseline <- map_dfr(ligand_order, function(drug) {
  col <- effect_col_for(drug)
  if (is.na(col)) return(NULL)
  tibble(drug = drug,
         syn_effect = mean(dms[[col]][dms$type == "synonymous"], na.rm = TRUE))
})
midpoint_val <- mean(syn_baseline$syn_effect, na.rm = TRUE)

## ── GPCRdb labels & wt residue ──────────────────────────────────────────────
gpcrdb_tbl <- read_csv("/Users/mkh/GitHub/mor_dms_analysis/annotations/GPCRdb_OPRM1_table.csv",
                       show_col_types = FALSE) %>%
  select(position = pos, wt_aa, gpcrdb = GPCRdb, sse = SSE) %>%
  mutate(gpcrdb = if_else(is.na(gpcrdb) | gpcrdb == "", sse, gpcrdb),
         gpcrdb = replace_na(gpcrdb, ""),
         ## GPCRdb generic numbers use a true multiplication sign (1x29 -> 1<d7>29);
         ## only the digit-x-digit form is swapped, so SSE fallbacks are untouched.
         gpcrdb = str_replace(gpcrdb, "(?<=[0-9])x(?=[0-9])", "×"))

## ── Theme (matches gprotein_contact_map.R) ──────────────────────────────────
lw <- 0.5 / 2.835
base_theme <- theme_classic(base_size = 6) +
  theme(
    text             = element_text(color = "black", size = 6),
    axis.text        = element_text(color = "black", size = 6),
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title       = element_text(color = "black", size = 6),
    axis.line        = element_blank(),
    axis.ticks       = element_line(color = "black", linewidth = lw),
    panel.grid.major = element_line(color = "#e8e8e8", linewidth = lw * 0.5),
    plot.title       = element_text(color = "black", size = 6, hjust = 0.5),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.text      = element_text(size = 6),
    legend.title     = element_text(size = 6),
    legend.key.size  = unit(0.25, "cm")
  )

## ── Plot builder ────────────────────────────────────────────────────────────
make_orthosteric_map <- function(shell_label, title) {

  positions <- shell_long %>%
    filter(shell == shell_label) %>%
    pull(position) %>% unique() %>% sort()
  positions <- setdiff(positions, c(142, 219))

  if (length(positions) == 0) {
    message("No positions for shell = ", shell_label)
    return(NULL)
  }

  plot_df <- expand_grid(position = positions, drug = ligand_order) %>%
    left_join(mean_effect, by = c("position", "drug")) %>%
    left_join(gpcrdb_tbl,  by = "position") %>%
    mutate(
      receptor_label = paste0(toupper(wt_aa), position),
      drug = factor(drug, levels = rev(ligand_order))
    )

  rec_meta <- plot_df %>%
    distinct(position, receptor_label, gpcrdb) %>%
    arrange(position)
  rec_order  <- rec_meta$receptor_label
  gpcrdb_vec <- rec_meta$gpcrdb

  plot_df <- plot_df %>%
    mutate(receptor_label = factor(receptor_label, levels = rec_order))

  ## Mark which (position, drug) cells are actually in this shell
  shell_cells <- shell_long %>%
    filter(shell == shell_label, position %in% positions) %>%
    mutate(receptor_label = factor(
      paste0(toupper(gpcrdb_tbl$wt_aa[match(position, gpcrdb_tbl$position)]), position),
      levels = rec_order),
      drug = factor(drug, levels = rev(ligand_order)))

  max_abs <- max(abs(plot_df$mean_effect - midpoint_val), na.rm = TRUE)

  ggplot(plot_df, aes(x = receptor_label, y = drug, fill = mean_effect)) +
    geom_tile(color = "white", linewidth = lw * 0.5) +
    geom_point(data = shell_cells, inherit.aes = FALSE,
               aes(x = receptor_label, y = drug),
               size = 0.6, shape = 16, color = "black") +
    scale_fill_gradient2(
      low      = "#9b1c1d",
      mid      = "white",
      high     = "#3779b9",
      midpoint = midpoint_val,
      limits   = c(midpoint_val - max_abs, midpoint_val + max_abs),
      na.value = "#dddddd",
      name     = "Mean effect"
    ) +
    scale_x_discrete(
      drop = FALSE,
      expand = expansion(add = 0.5),
      sec.axis = dup_axis(name = "GPCRdb", labels = gpcrdb_vec)
    ) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    labs(title = title, x = "MOR residue", y = NULL) +
    base_theme +
    theme(
      axis.text.x.top    = element_text(color = "black", size = 5, angle = 90,
                                        hjust = 0, vjust = 0.5),
      axis.ticks.x.top   = element_line(color = "black", linewidth = lw),
      axis.title.x.top   = element_text(color = "black", size = 6),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = lw)
    )
}

## ── First shell ─────────────────────────────────────────────────────────────
p_first <- make_orthosteric_map("first_shell",
                                "Orthosteric pocket — first-shell residues")
ggsave(file.path(OUT_DIR, "orthosteric_contact_map_first_shell.pdf"),
       p_first, width = 180, height = 88, units = "mm")
message("Saved: orthosteric_contact_map_first_shell.pdf")

## ── Second shell ────────────────────────────────────────────────────────────
p_second <- make_orthosteric_map("second_shell", NULL)
ggsave(file.path(OUT_DIR, "orthosteric_contact_map_second_shell.pdf"),
       p_second, width = 165, height = 50, units = "mm")
message("Saved: orthosteric_contact_map_second_shell.pdf")

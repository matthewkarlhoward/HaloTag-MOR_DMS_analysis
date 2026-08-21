## Compute per-position % LOF / % GOF for each drug from the composite DMS scores.
## Threshold = synonymous mean ± 1 SD per drug.
##
## Inputs:  dms_scores/composite_dms_scores.csv
## Outputs: plots/lof_gof/lof_gof_scores.csv         (per drug × position)
##          plots/lof_gof/lof_gof_group_scores.csv   (averaged within drug groups)
##
## Run: /usr/local/bin/Rscript plots/lof_gof/compute_lof_gof.R

library(tidyverse)

REPO <- "/Users/mkh/GitHub/mor_dms_analysis"
OUT  <- file.path(REPO, "plots/lof_gof")

dms <- read_csv(file.path(REPO, "dms_scores/composite_dms_scores.csv"),
                show_col_types = FALSE)

all_drugs <- c("Buprenorphine", "Butorphanol", "C6guano", "Carfentanil",
               "DAMGO", "Fentanyl", "MP", "Methadone", "Morphine", "Nalbuphine",
               "Naloxone", "Naltrexone", "PZM21", "SR17018", "TRV130")

synonymous <- dms %>% filter(type == "synonymous")
missense   <- dms %>% filter(type == "missense")

# ── Per-drug LOF / GOF ────────────────────────────────────────────────────────
lof_gof <- map_dfr(all_drugs, function(drug) {
  col <- paste0(drug, "_effect")

  syn_vals <- synonymous %>% pull(all_of(col))
  syn_mean <- mean(syn_vals, na.rm = TRUE)
  syn_sd   <- sd(syn_vals,   na.rm = TRUE)

  lof_thresh <- syn_mean - syn_sd
  gof_thresh <- syn_mean + syn_sd

  missense %>%
    filter(!is.na(.data[[col]])) %>%
    group_by(position) %>%
    summarise(
      n_variants = n(),
      pct_lof    = mean(.data[[col]] < lof_thresh) * 100,
      pct_gof    = mean(.data[[col]] > gof_thresh) * 100,
      .groups    = "drop"
    ) %>%
    mutate(drug      = drug,
           lof_thresh = lof_thresh,
           gof_thresh = gof_thresh)
})

write_csv(lof_gof, file.path(OUT, "lof_gof_scores.csv"))
message("Saved: lof_gof_scores.csv")
message(sprintf("  %d rows  |  %d positions  |  %d drugs",
                nrow(lof_gof),
                n_distinct(lof_gof$position),
                n_distinct(lof_gof$drug)))

# ── Ligand groups ─────────────────────────────────────────────────────────────
drug_groups <- tribble(
  ~drug,            ~group,
  "Naltrexone",     "Naltrexone / Naloxone",
  "Naloxone",       "Naltrexone / Naloxone",
  "Nalbuphine",     "Nalbuphine / Buprenorphine / Butorphanol / MP",
  "Buprenorphine",  "Nalbuphine / Buprenorphine / Butorphanol / MP",
  "Butorphanol",    "Nalbuphine / Buprenorphine / Butorphanol / MP",
  "MP",             "Nalbuphine / Buprenorphine / Butorphanol / MP",
  "TRV130",         "TRV130 / PZM21 / Methadone",
  "PZM21",          "TRV130 / PZM21 / Methadone",
  "Methadone",      "TRV130 / PZM21 / Methadone",
  "C6guano",        "C6guano / Morphine",
  "Morphine",       "C6guano / Morphine",
  "SR17018",        "SR17018 / DAMGO / Fentanyl / Carfentanil",
  "DAMGO",          "SR17018 / DAMGO / Fentanyl / Carfentanil",
  "Fentanyl",       "SR17018 / DAMGO / Fentanyl / Carfentanil",
  "Carfentanil",    "SR17018 / DAMGO / Fentanyl / Carfentanil"
)

lof_gof_groups <- lof_gof %>%
  left_join(drug_groups, by = "drug") %>%
  group_by(group, position) %>%
  summarise(mean_pct_lof = mean(pct_lof, na.rm = TRUE),
            mean_pct_gof = mean(pct_gof, na.rm = TRUE),
            n_drugs      = n(),
            .groups      = "drop")

write_csv(lof_gof_groups, file.path(OUT, "lof_gof_group_scores.csv"))
message("Saved: lof_gof_group_scores.csv")
message(sprintf("  %d rows  |  %d positions  |  %d groups",
                nrow(lof_gof_groups),
                n_distinct(lof_gof_groups$position),
                n_distinct(lof_gof_groups$group)))

lof_gof_groups %>%
  group_by(group) %>%
  summarise(mean_lof = round(mean(mean_pct_lof), 2),
            mean_gof = round(mean(mean_pct_gof), 2)) %>%
  print(n = Inf)

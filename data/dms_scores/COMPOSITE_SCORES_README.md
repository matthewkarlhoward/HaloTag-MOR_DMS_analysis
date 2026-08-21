# Composite DMS Score Sheet

## Overview
This composite score sheet (`composite_dms_scores.csv`) combines DMS scores from multiple assays into a single wide-format dataset. Each row represents a unique variant (identified by HGVS notation), and columns contain scores and parameters from different assays.

1. **Fentanyl DRC** (`dms_scores/drc/Fentanyl_wide_hill_ordinal_df.csv`)
     - concentration-specific columns: `MOR_Fentanyl_-12M_rescaled_score`, `MOR_Fentanyl_-12M_rescaled_sd`, etc. (8 concentrations: -12M to -5M)
     - columns with `Fentanyl_` prefix for curve parameters
   - Curve characteristics: curve_type, ec50, slope, emin, emax

2. **Morphine DRC** (`dms_scores/drc/Morphine_wide_hill_ordinal_df.csv`)
   - Same parameter structure as Fentanyl

3. **DAMGO DRC** (`dms_scores/drc/DAMGO_drc.csv`)
     - 12 concentration-specific columns: `DAMGO_mor_dms_6_rescaled_score`, `DAMGO_mor_dms_6_rescaled_sd`, etc. (6 concentrations: dms_1 through dms_6)
     - 12 columns for curve parameters
   - Parameters: curve_type, range, emax, emin, ec50

4. **Surface Expression** (`dms_scores/surface/surface_expression_scores.tsv`)
   - 6 columns with `Surface_` prefix
   - Effect sizes and statistics: effect, effect_se, lfsr, pos_effect, pos_sd, discovery05

5. **Boltz2** (`boltz2_data/mor_filtered_scored.csv`)
   - 17 columns with `Boltz2_` prefix and `_affinity` suffix
   - Protein structural affinity predictions for each drug
   - Drugs: Buprenorphine, Butorphanol, C6guano, DAMGO, FSK, Fentanyl, MP, Methadone, Morphine, Nalbuphine, Naloxone, Naltrexone, PZM21, SR17018, TRV130, carfentanil

6. **Multidrug Panel** (`dms_scores/multidrug/all_drug_scores_mapped_normalized.tsv`)
   - For each drug: effect, effect_se, and lfsr values
   - Drugs included: Buprenorphine, Butorphanol, C6guano, DAMGO, FSK, Fentanyl, MP, Methadone, Morphine, Nalbuphine, Naloxone, Naltrexone, PZM21, SR17018, TRV130, carfentanil

## Annotations
The composite dataset includes several annotation sources merged by position or HGVS:

7. **AlphaMissense** (`annotations/OPRM1_alphamissense.csv`)
   - 2 columns: `alphamissense_score`, `alphamissense_class`
   - Pathogenicity predictions for missense variants
   - Merged on `hgvs` identifier

8. **MTR (Missense Tolerance Ratio)** (`annotations/OPRM1_MTR.csv`)
   - 1 column: `MTR`
   - Regional constraint metric by position
   - Merged on `position`
   - From Regeneron 1M exomes database

9. **Mouse Homolog Position** (`annotations/oprm_human_mouse_pos.csv`)
   - 1 column: `pos_mouse`
   - Maps human positions to mouse OPRM1 positions
   - Merged on `position`

10. **GPCRdb** (`annotations/GPCRdb_OPRM1_table.csv`)
    - 3 columns: `GPCRdb`, `SSE`, `wt_aa`
    - GPCR database numbering and secondary structure elements
    - Merged on `position`

11. **Motifs** (`annotations/motifs.csv`)
    - 3 columns: `Gprotein`, `Sodium`, `Switch`
    - Boolean indicators for functional motif membership
    - Merged on `GPCRdb` identifier
    - taken from GPCRdb

## Column Structure

### Variant Metadata (5 columns)
- `hgvs`: HGVS notation for variant (e.g., p.(D2A))
- `type`: Variant type (missense, deletion, insertion, etc.)
- `position`: Amino acid position
- `wildtype`: Wild-type amino acid
- `mutation`: Mutant amino acid

### Assay-Specific Columns (151 columns)
All assay-specific columns are prefixed with the assay/drug name followed by underscore and the parameter name.

Examples:
- **Concentration-specific scores**:
  - `MOR_Fentanyl_-12M_rescaled_score`: Fentanyl score at 10^-12 M concentration
  - `MOR_Fentanyl_-12M_rescaled_sd`: Standard deviation for above
  - `MOR_Morphine_-115M_rescaled_score`: Morphine score at 10^-11.5 M concentration
  - `DAMGO_mor_dms_6_rescaled_score`: DAMGO score at concentration 6 (highest, 10µM,) - each concentration is a log step down from 10µM: 1µM, 0.1µM, etc. Note that there were only 5 drug concentrations here, so the curve paramaters are a bit wonky and less constrained that for the 8-point curves we have on Morphine and Fentanyl.
  - `DAMGO_mor_dms_1_rescaled_score`: DAMGO score at concentration 1 (lowest 0M)
- **Curve parameters**:
  - `Fentanyl_ec50`: EC50 value from Fentanyl dose-response curve
  - `Morphine_score_avg`: Average score from Morphine assay
  - `DAMGO_curve_type`: Curve classification for DAMGO
  - `DAMGO_ec50`: EC50 value from DAMGO dose-response curve
- **Other assays**:
  - `Surface_effect`: Effect size for surface expression
  - `Boltz2_Fentanyl_affinity`: Boltz2 predicted affinity for Fentanyl
  - `Boltz2_Morphine_affinity`: Boltz2 predicted affinity for Morphine
  - `Methadone_effect`: Effect size for Methadone from multidrug panel

### Annotation Columns (10 columns)
- `alphamissense_score`: AlphaMissense pathogenicity score (0-1, higher = more pathogenic)
- `alphamissense_class`: Classification (benign, ambiguous, pathogenic)
- `MTR`: Missense Tolerance Ratio (regional constraint metric)
- `pos_mouse`: Corresponding mouse OPRM1 position
- `GPCRdb`: GPCRdb generic numbering
- `SSE`: Secondary structure element
- `wt_aa`: Wild-type amino acid from GPCRdb
- `Gprotein`: Boolean — part of G-protein coupling motif
- `Sodium`: Boolean — part of sodium binding motif
- `Switch`: Boolean — part of molecular switch motif


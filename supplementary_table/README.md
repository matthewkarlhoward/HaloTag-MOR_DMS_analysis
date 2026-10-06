# Consolidated supplementary data table

`MOR_DMS_supplementary_data.xlsx` — one workbook, nine sheets plus a `00_README`
index. `csv/` holds the same nine sheets as flat CSVs for programmatic use.
Rebuild with `python3 build_supplementary_table.py` (reads only from `../data/`).

| Sheet | Grain | Rows × cols |
|---|---|---|
| `01_variant_scores` | variant (hgvs) | 10,178 × 109 |
| `02_dose_response_scores` | variant (hgvs) | 10,178 × 49 |
| `03_position_summary` | receptor position | 399 × 42 |
| `04_position_lof_gof` | position × ligand | 5,985 × 8 |
| `05_ligand_summary` | ligand | 15 × 20 |
| `06_ligand_contacts` | position × ligand | 4,517 × 6 |
| `07_validation_pharmacology` | variant × ligand (fitted) | 109 × 20 |
| `08_screen_samples` | sorted bin | 532 × 12 |
| `09_validation_trupath_points` | individual replicate point | 6,335 × 11 |

Grains are deliberately separate. Variants, positions, ligands and sorted samples are
different units of observation; flattening them into one file would either duplicate
rows or strand columns.

---

## Conventions that apply throughout

**Sign.** The reporter is signal-down: μOR activation *inhibits* cAMP and *lowers* eGFP.
All DMS scores are centred at zero for synonymous variants, with **negative = loss of
signalling (or loss of surface expression)** and **positive = gain**.

**Emax scale.** Hill `emax` is on an activity scale normalised per ligand so that the
mean synonymous variant has `emin = 0` and `emax = 1`. It is *not* a percentage.

**EC50.** Reported as log10 M throughout (`logEC50_M`), never as a molar concentration.

**Positions.** Human OPRM1 numbering, 1–400. `pos_mouse` gives the mouse equivalent
(human − 2). `GPCRdb` is the Ballesteros–Weinstein / GPCRdb generic number.

**Ligand names.** Normalised to: TRV130 (not Oliceridine), Carfentanil, MP (mitragynine
pseudoindoxyl), C6guano, SR17018.

**Efficacy classes**, from the unsupervised clustering in figure 5b and concordant with
TRUPATH Emax: antagonist (naltrexone, naloxone); weak (nalbuphine, buprenorphine,
butorphanol, MP); intermediate (TRV130, PZM21, methadone); strong (C6guano, morphine,
SR17018, fentanyl, carfentanil, DAMGO).

**Ligand-free conditions.** There are two, and they are not interchangeable:
`NoLigand_effect` (`FSK_effect` upstream) is the forskolin-only arm of the 15-ligand
saturating panel; `NoLigand_DAMGOscreen_score` is the forskolin-only arm of the DAMGO
concentration-response screen (`DAMGO_mor_dms_6` upstream, which is **not** a DAMGO
concentration despite the column name).

---

## 01_variant_scores

| Column(s) | Meaning |
|---|---|
| `hgvs` | Variant identifier, the join key across every sheet |
| `type` | missense / synonymous / insertion / deletion |
| `position`, `wildtype`, `mutation` | Human position, WT residue, substituted residue |
| `GPCRdb`, `SSE`, `pos_mouse` | Generic number, secondary structure element, mouse position |
| `surface_effect`, `_effect_se`, `surface_lfsr` | HaloTag/JF635i surface expression score, standard error, local false sign rate |
| `<Ligand>_effect`, `_effect_se`, `_lfsr` | Saturating (10 µM) signalling score per ligand, from Lilace |
| `NoLigand_effect`, `_effect_se`, `_lfsr` | Forskolin-only arm of the same panel |
| `<Morphine\|Fentanyl\|DAMGO>_hill_emin` | Fitted basal activity |
| `..._hill_logEC50_M` | Fitted potency, log10 M |
| `..._hill_emax` | Fitted efficacy on the synonymous-normalised activity scale |
| `..._hill_r2`, `_n_points`, `_n_outliers_dropped` | Fit quality and how many points entered the refit |
| `..._hill_curve_type` | `sigmoid` where both asymptotes were reached; other values mark flat or unconstrained fits |
| `..._hill_fit_ok` | Convergence flag. **Filter on `curve_type == "sigmoid"` and `fit_ok` before using EC50 or Emax** |
| `..._op_logtau`, `_op_logKA`, `_op_logEC50_M`, `_op_emax_obs`, `_op_rmse` | Black–Leff operational fit |
| `Morphine_op_log10_rho` | Expression-corrected per-receptor coupling, log10, WT = 0 |
| `..._op_saturating`, `_op_KA_unidentified` | Flags: the curve saturated, so K_A is not identifiable |
| `PCA_variant_PC1..3` | Variant-level PCA scores (PC1 49.1%, PC2 10.2% of variance) |
| `alphamissense_score`, `alphamissense_class`, `MTR` | External constraint annotation, not generated here |

Only morphine carries `log10_rho`; fentanyl and DAMGO were fit with the operational model
for concordance only (supplemental figure 5), without the expression decomposition.

## 02_dose_response_scores

`<Ligand>_score_logM<c>` and `<Ligand>_sd_logM<c>` for each measured concentration:
morphine −11.5 to −4.5 in 1-log steps, fentanyl −12 to −5, DAMGO −9 to −5. Plus
`NoLigand_DAMGOscreen_score` / `_sd`. These are the Lilace rescaled scores that the Hill
fits in sheet 01 were fit to.

## 03_position_summary

Position identity and annotation; `n_missense`; `<Ligand>_mean_missense_effect` (the
position average used in most figures); `surface_mean_missense_effect`;
`<Ligand>_mean_logEC50_M` and `_mean_emax` for the three DRC ligands; position-level PCA
scores (`PCA_position_PC*`, PC1 66.8%, PC2 15.6%) and the per-position mean of the
variant-level PCA (`PCA_variantmean_PC*`, the values painted on the structure in figure
5g); `EC50bias_*` and `Emaxbias_*` (morphine minus fentanyl, with FDR-corrected q and a
significance flag — figures 4c and 4e); and `gprotein_contact_class`
(Active / Intermediate / Both / blank) from the first shells of the active and
latent-state structures.

## 04_position_lof_gof

`pct_lof` and `pct_gof` are the fraction of missense variants at that position falling
beyond `lof_thresh` / `gof_thresh`, which are the synonymous mean ± 2 SD for that ligand.
`efficacy_class` is carried alongside for convenience.

## 05_ligand_summary

`camp_*` and `tp_*` are the wild-type pooled-cAMP and TRUPATH Gi1 parameters as reported
in supplemental figure 1d: `logec50`, `emax` (% DAMGO), their SEMs, and n. **The two Emax
values are computed differently** — the cAMP value is `Span`/DAMGO `Span` × 100 from a
three-parameter Hill fit with slope fixed at 1, while the TRUPATH value is the peak of
the per-dose mean activity as a percentage of the DAMGO peak, not a fit parameter.

`PCA_position_loading_PC*` and `PCA_variant_loading_PC*` are the ligand loadings from the
two PCA fits. Figure 5f uses the variant-level PC2 loading. Loadings are sign-aligned to
the position-level fit.

SR17018 has no wild-type concentration-response data in either assay (NA), and
naltrexone and C6guano have TRUPATH only.

## 06_ligand_contacts

One row per position × ligand where a contact was assigned, from experimental structures
where available and Chai-1 predictions otherwise (`structure` names the source). First
shell = **side-chain atoms within 5.0 Å** of any ligand atom; second shell = not first
shell and **side-chain atoms within 4.5 Å** of a first-shell residue's side chain. Where a
ligand has several structures the closest approach is kept. This matches the definition
given in the manuscript Methods.

## 07_validation_pharmacology

Fitted parameters for **all** low-throughput TRUPATH Gi1 BRET in the paper: the
single-mutant series behind figure 7d–f and the A119L double-mutant series behind
figures 8d/f/g and supplemental figure 11a. 109 curves over the four ligands and twenty
variants that appear in those figures.
The individual points behind every one of these fits are in sheet `09`.

`pipeline` says which analysis produced the row, and this matters because **the two
pipelines normalise differently**:

| `pipeline` | Normalisation | Figures |
|---|---|---|
| `singles` | each curve divided by its own fitted no-drug plateau | 7d–f |
| `doubles` | the two runs co-scaled on their shared WT and A119L arms, not renormalised to DAMGO | 8d, 8f, 8g, S11a |

`activation_window` = −Span, the depth of the signal-down curve.
`pct_activation_vs_WT_DAMGO` is computed **within each source**, against that source's own
WT + DAMGO window, so values are comparable inside a pipeline but should not be compared
across pipelines without care. `responsive` is the extra-sum-of-squares F test against a
flat line; where it is False the window is reported as 0 for the doubles. A small negative
percentage means a significant but upward, non-activating curve, not negative efficacy.

> **Scope.** Both sheets are restricted to what the paper plots: DAMGO, PZM21 and
> nalbuphine (figure 7d–f) plus naloxone (figures 8d/f/g, S11a). The singles pipeline
> also measured buprenorphine, butorphanol, fentanyl, MP and morphine against these
> variants, and two further substitutions (I280D, I280P), none of which appear in any
> figure; those are not deposited.

> **One raw run appears twice, under both normalisations.** The 20260427 experiment is
> `source = doubles_raw` in the singles pipeline and part of `doubles1+doubles2` in the
> doubles pipeline. These are the same wells analysed two ways, not independent
> measurements. Filter on `pipeline` before aggregating.

## 08_screen_samples

Per-sorted-bin metadata for every DMS screen: assay, date, replicate, ligand and
concentration, bin, and mean/median sequencing coverage. Supports supplemental figures
2h–i and documents which raw count file each bin corresponds to.

## 09_validation_trupath_points

Every individual replicate point behind sheet `07`, one row per
pipeline × source × ligand × variant × concentration × replicate. 6,335 points.

| Column | Meaning |
|---|---|
| `pipeline`, `source`, `run` | Which analysis and which experimental run |
| `ligand`, `variant`, `is_double` | What was measured |
| `log10_conc_M`, `replicate` | Dose and biological replicate |
| `raw_bret` | Unnormalised BRET ratio (515/410). **Singles only** |
| `normalized` | The value plotted. Singles: raw ÷ `divisor`. Doubles: the co-scaled response extracted from Prism |
| `divisor` | The fitted no-drug plateau used to normalise that curve. Singles only |

`raw_bret` and `divisor` are empty for the doubles because those points were extracted
from the Prism projects after normalisation; the unnormalised values are not recoverable
from the deposited files.

---

## Decisions worth confirming before submission

1. **MCAM is excluded.** The analysis master table carries a full set of MCAM
   (methocinnamox) columns, but MCAM is not one of the 15 ligands in this paper.
   Excluded from every sheet. Reinstate only if it is being added to the manuscript.
2. **Boltz2 predicted affinities are excluded.** Present in the master table, not used
   in any figure. Easy to add back as a block in sheet 01 if wanted.
3. **Legacy Lilace ordinal curve parameters are excluded.** The master table also holds
   an earlier set of `Morphine_ec50` / `_emax` / `_slope` columns from the ordinal
   Lilace fit. The paper reports the robust three-parameter refit, so only that is
   carried here, to avoid readers picking the wrong column.
4. **Sheet 07 covers the doubles only.** The single-mutant BRET behind figures 7d–f
   (K100D/N, V175E/N, R278D) and the surface-expression flow data behind supplemental
   figure 10b live in Prism projects that are not in this repository. If those points
   are exported to CSV in the same shape as
   `data/pharmacology/doubles/doubles_merged_points.csv`, they drop straight into this
   sheet and the table becomes complete.
5. **Workbook size is 16 MB.** Under most journal limits but worth checking against the
   target journal's supplementary file cap; the CSVs can ship instead or in addition.

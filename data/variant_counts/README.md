# variant_counts/

Per-variant read counts for every sorted bin in every screen reported in the paper.
These are the primary data: Lilace consumes them to produce the scores in
`../dms_scores/`, which everything else in this repository is derived from.

**532 count files across five screens.** One file per sorted bin, per replicate, per
condition. Each row of `../../supplementary_table/csv/08_screen_samples.csv` corresponds
to exactly one file here, and its `variant_counts_filename` column is the join key.

| Folder | Screen | Files | Size |
|---|---|---|---|
| `surface/` | HaloTag/JF635i surface expression | 12 | 4.4 MB |
| `damgo_drc/` | DAMGO concentration-response | 72 | 26 MB |
| `morphine_drc/` | Morphine concentration-response | 108 | 39 MB |
| `fentanyl_drc/` | Fentanyl concentration-response | 108 | 39 MB |
| `multidrug/` | 15-ligand panel at 10 µM saturating | 232 | 85 MB |
| `_metadata/` | Sequencing sample sheets and cell counts | 6 | 240 KB |

## File format

Every count file is a CSV with one row per observed library variant:

| Column | Meaning |
|---|---|
| `count` | Reads observed for this variant in this bin |
| `pos` | Position in µOR (human numbering) |
| `mutation_type` | `M` missense, `S` synonymous, `I` insertion, `D` deletion |
| `name` | Short variant label, e.g. `A10del` |
| `codon` | Codon called, where applicable |
| `mutation` | Substitution or indel identifier |
| `length` | Indel length in amino acids (1 for substitutions) |
| `hgvs` | HGVS protein notation — **the join key** to every other table in this repository |

Roughly 10,400 rows per file, matching the designed library.

## Naming conventions

Two conventions are present, reflecting when each screen was run.

**Surface, DAMGO** — `MOR_<condition>_R<replicate>_B<bin>.csv`
`B0`–`B3` are the four sort bins, each 25% of the population.
For `damgo_drc/`, `MOR_1` … `MOR_6` are the six DAMGO arms, where `MOR_1` is 10 µM
(−5 log M) descending to `MOR_5` at −9 log M, and **`MOR_6` is the ligand-free
(forskolin-only) arm** — not a DAMGO concentration, despite the name.

**Morphine, fentanyl, multidrug** — `a<date>_rep<n>_<ligand>_<conc>_sample_<n>_bin_<a-d>_S<n>_merged.csv`
Here bins are lettered `a`–`d` rather than numbered `B0`–`B3`, and `FSK_0M` marks the
ligand-free arm.

Because the two conventions disagree, use `08_screen_samples.csv` rather than parsing
filenames — it carries assay, date, replicate, ligand, concentration and bin as proper
columns, already reconciled.

## Upstream of these files

FASTQ → adapter removal (BBDuk) → paired-read error correction (BBMerge) → mapping
(BBMap, 15-mers) → variant calling (GATK `AnalyzeSaturationMutagenesis`) → filtering to
designed variants. Run through the **Dumpling** pipeline, released as part of Rosace.
Raw FASTQ are deposited separately; see the paper's Data availability section.

## What is deliberately not here

Counts for screens **not reported in this paper** — internalization, abundance, and
several additional arms — are excluded, as is `dms_data/` from the analysis repository,
which holds the same reads in a different, partly redundant layout. The five folders
above are the complete set behind every published figure.

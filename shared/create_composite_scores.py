#!/usr/bin/env python3
"""
Create the master composite DMS score sheet (wide format, one row per variant).

Combines:
  - Dose-response curves (Fentanyl, Morphine, DAMGO): EC50, Emax, Emin, slope,
    fit parameters, and per-concentration rescaled scores
  - Multidrug panel raw effects (mapped=='raw' AND lilace_run==drug, sign flipped)
  - Surface expression scores
  - Boltz2 affinity predictions (per drug)
  - Annotations: AlphaMissense, MTR, mouse position, GPCRdb, motifs
  - Ligand distance and shell info (per drug, per position) from /structures

Output: dms_scores/composite_dms_scores.csv
"""

import re
import pandas as pd
import numpy as np

# ---------- File paths ----------
fentanyl_file  = 'dms_scores/drc/Fentanyl_wide_hill_ordinal_df.csv'
morphine_file  = 'dms_scores/drc/Morphine_wide_hill_ordinal_df.csv'
damgo_file     = 'dms_scores/drc/DAMGO_drc.csv'
surface_file   = 'dms_scores/surface/surface_expression_scores.tsv'
multidrug_file = 'dms_scores/multidrug/all_drug_scores_mapped_normalized.tsv'
boltz2_file    = 'boltz2_data/mor_filtered_scored.csv'

alphamissense_file = 'annotations/OPRM1_alphamissense.csv'
mtr_file           = 'annotations/OPRM1_MTR.csv'
mouse_pos_file     = 'annotations/oprm_human_mouse_pos.csv'
gpcrdb_file        = 'annotations/GPCRdb_OPRM1_table.csv'
motifs_file        = 'annotations/motifs.csv'

ligand_dist_files = [
    'structures/processed/ligand_distances/experimental_and_chai_combined.csv',
    'structures/processed/ligand_distances/all_structures_ligand_distances.csv',
]

# DMS-tested ligands we want distance/shell coverage for.
# FSK (forskolin control) has no structure.
TESTED_LIGANDS = {
    'Buprenorphine', 'Butorphanol', 'C6guano', 'DAMGO', 'Fentanyl',
    'Methadone', 'Mitragynine_Pseudoindoxyl', 'Morphine', 'Nalbuphine',
    'Naloxone', 'Naltrexone', 'TRV130', 'PZM21', 'SR17018', 'Carfentanil',
}

# Canonical drug name remapping applied to all sources before pivoting
DRUG_RENAME = {
    'Oliceridine': 'TRV130',
    'C6-guano':    'C6guano',
    'carfentanil': 'Carfentanil',
}

output_file = 'dms_scores/composite_dms_scores.csv'

# ---------- Dose-response curves ----------
print("Loading dose-response curves...")

metadata_cols = ['hgvs', 'type', 'position', 'wildtype', 'mutation']
drc_param_cols = ['a', 'b', 'm', 's', 'convergence', 'score_avg', 'sd_avg',
                  'curve_type', 'ec50', 'slope', 'emin', 'emax']

print("  Fentanyl...")
fentanyl_df = pd.read_csv(fentanyl_file)
f_score_cols = [c for c in fentanyl_df.columns if 'MOR_Fentanyl' in c and 'rescaled' in c]
fentanyl_df = fentanyl_df[metadata_cols + f_score_cols + drc_param_cols].copy()
fentanyl_df.columns = metadata_cols + f_score_cols + [f'Fentanyl_{c}' for c in drc_param_cols]

print("  Morphine...")
morphine_df = pd.read_csv(morphine_file)
m_score_cols = [c for c in morphine_df.columns if 'MOR_Morphine' in c and 'rescaled' in c]
morphine_df = morphine_df[['hgvs'] + m_score_cols + drc_param_cols].copy()
morphine_df.columns = ['hgvs'] + m_score_cols + [f'Morphine_{c}' for c in drc_param_cols]

print("  DAMGO...")
damgo_df = pd.read_csv(damgo_file)
# NOTE: the raw DAMGO_drc.csv has `emin` and `emax` labels swapped relative
# to the Morphine / Fentanyl convention. Swap them here so all three drugs
# use the same convention (emin = baseline, emax = drug-induced extreme).
damgo_df = damgo_df.rename(columns={'emin': '__tmp', 'emax': 'emin'}) \
                   .rename(columns={'__tmp': 'emax'})
d_score_cols = [c for c in damgo_df.columns if 'mor_dms' in c and 'rescaled' in c]
d_param_cols = ['curve_type', 'range', 'emax', 'emin', 'ec50', 'slope',
                'score_avg', 'sigmoid_PC1', 'sigmoid_PC2']
damgo_df = damgo_df[['hgvs'] + d_score_cols + d_param_cols].copy()
damgo_df.columns = ['hgvs'] + [f'DAMGO_{c}' for c in d_score_cols] + [f'DAMGO_{c}' for c in d_param_cols]

# ---------- Surface expression ----------
print("Loading surface expression...")
surface_df = pd.read_csv(surface_file, sep='\t')
surface_df = surface_df[['variant', 'effect', 'effect_se', 'lfsr',
                         'pos_effect', 'pos_sd', 'discovery05']].copy()
surface_df.columns = ['hgvs', 'Surface_effect', 'Surface_effect_se', 'Surface_lfsr',
                      'Surface_pos_effect', 'Surface_pos_sd', 'Surface_discovery05']

# ---------- Boltz2 affinity ----------
print("Loading Boltz2 affinity...")
boltz2_df = pd.read_csv(boltz2_file)[['variant', 'drug', 'Boltz2_affinity']].copy()
boltz2_df.columns = ['hgvs', 'drug', 'Boltz2_affinity']
boltz2_df['drug'] = boltz2_df['drug'].replace(DRUG_RENAME)
boltz2_df = boltz2_df.groupby(['hgvs', 'drug'], as_index=False)['Boltz2_affinity'].mean()
boltz2_wide = boltz2_df.pivot(index='hgvs', columns='drug', values='Boltz2_affinity')
boltz2_wide.columns = [f'Boltz2_{drug}_affinity' for drug in boltz2_wide.columns]
boltz2_wide = boltz2_wide.reset_index()

# ---------- Multidrug panel (raw effects, sign flipped) ----------
print("Loading multidrug panel...")
multidrug_df = pd.read_csv(multidrug_file, sep='\t')
multidrug_df = multidrug_df[
    (multidrug_df['mapped'] == 'raw') &
    (multidrug_df['lilace_run'] == multidrug_df['drug'])
].copy()
multidrug_df = multidrug_df[['hgvs', 'drug', 'effect', 'effect_se', 'lfsr']].copy()
multidrug_df['drug'] = multidrug_df['drug'].replace(DRUG_RENAME)
multidrug_df['effect'] = -multidrug_df['effect']  # sign flip
multidrug_df = multidrug_df.groupby(['hgvs', 'drug'], as_index=False).agg({
    'effect': 'mean', 'effect_se': 'mean', 'lfsr': 'mean'
})
multidrug_wide = multidrug_df.pivot(index='hgvs', columns='drug',
                                    values=['effect', 'effect_se', 'lfsr'])
multidrug_wide.columns = [f'{drug}_{metric}' for metric, drug in multidrug_wide.columns]
multidrug_wide = multidrug_wide.reset_index()

# ---------- Annotations ----------
print("Loading annotations...")

alphamissense_df = pd.read_csv(alphamissense_file, encoding='utf-8-sig')
alphamissense_df.columns = ['hgvs', 'alphamissense_score', 'alphamissense_class']

mtr_df = pd.read_csv(mtr_file, encoding='utf-8-sig')
mtr_df.columns = ['MTR', 'position']

mouse_pos_df = pd.read_csv(mouse_pos_file, encoding='utf-8-sig')
mouse_pos_df.columns = ['pos_mouse', 'position']
mouse_pos_df = mouse_pos_df.dropna()

gpcrdb_df = pd.read_csv(gpcrdb_file, encoding='utf-8-sig')
gpcrdb_df.columns = ['GPCRdb', 'position', 'SSE', 'wt_aa']

motifs_df = pd.read_csv(motifs_file, encoding='utf-8-sig')
motifs_df.columns = ['GPCRdb', 'Gprotein', 'Sodium', 'Switch']
for col in ['Gprotein', 'Sodium', 'Switch']:
    motifs_df[col] = motifs_df[col].notna()
motifs_df = motifs_df.groupby('GPCRdb', as_index=False).agg(
    {'Gprotein': 'max', 'Sodium': 'max', 'Switch': 'max'}
)

# ---------- Ligand distance / shell ----------
# Union of available distance files (later files override earlier on drug+position).
# Provides per-drug per-position distance, sidechain distance, and shell label.
print("Loading ligand distance / shell data...")
ligand_keep = ['pdb_id', 'drug_name', 'residue_number_human',
               'distance_angstrom', 'distance_sidechain_angstrom', 'shell']
ligand_parts = []
for f in ligand_dist_files:
    df = pd.read_csv(f, usecols=lambda c: c in ligand_keep)
    ligand_parts.append(df)
ligand_long = pd.concat(ligand_parts, ignore_index=True)
ligand_long = ligand_long.dropna(subset=['drug_name', 'residue_number_human'])

# Apply canonical drug name remapping
ligand_long['drug_name'] = ligand_long['drug_name'].replace(DRUG_RENAME)

# Carfentanil has no experimental structure; use Lofentanil/7t2h as proxy
carfentanil_proxy = ligand_long[
    (ligand_long['drug_name'] == 'Lofentanil') &
    (ligand_long.get('pdb_id', '7t2h').astype(str).str.lower() == '7t2h')
].copy() if 'pdb_id' in ligand_long.columns else ligand_long[
    ligand_long['drug_name'] == 'Lofentanil'
].copy()
carfentanil_proxy['drug_name'] = 'Carfentanil'
ligand_long = pd.concat([ligand_long, carfentanil_proxy], ignore_index=True)

# Filter to DMS-tested ligands only
ligand_long = ligand_long[ligand_long['drug_name'].isin(TESTED_LIGANDS)]

# Within each drug+position, take minimum distance (closest contact across PDBs)
ligand_long = ligand_long.sort_values('distance_angstrom').drop_duplicates(
    subset=['drug_name', 'residue_number_human'], keep='first'
)
ligand_wide = ligand_long.pivot(
    index='residue_number_human', columns='drug_name',
    values=['distance_angstrom', 'distance_sidechain_angstrom', 'shell']
)
ligand_wide.columns = [f'{drug}_{metric}' for metric, drug in ligand_wide.columns]
ligand_wide = ligand_wide.reset_index().rename(columns={'residue_number_human': 'position'})
print(f"  Ligand drugs: {sorted(ligand_long['drug_name'].unique())}")

# ---------- Merge ----------
print("\nMerging...")
composite_df = fentanyl_df.copy()
composite_df = composite_df.merge(morphine_df, on='hgvs', how='outer')
composite_df = composite_df.merge(damgo_df, on='hgvs', how='outer')
composite_df = composite_df.merge(surface_df, on='hgvs', how='outer')
composite_df = composite_df.merge(boltz2_wide, on='hgvs', how='outer')
composite_df = composite_df.merge(multidrug_wide, on='hgvs', how='outer')

# Backfill metadata from hgvs notation for variants not in fentanyl base
print("Backfilling metadata from HGVS...")
missing_mask = composite_df['type'].isna()
print(f"  {missing_mask.sum()} rows with missing metadata")

for idx in composite_df[missing_mask].index:
    hgvs = composite_df.at[idx, 'hgvs']
    m = re.match(r'p\.\(([A-Z])(\d+)([A-Z])\)', hgvs)
    if m:
        wt, pos, mut = m.groups()
        composite_df.at[idx, 'type'] = 'synonymous' if wt == mut else 'missense'
        composite_df.at[idx, 'position'] = int(pos)
        composite_df.at[idx, 'wildtype'] = wt
        composite_df.at[idx, 'mutation'] = mut
        continue
    m = re.match(r'p\.\(([A-Z])(\d+)del\)', hgvs)
    if m:
        wt, pos = m.groups()
        composite_df.at[idx, 'type'] = 'deletion'
        composite_df.at[idx, 'position'] = int(pos)
        composite_df.at[idx, 'wildtype'] = wt
        composite_df.at[idx, 'mutation'] = 'del'
        continue
    m = re.match(r'p\.\(([A-Z])(\d+)_([A-Z])(\d+)del\)', hgvs)
    if m:
        wt1, pos1, wt2, pos2 = m.groups()
        composite_df.at[idx, 'type'] = 'deletion'
        composite_df.at[idx, 'position'] = int(pos1)
        composite_df.at[idx, 'wildtype'] = wt1
        composite_df.at[idx, 'mutation'] = f'del{int(pos2) - int(pos1)}'
        continue
    m = re.match(r'p\.\(([A-Z])(\d+)_([A-Z])(\d+)ins([A-Z]+)\)', hgvs)
    if m:
        wt1, pos1, wt2, pos2, ins_seq = m.groups()
        composite_df.at[idx, 'type'] = 'insertion'
        composite_df.at[idx, 'position'] = int(pos1)
        composite_df.at[idx, 'wildtype'] = wt1
        composite_df.at[idx, 'mutation'] = f'ins{ins_seq}'
        continue
    print(f"  WARNING: Could not parse hgvs: {hgvs}")

still_missing = composite_df['type'].isna().sum()
print(f"  {missing_mask.sum() - still_missing} backfilled, {still_missing} still missing")

# Annotation merges (position- or hgvs-based)
print("Merging annotations...")
composite_df = composite_df.merge(alphamissense_df, on='hgvs', how='left')
composite_df = composite_df.merge(mtr_df, on='position', how='left')
composite_df = composite_df.merge(mouse_pos_df, on='position', how='left')
composite_df = composite_df.merge(gpcrdb_df[['position', 'GPCRdb', 'SSE', 'wt_aa']],
                                  on='position', how='left')
composite_df = composite_df.merge(motifs_df, on='GPCRdb', how='left')
for col in ['Gprotein', 'Sodium', 'Switch']:
    composite_df[col] = composite_df[col].fillna(False)

# Ligand distance / shell merge (position-based)
print("Merging ligand distance / shell...")
composite_df = composite_df.merge(ligand_wide, on='position', how='left')

# ---------- Save ----------
composite_df = composite_df.sort_values(['position', 'hgvs']).reset_index(drop=True)
print(f"\nSaving to: {output_file}")
composite_df.to_csv(output_file, index=False)

print(f"\nComplete!")
print(f"  Total variants: {len(composite_df)}")
print(f"  Total columns: {len(composite_df.columns)}")
print(f"\nColumn groups:")
print(f"  Metadata: hgvs, type, position, wildtype, mutation")
print(f"  Fentanyl DRC: {len([c for c in composite_df.columns if c.startswith('Fentanyl_')])}")
print(f"  Morphine DRC: {len([c for c in composite_df.columns if c.startswith('Morphine_')])}")
print(f"  DAMGO DRC: {len([c for c in composite_df.columns if c.startswith('DAMGO_')])}")
print(f"  Surface: {len([c for c in composite_df.columns if c.startswith('Surface_')])}")
print(f"  Boltz2: {len([c for c in composite_df.columns if c.startswith('Boltz2_')])}")
print(f"  Ligand distance/shell: {len([c for c in composite_df.columns if 'distance_angstrom' in c or c.endswith('_shell')])}")

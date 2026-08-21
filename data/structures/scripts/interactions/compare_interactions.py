#!/usr/bin/env python3
"""
Compare non-bonded interactions across MOR structure states.

Three analyses:
  1. Active vs Inactive: Per-residue comparison of interaction frequencies
     - Mean number of each interaction type per residue
     - Identifies residues that gain/lose interactions upon activation
  2. Within-Active variance: How variable is each residue's interaction count
  3. Within-Inactive variance: Same for inactive state

Also analyses specific contact pairs:
  - Which residue-residue contacts are state-specific?
  - Which salt bridges / H-bonds are conserved vs variable?

Input: per_residue_interaction_counts.csv and all_*.csv from extract_interactions.py
Output: comparison CSVs and summary report
"""

import os
import csv
import math
import numpy as np
from collections import defaultdict, Counter

from pathlib import Path
REPO_ROOT = Path(__file__).resolve().parents[3]
CIF_DIR = REPO_ROOT / "structures" / "raw" / "experimental"

# Structures excluded from state comparisons (per-structure data kept)
EXCLUDED_FROM_COMPARISON = {"9WSV", "9WSX"}


def load_per_residue_counts(csv_file):
    """
    Load per-residue interaction counts, excluding arrestin-bound structures.
    Returns: {resnum: {"active": [list of count dicts], "inactive": [list of count dicts]}}
    """
    data = defaultdict(lambda: {"active": [], "inactive": []})

    with open(csv_file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["pdb_id"] in EXCLUDED_FROM_COMPARISON:
                continue
            resnum = int(row["residue_number_human"])
            state = row["state"]
            counts = {
                "pdb_id": row["pdb_id"],
                "hbond_total": int(row["hbond_total"]),
                "hbond_bb_bb": int(row["hbond_bb_bb"]),
                "hbond_bb_sc": int(row["hbond_bb_sc"]),
                "hbond_sc_sc": int(row["hbond_sc_sc"]),
                "salt_bridge": int(row["salt_bridge"]),
                "hydrophobic": int(row["hydrophobic"]),
            }
            data[resnum][state].append(counts)

    return data


def load_contact_pairs(csv_file, interaction_type):
    """
    Load pairwise contacts. Returns list of dicts.
    """
    pairs = []
    with open(csv_file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pairs.append(row)
    return pairs


def mean_std(values):
    """Compute mean and std of a list of numbers."""
    if not values:
        return None, None
    m = np.mean(values)
    s = np.std(values, ddof=1) if len(values) > 1 else 0
    return round(float(m), 3), round(float(s), 3)


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = str(REPO_ROOT / "structures" / "processed" / "interactions")

    counts_file = os.path.join(results_dir, "per_residue_interaction_counts.csv")
    if not os.path.exists(counts_file):
        print(f"ERROR: {counts_file} not found. Run extract_interactions.py first.")
        return

    print("Loading per-residue interaction counts...")
    data = load_per_residue_counts(counts_file)
    print(f"  {len(data)} residue positions loaded")

    interaction_types = ["hbond_total", "hbond_bb_bb", "hbond_bb_sc", "hbond_sc_sc",
                         "salt_bridge", "hydrophobic"]

    # ── Analysis 1: Active vs Inactive per-residue comparison ────────────────
    print("\n1. Active vs Inactive comparison...")
    comparison_rows = []

    for resnum in sorted(data.keys()):
        d = data[resnum]
        row = {"residue_number_human": resnum}

        for itype in interaction_types:
            active_vals = [c[itype] for c in d["active"]]
            inactive_vals = [c[itype] for c in d["inactive"]]

            act_mean, act_std = mean_std(active_vals)
            ina_mean, ina_std = mean_std(inactive_vals)

            row[f"{itype}_active_mean"] = act_mean if act_mean is not None else ""
            row[f"{itype}_active_std"] = act_std if act_std is not None else ""
            row[f"{itype}_active_n"] = len(active_vals)
            row[f"{itype}_inactive_mean"] = ina_mean if ina_mean is not None else ""
            row[f"{itype}_inactive_std"] = ina_std if ina_std is not None else ""
            row[f"{itype}_inactive_n"] = len(inactive_vals)

            if act_mean is not None and ina_mean is not None:
                row[f"{itype}_delta"] = round(act_mean - ina_mean, 3)
            else:
                row[f"{itype}_delta"] = ""

        comparison_rows.append(row)

    # Build field names
    comp_fields = ["residue_number_human"]
    for itype in interaction_types:
        comp_fields.extend([
            f"{itype}_active_mean", f"{itype}_active_std", f"{itype}_active_n",
            f"{itype}_inactive_mean", f"{itype}_inactive_std", f"{itype}_inactive_n",
            f"{itype}_delta"
        ])

    comp_file = os.path.join(results_dir, "active_vs_inactive_interactions.csv")
    with open(comp_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=comp_fields)
        writer.writeheader()
        writer.writerows(comparison_rows)
    print(f"  Written: {comp_file}")

    # ── Analysis 2 & 3: Within-state variance ────────────────────────────────
    print("\n2. Within-state interaction variance...")

    for state_name in ["active", "inactive"]:
        var_rows = []
        for resnum in sorted(data.keys()):
            d = data[resnum]
            entries = d[state_name]
            if not entries:
                continue

            row = {"residue_number_human": resnum, "n_structures": len(entries)}
            for itype in interaction_types:
                vals = [c[itype] for c in entries]
                m, s = mean_std(vals)
                row[f"{itype}_mean"] = m if m is not None else ""
                row[f"{itype}_std"] = s if s is not None else ""
                # Coefficient of variation (only if mean > 0)
                if m and m > 0 and s is not None:
                    row[f"{itype}_cv"] = round(s / m, 3)
                else:
                    row[f"{itype}_cv"] = ""

            var_rows.append(row)

        var_fields = ["residue_number_human", "n_structures"]
        for itype in interaction_types:
            var_fields.extend([f"{itype}_mean", f"{itype}_std", f"{itype}_cv"])

        var_file = os.path.join(results_dir, f"{state_name}_within_state_interactions.csv")
        with open(var_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=var_fields)
            writer.writeheader()
            writer.writerows(var_rows)
        print(f"  Written: {var_file}")

    # ── Analysis 4: Contact pair conservation ────────────────────────────────
    print("\n3. Analyzing contact pair conservation...")

    # Load pairwise contacts
    for contact_type, filename in [
        ("hbond_sc_sc", "all_hbonds.csv"),
        ("salt_bridge", "all_salt_bridges.csv"),
        ("hydrophobic", "all_hydrophobic_contacts.csv"),
    ]:
        master_file = os.path.join(results_dir, filename)
        if not os.path.exists(master_file):
            print(f"  WARNING: {filename} not found, skipping {contact_type}")
            continue

        # Count how many structures each residue pair appears in, per state
        pair_counts = defaultdict(lambda: {"active": set(), "inactive": set()})
        active_pdbs = set()
        inactive_pdbs = set()

        with open(master_file, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row["pdb_id"] in EXCLUDED_FROM_COMPARISON:
                    continue
                # For H-bonds, only consider sc-sc
                if contact_type == "hbond_sc_sc" and row.get("hb_type", "") != "sc-sc":
                    continue

                r1 = int(row["res1_num"])
                r2 = int(row["res2_num"])
                pair_key = (min(r1, r2), max(r1, r2))
                pdb_id = row["pdb_id"]
                state = row["state"]

                pair_counts[pair_key][state].add(pdb_id)
                if state == "active":
                    active_pdbs.add(pdb_id)
                else:
                    inactive_pdbs.add(pdb_id)

        n_active = len(active_pdbs)
        n_inactive = len(inactive_pdbs)

        # Write contact pair conservation
        pair_rows = []
        for (r1, r2), counts in sorted(pair_counts.items()):
            n_act = len(counts["active"])
            n_ina = len(counts["inactive"])
            act_freq = round(n_act / n_active, 3) if n_active > 0 else 0
            ina_freq = round(n_ina / n_inactive, 3) if n_inactive > 0 else 0

            # Classification
            if act_freq >= 0.5 and ina_freq >= 0.5:
                conservation = "conserved"
            elif act_freq >= 0.3 and ina_freq < 0.2:
                conservation = "active_specific"
            elif ina_freq >= 0.3 and act_freq < 0.2:
                conservation = "inactive_specific"
            else:
                conservation = "variable"

            pair_rows.append({
                "res1_num": r1,
                "res2_num": r2,
                "active_count": n_act,
                "active_freq": act_freq,
                "inactive_count": n_ina,
                "inactive_freq": ina_freq,
                "delta_freq": round(act_freq - ina_freq, 3),
                "conservation": conservation,
            })

        pair_file = os.path.join(results_dir, f"{contact_type}_pair_conservation.csv")
        with open(pair_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=[
                "res1_num", "res2_num",
                "active_count", "active_freq",
                "inactive_count", "inactive_freq",
                "delta_freq", "conservation"
            ])
            writer.writeheader()
            # Sort by absolute delta_freq descending
            pair_rows.sort(key=lambda r: abs(r["delta_freq"]), reverse=True)
            writer.writerows(pair_rows)
        print(f"  {pair_file}: {len(pair_rows)} unique pairs")

    # ── Summary report ───────────────────────────────────────────────────────
    print("\n4. Generating summary report...")
    summary_file = os.path.join(results_dir, "comparison_summary.txt")

    with open(summary_file, "w") as f:
        f.write("Interaction Comparison Summary\n")
        f.write("=" * 70 + "\n\n")

        # ── Active vs Inactive: top changing residues ──
        f.write("ACTIVE vs INACTIVE — Top residues by change in interaction count\n")
        f.write("-" * 70 + "\n\n")

        for itype in ["hbond_sc_sc", "salt_bridge", "hydrophobic"]:
            # Filter rows with data in both states
            rows_with_delta = [r for r in comparison_rows
                               if r[f"{itype}_delta"] != ""
                               and r[f"{itype}_active_n"] >= 5
                               and r[f"{itype}_inactive_n"] >= 2]
            rows_with_delta.sort(key=lambda r: abs(r[f"{itype}_delta"]), reverse=True)

            f.write(f"\n{itype.upper()} — Top 20 residues by |Δmean|:\n")
            f.write(f"{'Res#':>5s}  {'Act_mean':>8s} {'Act_std':>7s}  "
                    f"{'Ina_mean':>8s} {'Ina_std':>7s}  {'Delta':>7s}\n")
            f.write("-" * 55 + "\n")
            for r in rows_with_delta[:20]:
                f.write(f"{r['residue_number_human']:>5d}  "
                        f"{str(r[f'{itype}_active_mean']):>8s} {str(r[f'{itype}_active_std']):>7s}  "
                        f"{str(r[f'{itype}_inactive_mean']):>8s} {str(r[f'{itype}_inactive_std']):>7s}  "
                        f"{str(r[f'{itype}_delta']):>7s}\n")

        # ── Contact pair conservation ──
        f.write(f"\n\n{'='*70}\n")
        f.write("STATE-SPECIFIC CONTACT PAIRS\n")
        f.write("-" * 70 + "\n\n")

        for contact_type in ["hbond_sc_sc", "salt_bridge", "hydrophobic"]:
            pair_file = os.path.join(results_dir, f"{contact_type}_pair_conservation.csv")
            if not os.path.exists(pair_file):
                continue

            pairs = []
            with open(pair_file, "r") as pf:
                reader = csv.DictReader(pf)
                for row in reader:
                    pairs.append(row)

            active_spec = [p for p in pairs if p["conservation"] == "active_specific"]
            inactive_spec = [p for p in pairs if p["conservation"] == "inactive_specific"]
            conserved = [p for p in pairs if p["conservation"] == "conserved"]

            f.write(f"\n{contact_type.upper()}:\n")
            f.write(f"  Total unique pairs: {len(pairs)}\n")
            f.write(f"  Conserved (both states ≥50%): {len(conserved)}\n")
            f.write(f"  Active-specific (≥30% active, <20% inactive): {len(active_spec)}\n")
            f.write(f"  Inactive-specific (≥30% inactive, <20% active): {len(inactive_spec)}\n\n")

            if active_spec:
                f.write(f"  Top active-specific pairs:\n")
                f.write(f"  {'Res1':>5s} {'Res2':>5s}  {'Act_freq':>8s} {'Ina_freq':>8s}\n")
                f.write("  " + "-" * 35 + "\n")
                for p in active_spec[:15]:
                    f.write(f"  {p['res1_num']:>5s} {p['res2_num']:>5s}  "
                            f"{p['active_freq']:>8s} {p['inactive_freq']:>8s}\n")

            if inactive_spec:
                f.write(f"\n  Top inactive-specific pairs:\n")
                f.write(f"  {'Res1':>5s} {'Res2':>5s}  {'Act_freq':>8s} {'Ina_freq':>8s}\n")
                f.write("  " + "-" * 35 + "\n")
                for p in inactive_spec[:15]:
                    f.write(f"  {p['res1_num']:>5s} {p['res2_num']:>5s}  "
                            f"{p['active_freq']:>8s} {p['inactive_freq']:>8s}\n")

        # ── Within-state variance highlights ──
        f.write(f"\n\n{'='*70}\n")
        f.write("WITHIN-STATE VARIANCE HIGHLIGHTS\n")
        f.write("-" * 70 + "\n\n")

        for state_name in ["active", "inactive"]:
            var_file = os.path.join(results_dir, f"{state_name}_within_state_interactions.csv")
            if not os.path.exists(var_file):
                continue

            var_rows = []
            with open(var_file, "r") as vf:
                reader = csv.DictReader(vf)
                for row in reader:
                    var_rows.append(row)

            f.write(f"\n{state_name.upper()} state — Most variable residues (by hbond_sc_sc std):\n")
            sc_rows = [r for r in var_rows if r["hbond_sc_sc_std"] and float(r["hbond_sc_sc_std"]) > 0]
            sc_rows.sort(key=lambda r: float(r["hbond_sc_sc_std"]), reverse=True)

            f.write(f"  {'Res#':>5s}  {'Mean':>6s} {'Std':>6s} {'n':>3s}\n")
            f.write("  " + "-" * 30 + "\n")
            for r in sc_rows[:15]:
                f.write(f"  {r['residue_number_human']:>5s}  "
                        f"{r['hbond_sc_sc_mean']:>6s} {r['hbond_sc_sc_std']:>6s} "
                        f"{r['n_structures']:>3s}\n")

    print(f"  Written: {summary_file}")
    print("\nDone! All comparison results in results/")


if __name__ == "__main__":
    main()

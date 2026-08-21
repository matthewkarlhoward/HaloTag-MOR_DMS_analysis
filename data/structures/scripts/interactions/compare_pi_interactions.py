#!/usr/bin/env python3
"""
Compare pi interaction counts between active and inactive MOR structures.

Loads per_residue_pi_counts.csv produced by extract_pi_interactions.py and
computes, for each residue position, the mean count across active structures,
the mean count across inactive structures, and the delta (active - inactive)
for each pi interaction type.

Output: active_vs_inactive_pi_interactions.csv
"""

import csv
import os
from collections import defaultdict

from pathlib import Path
REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = REPO_ROOT / "structures" / "processed" / "pi_interactions"


def load_per_residue_counts(filepath):
    """
    Load per_residue_pi_counts.csv.
    Returns: dict keyed by (pdb_id, residue_number_human) ->
             {"state": str, "pi_pi_count": int, "pi_cation_count": int, "pi_ligand_count": int}
    """
    records = []
    with open(filepath, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            records.append({
                "pdb_id": row["pdb_id"],
                "state": row["state"],
                "residue_number_human": int(row["residue_number_human"]),
                "pi_pi_count": int(row["pi_pi_count"]),
                "pi_cation_count": int(row["pi_cation_count"]),
                "pi_ligand_count": int(row["pi_ligand_count"]),
            })
    return records


def compute_active_vs_inactive(records):
    """
    For each residue position, accumulate counts from active and inactive
    structures separately, then compute means and delta.

    Returns: list of dicts with per-residue comparison statistics,
             sorted by residue_number_human.
    """
    # data[resnum][state] = list of (pi_pi, pi_cation, pi_ligand) tuples
    data = defaultdict(lambda: {"active": [], "inactive": []})

    for rec in records:
        resnum = rec["residue_number_human"]
        state = rec["state"]
        if state not in ("active", "inactive"):
            continue
        data[resnum][state].append((
            rec["pi_pi_count"],
            rec["pi_cation_count"],
            rec["pi_ligand_count"],
        ))

    results = []
    for resnum in sorted(data.keys()):
        active_vals   = data[resnum]["active"]
        inactive_vals = data[resnum]["inactive"]

        def safe_mean(vals, idx):
            if not vals:
                return None
            return round(sum(v[idx] for v in vals) / len(vals), 4)

        def safe_delta(active_mean, inactive_mean):
            if active_mean is None or inactive_mean is None:
                return None
            return round(active_mean - inactive_mean, 4)

        pi_pi_active   = safe_mean(active_vals,   0)
        pi_pi_inactive = safe_mean(inactive_vals, 0)
        pi_pi_delta    = safe_delta(pi_pi_active, pi_pi_inactive)

        pi_cation_active   = safe_mean(active_vals,   1)
        pi_cation_inactive = safe_mean(inactive_vals, 1)
        pi_cation_delta    = safe_delta(pi_cation_active, pi_cation_inactive)

        pi_ligand_active   = safe_mean(active_vals,   2)
        pi_ligand_inactive = safe_mean(inactive_vals, 2)
        pi_ligand_delta    = safe_delta(pi_ligand_active, pi_ligand_inactive)

        results.append({
            "residue_number_human": resnum,
            "pi_pi_active_mean":      pi_pi_active,
            "pi_pi_inactive_mean":    pi_pi_inactive,
            "pi_pi_delta":            pi_pi_delta,
            "pi_cation_active_mean":  pi_cation_active,
            "pi_cation_inactive_mean":pi_cation_inactive,
            "pi_cation_delta":        pi_cation_delta,
            "pi_ligand_active_mean":  pi_ligand_active,
            "pi_ligand_inactive_mean":pi_ligand_inactive,
            "pi_ligand_delta":        pi_ligand_delta,
        })

    return results


def main():
    counts_file = RESULTS_DIR / "per_residue_pi_counts.csv"
    output_file = RESULTS_DIR / "active_vs_inactive_pi_interactions.csv"

    if not counts_file.exists():
        raise FileNotFoundError(
            f"Input file not found: {counts_file}\n"
            "Run extract_pi_interactions.py first."
        )

    print(f"Loading: {counts_file}")
    records = load_per_residue_counts(str(counts_file))
    print(f"  Loaded {len(records)} residue-structure records")

    n_active   = sum(1 for r in records if r["state"] == "active")
    n_inactive = sum(1 for r in records if r["state"] == "inactive")
    print(f"  Active records:   {n_active}")
    print(f"  Inactive records: {n_inactive}")

    print("\nComputing active vs inactive means and deltas...")
    results = compute_active_vs_inactive(records)
    print(f"  Residue positions compared: {len(results)}")

    fieldnames = [
        "residue_number_human",
        "pi_pi_active_mean",
        "pi_pi_inactive_mean",
        "pi_pi_delta",
        "pi_cation_active_mean",
        "pi_cation_inactive_mean",
        "pi_cation_delta",
        "pi_ligand_active_mean",
        "pi_ligand_inactive_mean",
        "pi_ligand_delta",
    ]

    with open(str(output_file), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"\nWrote: {output_file}")
    print("Done!")


if __name__ == "__main__":
    main()

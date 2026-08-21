#!/usr/bin/env python3
"""
Compare chi angles across MOR structure states.

Three analyses:
  1. Active vs Inactive: Per-residue comparison of chi angle distributions
     - Mean circular difference in χ1, χ2
     - Rotamer state frequencies (g+, g-, t) per state
     - Identifies residues with significant rotameric shifts
  2. Within-Active variance: Circular variance of chi angles among active structures
  3. Within-Inactive variance: Circular variance of chi angles among inactive structures

Uses circular statistics (angles are periodic: -180° wraps to +180°).

Input: all_chi_angles_combined.csv from extract_chi_angles.py
Output: comparison CSVs and summary report
"""

import os
import csv
import math
import numpy as np
from collections import defaultdict

from pathlib import Path
REPO_ROOT = Path(__file__).resolve().parents[3]
CIF_DIR = REPO_ROOT / "structures" / "raw" / "experimental"


def circular_mean(angles_deg):
    """Compute circular mean of angles (in degrees). Returns degrees in [-180, 180]."""
    if not angles_deg:
        return None
    angles_rad = [math.radians(a) for a in angles_deg]
    sin_sum = sum(math.sin(a) for a in angles_rad)
    cos_sum = sum(math.cos(a) for a in angles_rad)
    mean_rad = math.atan2(sin_sum / len(angles_rad), cos_sum / len(angles_rad))
    return round(math.degrees(mean_rad), 2)


def circular_variance(angles_deg):
    """
    Compute circular variance (0 = no spread, 1 = uniform/max spread).
    V = 1 - R, where R is the mean resultant length.
    """
    if len(angles_deg) < 2:
        return None
    angles_rad = [math.radians(a) for a in angles_deg]
    sin_sum = sum(math.sin(a) for a in angles_rad)
    cos_sum = sum(math.cos(a) for a in angles_rad)
    R = math.sqrt(sin_sum**2 + cos_sum**2) / len(angles_rad)
    return round(1 - R, 4)


def circular_std(angles_deg):
    """
    Circular standard deviation in degrees.
    σ = sqrt(-2 * ln(R)) converted to degrees.
    """
    if len(angles_deg) < 2:
        return None
    angles_rad = [math.radians(a) for a in angles_deg]
    sin_sum = sum(math.sin(a) for a in angles_rad)
    cos_sum = sum(math.cos(a) for a in angles_rad)
    R = math.sqrt(sin_sum**2 + cos_sum**2) / len(angles_rad)
    if R < 1e-10:
        return 180.0  # Essentially uniform
    if R > 1.0:
        R = 1.0
    return round(math.degrees(math.sqrt(-2 * math.log(R))), 2)


def circular_diff(angle1, angle2):
    """Shortest angular difference between two angles (in degrees). Returns [-180, 180]."""
    if angle1 is None or angle2 is None:
        return None
    diff = angle1 - angle2
    while diff > 180:
        diff -= 360
    while diff < -180:
        diff += 360
    return round(diff, 2)


def classify_rotamer(angle):
    """Classify χ1 into g+, g-, t."""
    if angle is None:
        return None
    while angle > 180:
        angle -= 360
    while angle < -180:
        angle += 360
    if -120 <= angle < 0:
        return "g-"
    elif 0 <= angle < 120:
        return "g+"
    else:
        return "t"


# Structures excluded from state comparisons (per-structure data kept)
EXCLUDED_FROM_COMPARISON = {"9WSV", "9WSX"}


def load_combined_data(csv_file):
    """
    Load combined chi angle data, excluding arrestin-bound structures.
    Returns:
        dict: {resnum: {"resname": str, "active": [{chi1, chi2, chi3, chi4, pdb_id}, ...],
                                         "inactive": [...]}}
    """
    data = defaultdict(lambda: {"resname": None, "active": [], "inactive": []})

    with open(csv_file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["pdb_id"] in EXCLUDED_FROM_COMPARISON:
                continue
            resnum = int(row["residue_number_human"])
            state = row["state"]
            resname = row["residue_name"]

            data[resnum]["resname"] = resname

            entry = {
                "pdb_id": row["pdb_id"],
                "chi1": float(row["chi1"]) if row["chi1"] else None,
                "chi2": float(row["chi2"]) if row["chi2"] else None,
                "chi3": float(row["chi3"]) if row["chi3"] else None,
                "chi4": float(row["chi4"]) if row["chi4"] else None,
            }
            data[resnum][state].append(entry)

    return data


def compute_rotamer_frequencies(entries, chi_key="chi1"):
    """Compute rotamer state frequencies from a list of entries."""
    counts = {"g+": 0, "g-": 0, "t": 0}
    total = 0
    for e in entries:
        rot = classify_rotamer(e[chi_key])
        if rot:
            counts[rot] += 1
            total += 1
    if total == 0:
        return {"g+": 0, "g-": 0, "t": 0, "n": 0}
    freqs = {k: round(v / total, 3) for k, v in counts.items()}
    freqs["n"] = total
    return freqs


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = str(REPO_ROOT / "structures" / "processed" / "rotamers")
    combined_csv = os.path.join(results_dir, "all_chi_angles_combined.csv")

    if not os.path.exists(combined_csv):
        print(f"ERROR: {combined_csv} not found. Run extract_chi_angles.py first.")
        return

    print("Loading combined chi angle data...")
    data = load_combined_data(combined_csv)
    print(f"  {len(data)} residue positions loaded")

    # Count structures
    active_pdbs = set()
    inactive_pdbs = set()
    for resnum, d in data.items():
        for e in d["active"]:
            active_pdbs.add(e["pdb_id"])
        for e in d["inactive"]:
            inactive_pdbs.add(e["pdb_id"])
    print(f"  Active structures: {len(active_pdbs)}")
    print(f"  Inactive structures: {len(inactive_pdbs)}")

    # ── Analysis 1: Active vs Inactive comparison ────────────────────────────
    print("\n1. Computing active vs inactive comparison...")
    comparison_rows = []

    for resnum in sorted(data.keys()):
        d = data[resnum]
        resname = d["resname"]

        for chi_idx in ["chi1", "chi2"]:
            active_angles = [e[chi_idx] for e in d["active"] if e[chi_idx] is not None]
            inactive_angles = [e[chi_idx] for e in d["inactive"] if e[chi_idx] is not None]

            if not active_angles and not inactive_angles:
                continue

            active_mean = circular_mean(active_angles) if active_angles else None
            inactive_mean = circular_mean(inactive_angles) if inactive_angles else None
            active_cstd = circular_std(active_angles) if len(active_angles) >= 2 else None
            inactive_cstd = circular_std(inactive_angles) if len(inactive_angles) >= 2 else None
            mean_diff = circular_diff(active_mean, inactive_mean)

            # Rotamer frequencies (chi1 and chi2)
            active_rot = compute_rotamer_frequencies(d["active"], chi_idx)
            inactive_rot = compute_rotamer_frequencies(d["inactive"], chi_idx)
            active_dominant = max(["g+", "g-", "t"], key=lambda k: active_rot[k]) if active_rot["n"] > 0 else ""
            inactive_dominant = max(["g+", "g-", "t"], key=lambda k: inactive_rot[k]) if inactive_rot["n"] > 0 else ""
            rotamer_shift = "YES" if (active_dominant and inactive_dominant and active_dominant != inactive_dominant) else "NO"

            comparison_rows.append({
                "residue_number_human": resnum,
                "residue_name": resname,
                "chi_angle": chi_idx,
                "active_n": len(active_angles),
                "inactive_n": len(inactive_angles),
                "active_circular_mean": active_mean if active_mean is not None else "",
                "inactive_circular_mean": inactive_mean if inactive_mean is not None else "",
                "mean_difference": mean_diff if mean_diff is not None else "",
                "abs_mean_difference": abs(mean_diff) if mean_diff is not None else "",
                "active_circular_std": active_cstd if active_cstd is not None else "",
                "inactive_circular_std": inactive_cstd if inactive_cstd is not None else "",
                "active_dominant_rotamer": active_dominant,
                "inactive_dominant_rotamer": inactive_dominant,
                "rotamer_shift": rotamer_shift,
                "active_gp_freq": active_rot["g+"],
                "active_gm_freq": active_rot["g-"],
                "active_t_freq": active_rot["t"],
                "inactive_gp_freq": inactive_rot["g+"],
                "inactive_gm_freq": inactive_rot["g-"],
                "inactive_t_freq": inactive_rot["t"],
            })

    comp_file = os.path.join(results_dir, "active_vs_inactive_chi_comparison.csv")
    with open(comp_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "residue_number_human", "residue_name", "chi_angle",
            "active_n", "inactive_n",
            "active_circular_mean", "inactive_circular_mean",
            "mean_difference", "abs_mean_difference",
            "active_circular_std", "inactive_circular_std",
            "active_dominant_rotamer", "inactive_dominant_rotamer", "rotamer_shift",
            "active_gp_freq", "active_gm_freq", "active_t_freq",
            "inactive_gp_freq", "inactive_gm_freq", "inactive_t_freq",
        ])
        writer.writeheader()
        writer.writerows(comparison_rows)
    print(f"  Written: {comp_file}")

    # ── Analysis 2 & 3: Within-state variance ────────────────────────────────
    print("\n2. Computing within-state chi angle variance...")

    for state_name, state_key in [("active", "active"), ("inactive", "inactive")]:
        variance_rows = []

        for resnum in sorted(data.keys()):
            d = data[resnum]
            resname = d["resname"]
            entries = d[state_key]

            for chi_idx in ["chi1", "chi2"]:
                angles = [e[chi_idx] for e in entries if e[chi_idx] is not None]
                if not angles:
                    continue

                c_mean = circular_mean(angles)
                c_var = circular_variance(angles)
                c_std = circular_std(angles)
                n = len(angles)

                # Rotamer frequencies
                rot_freq = compute_rotamer_frequencies(entries, chi_idx) if chi_idx == "chi1" else {"g+": "", "g-": "", "t": "", "n": ""}
                dominant = ""
                mixed = ""
                if chi_idx == "chi1" and rot_freq["n"] > 0:
                    dominant = max(["g+", "g-", "t"], key=lambda k: rot_freq[k])
                    # "mixed" if dominant rotamer < 70% of observations
                    max_freq = max(rot_freq["g+"], rot_freq["g-"], rot_freq["t"])
                    if isinstance(max_freq, (int, float)):
                        mixed = "YES" if max_freq < 0.7 else "NO"

                variance_rows.append({
                    "residue_number_human": resnum,
                    "residue_name": resname,
                    "chi_angle": chi_idx,
                    "n_structures": n,
                    "circular_mean": c_mean if c_mean is not None else "",
                    "circular_variance": c_var if c_var is not None else "",
                    "circular_std": c_std if c_std is not None else "",
                    "dominant_rotamer": dominant,
                    "mixed_rotamers": mixed,
                    "gp_freq": rot_freq["g+"],
                    "gm_freq": rot_freq["g-"],
                    "t_freq": rot_freq["t"],
                })

        var_file = os.path.join(results_dir, f"{state_name}_within_state_chi_variance.csv")
        with open(var_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=[
                "residue_number_human", "residue_name", "chi_angle",
                "n_structures", "circular_mean", "circular_variance", "circular_std",
                "dominant_rotamer", "mixed_rotamers",
                "gp_freq", "gm_freq", "t_freq",
            ])
            writer.writeheader()
            writer.writerows(variance_rows)
        print(f"  Written: {var_file}")

    # ── Summary report ───────────────────────────────────────────────────────
    print("\n3. Generating summary report...")
    summary_file = os.path.join(results_dir, "comparison_summary.txt")

    # Find top rotamer shifts
    chi1_shifts = [r for r in comparison_rows
                   if r["chi_angle"] == "chi1" and r["rotamer_shift"] == "YES"
                   and r["abs_mean_difference"] != ""]
    chi1_shifts.sort(key=lambda r: abs(float(r["abs_mean_difference"])), reverse=True)

    # Find largest absolute differences (chi1 only, with enough data)
    chi1_big_diff = [r for r in comparison_rows
                     if r["chi_angle"] == "chi1"
                     and r["abs_mean_difference"] != ""
                     and r["active_n"] >= 5 and r["inactive_n"] >= 2]
    chi1_big_diff.sort(key=lambda r: abs(float(r["abs_mean_difference"])), reverse=True)

    with open(summary_file, "w") as f:
        f.write("Chi Angle Comparison Summary\n")
        f.write("=" * 70 + "\n\n")

        # ── Active vs Inactive ──
        f.write("ACTIVE vs INACTIVE COMPARISON\n")
        f.write("-" * 70 + "\n\n")

        n_shifts = len(chi1_shifts)
        f.write(f"Residues with χ1 rotamer state shift (active vs inactive): {n_shifts}\n\n")

        if chi1_shifts:
            f.write("Top 20 rotamer shifts (sorted by |Δmean|):\n")
            f.write(f"{'Res#':>5s} {'Name':>4s}  {'Active':>7s} {'Inact':>7s}  "
                    f"{'ΔMean':>7s}  {'Act_rot':>7s} {'Ina_rot':>7s}\n")
            f.write("-" * 60 + "\n")
            for r in chi1_shifts[:20]:
                f.write(f"{r['residue_number_human']:>5}  {r['residue_name']:>4s}  "
                        f"{str(r['active_circular_mean']):>7s} {str(r['inactive_circular_mean']):>7s}  "
                        f"{str(r['mean_difference']):>7s}  "
                        f"{r['active_dominant_rotamer']:>7s} {r['inactive_dominant_rotamer']:>7s}\n")

        f.write(f"\n\nTop 30 largest |Δχ1 mean| (active_n≥5, inactive_n≥2):\n")
        f.write(f"{'Res#':>5s} {'Name':>4s}  {'Active':>7s} {'Inact':>7s}  "
                f"{'ΔMean':>7s}  {'|ΔMean|':>7s}  {'Act_std':>7s} {'Ina_std':>7s}\n")
        f.write("-" * 70 + "\n")
        for r in chi1_big_diff[:30]:
            f.write(f"{r['residue_number_human']:>5}  {r['residue_name']:>4s}  "
                    f"{str(r['active_circular_mean']):>7s} {str(r['inactive_circular_mean']):>7s}  "
                    f"{str(r['mean_difference']):>7s}  {str(r['abs_mean_difference']):>7s}  "
                    f"{str(r['active_circular_std']):>7s} {str(r['inactive_circular_std']):>7s}\n")

        # ── Within-state variance ──
        f.write(f"\n\n{'='*70}\n")
        f.write("WITHIN-STATE CHI ANGLE VARIANCE\n")
        f.write("-" * 70 + "\n\n")

        for state_name in ["active", "inactive"]:
            var_file = os.path.join(results_dir, f"{state_name}_within_state_chi_variance.csv")
            rows = []
            with open(var_file, "r") as vf:
                reader = csv.DictReader(vf)
                for row in reader:
                    rows.append(row)

            chi1_rows = [r for r in rows if r["chi_angle"] == "chi1" and r["circular_variance"]]
            chi1_rows.sort(key=lambda r: float(r["circular_variance"]), reverse=True)

            mixed_rows = [r for r in chi1_rows if r.get("mixed_rotamers") == "YES"]

            f.write(f"\n{state_name.upper()} state:\n")
            f.write(f"  Residues with χ1 data: {len(chi1_rows)}\n")
            f.write(f"  Residues with mixed rotamers (<70% dominant): {len(mixed_rows)}\n\n")

            f.write(f"  Top 20 most variable χ1 positions:\n")
            f.write(f"  {'Res#':>5s} {'Name':>4s}  {'CirVar':>7s} {'CirStd':>7s}  "
                    f"{'n':>3s}  {'g+':>5s} {'g-':>5s} {'t':>5s}  {'Mixed':>5s}\n")
            f.write("  " + "-" * 60 + "\n")
            for r in chi1_rows[:20]:
                f.write(f"  {r['residue_number_human']:>5}  {r['residue_name']:>4s}  "
                        f"{r['circular_variance']:>7s} {r['circular_std']:>7s}  "
                        f"{r['n_structures']:>3s}  "
                        f"{r['gp_freq']:>5s} {r['gm_freq']:>5s} {r['t_freq']:>5s}  "
                        f"{r.get('mixed_rotamers', ''):>5s}\n")

    print(f"  Written: {summary_file}")
    print("\nDone! All comparison results in results/")


if __name__ == "__main__":
    main()

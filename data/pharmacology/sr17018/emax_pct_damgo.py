#!/usr/bin/env python3
"""SR-17018 TruPath Emax as % DAMGO, matched on the 20260502 plate.

Applies the house TruPath recipe (see reference_trupath_emax_normalization):
  1. per replicate, divide the whole series by its own no-drug baseline -> norm
  2. activity = 1 - norm            (the assay is signal-down)
  3. Emax = MAX of the per-dose MEAN activity   (peak, not a fitted Top/Span)
  4. % DAMGO = Emax_ligand / Emax_DAMGO x 100

Replicate structure on this plate: the dose series runs across columns in pairs
(1-2 = 1e-4 ... 21-22 = 1e-14, 23-24 = vehicle), so each plate row carries TWO
replicate series, the odd and the even column of each pair. DAMGO has 3 rows
(6 replicates), buprenorphine and SR-17018 have 6 rows (12 replicates).

Buprenorphine is the cross-plate control: it also appears in the published
14-ligand panel at 60.2% DAMGO.

Reads only sr17018_trupath_wells.csv, written by extract_sr17018_trupath.py.
"""

import csv
import statistics as st
from pathlib import Path

HERE = Path(__file__).resolve().parent
WELLS = HERE / "sr17018_trupath_wells.csv"

PRISM_READ = "READ2"
LIGAND_ROWS = {"DAMGO": "AFK", "Buprenorphine": "BCGHLM", "SR17018": "DEIJNO"}
DOSES = [-4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14]
VEHICLE = 0
PUBLISHED_BUP_PCT_DAMGO = 60.2  # Data 4 tab, 20260117 panel


def load():
    """-> {ligand: {(row, parity): {dose: ratio}}} for the read the .prism uses."""
    out = {}
    with open(WELLS) as fh:
        for r in csv.DictReader(fh):
            if r["read"] != PRISM_READ or not r["ligand"] or r["ligand"] == "blank":
                continue
            col = int(r["plate_col"])
            key = (r["plate_row"], "odd" if col % 2 else "even")
            out.setdefault(r["ligand"], {}).setdefault(key, {})[
                int(r["log_dose_M"])] = float(r["bret_ratio"])
    return out


def emax(series, baseline_dose):
    """Peak of the per-dose mean activity, plus the per-replicate values at that peak.

    baseline_dose: a single dose used per-replicate as the no-drug anchor, or the
    string 'plateau' to use each replicate's mean over 1e-7..1e-14 (more robust,
    reported as a sensitivity check only).
    """
    activity = {}
    for key, curve in series.items():
        if baseline_dose == "plateau":
            base = st.mean(curve[d] for d in (-7, -8, -9, -10, -11, -12, -13, -14))
        else:
            base = curve[baseline_dose]
        for dose, val in curve.items():
            activity.setdefault(dose, {})[key] = 1.0 - val / base
    means = {d: st.mean(v.values()) for d, v in activity.items() if d in DOSES}
    peak_dose = max(means, key=means.get)
    return peak_dose, means[peak_dose], activity[peak_dose], means


def main():
    data = load()
    print("Replicate structure (read %s):" % PRISM_READ)
    for lig, series in data.items():
        print(f"  {lig:14s} rows {LIGAND_ROWS[lig]:7s} -> {len(series)} replicate series")

    for baseline, label in [(-13, "own 1e-13 (house recipe)"),
                            (-14, "own 1e-14, most dilute dose here"),
                            (VEHICLE, "own vehicle well"),
                            ("plateau", "own 1e-7..1e-14 mean")]:
        res = {lig: emax(s, baseline) for lig, s in data.items()}
        d_peak = res["DAMGO"][1]
        print(f"\n=== baseline: {label} ===")
        print(f"  {'ligand':14s} {'peak at':>9s} {'Emax':>9s} {'% DAMGO':>9s} {'SEM':>7s}")
        for lig in ("DAMGO", "Buprenorphine", "SR17018"):
            pd_, pm, per_rep, _ = res[lig]
            vals = list(per_rep.values())
            sem_pct = (st.stdev(vals) / len(vals) ** 0.5) / d_peak * 100
            print(f"  {lig:14s} {('1e%d' % pd_):>9s} {pm:9.4f} "
                  f"{pm / d_peak * 100:9.1f} {sem_pct:7.1f}")
        bup = res["Buprenorphine"][1] / d_peak * 100
        print(f"  cross-plate control: buprenorphine {bup:.1f}% here vs "
              f"{PUBLISHED_BUP_PCT_DAMGO}% published "
              f"({bup - PUBLISHED_BUP_PCT_DAMGO:+.1f} pts)")

    # write the sensitivity table
    out = []
    for baseline, label in [(-13, "own 1e-13 (house recipe)"),
                            (-14, "own 1e-14"), (VEHICLE, "own vehicle well"),
                            ("plateau", "own 1e-7..1e-14 mean")]:
        res = {lig: emax(s, baseline) for lig, s in data.items()}
        d_peak = res["DAMGO"][1]
        for lig in ("DAMGO", "Buprenorphine", "SR17018"):
            pd_, pm, per_rep, _ = res[lig]
            vals = list(per_rep.values())
            out.append({
                "baseline": label, "ligand": lig, "peak_log_dose_M": pd_,
                "emax_activity": round(pm, 6),
                "pct_damgo": round(pm / d_peak * 100, 1),
                "pct_damgo_sem": round(
                    (st.stdev(vals) / len(vals) ** 0.5) / d_peak * 100, 1),
                "n_replicates": len(vals),
            })
    with open(HERE / "sr17018_emax_pct_damgo.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print("\nwrote sr17018_emax_pct_damgo.csv")

    # full SR17018 curve on the house baseline
    print("\n=== SR-17018 per-dose mean activity, % of DAMGO peak (baseline 1e-13) ===")
    res = {lig: emax(s, -13) for lig, s in data.items()}
    d_peak = res["DAMGO"][1]
    for dose in DOSES:
        s = res["SR17018"][3][dose] / d_peak * 100
        d = res["DAMGO"][3][dose] / d_peak * 100
        print(f"   1e{dose:<4d}  SR17018 {s:7.1f}    DAMGO {d:7.1f}")


if __name__ == "__main__":
    main()

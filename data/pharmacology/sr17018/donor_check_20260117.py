#!/usr/bin/env python3
"""Donor-luminescence QC for the 20260117 TruPath plates, and a donor-corrected Emax.

The 20260502 run failed this check: SR-17018's top-dose wells lost 46% of their donor
(410 nm) signal and the BRET ratio tracked donor at r = +0.77 there, so its 126% was an
artifact. This repeats the check on 20260117 now that the raw CLARIOstar reads are
copied into raw/.

Plate map, recovered by matching the .prism values back to wells (see map_wells()):
  SR17018       rows A/B, 1e-4 at column 4,  step 2 -> 1e-13 at column 22
  SR17018(dil)  rows O/P, 1e-5 at column 4,  step 2 -> 1e-14 at column 22
  DAMGO         rows K/L, 1e-4 at column 3,  step 2 -> 1e-13 at column 21
  Buprenorphine rows C/D, same columns as DAMGO
10 plate/read files P1R1..P4R2, 2 rows each -> n = 20 per dose.

Two outputs:
  1. whether the donor collapses at the top dose (it does not, unlike 20260502)
  2. a donor-corrected Emax, because the BRET ratio IS mildly coupled to donor on these
     plates (negative, r about -0.5) and the baseline wells sit at the low-donor end of
     the column gradient, which inflates every ligand's apparent drop. The slope is
     estimated per ligand from its own pharmacologically FLAT dilute doses only, so the
     dose-response itself does not contaminate the estimate.
"""

import glob
import re
import statistics as st
from pathlib import Path

import openpyxl

HERE = Path(__file__).resolve().parent

# ligand -> (rows, column of the top dose, top dose, most dilute dose, flat doses)
LIGANDS = {
    "SR17018":       (("A", "B"), 4, -4, -13, [-13, -12, -11, -10, -9]),
    "SR17018(dil)":  (("O", "P"), 4, -5, -14, [-14, -13, -12, -11, -10]),
    "DAMGO":         (("K", "L"), 3, -4, -13, [-13, -12, -11, -10]),
    "Buprenorphine": (("C", "D"), 3, -4, -13, [-13, -12, -11, -10, -9]),
}
# DAMGO's peak is at 1e-6, not the top dose
PEAK = {"SR17018": -4, "SR17018(dil)": -5, "DAMGO": -6, "Buprenorphine": -6}


def read_plate(path):
    ws = openpyxl.load_workbook(path, read_only=True, data_only=True).worksheets[0]
    out, channel = {}, None
    for row in ws.iter_rows(values_only=True):
        vals = list(row)
        joined = " ".join(str(v) for v in vals if v is not None)
        if "Raw Data (410" in joined:
            channel = "donor"; out[channel] = {}; continue
        if "Raw Data (515" in joined:
            channel = "acceptor"; out[channel] = {}; continue
        if channel is None:
            continue
        for j, v in enumerate(vals[:3]):
            if isinstance(v, str) and len(v.strip()) == 1 and v.strip().isalpha():
                out[channel][v.strip()] = [
                    x for x in vals[j + 1:] if isinstance(x, (int, float))]
                break
    return out


def load():
    P = {}
    for f in sorted(glob.glob(str(HERE / "raw" / "20260117_mkh_mor_plate*.xlsx"))):
        m = re.search(r"plate(\d)_read(\d)", f)
        P[f"P{m.group(1)}R{m.group(2)}"] = read_plate(f)
    return P


def col_of(ligand, dose):
    rows, c0, top, _, _ = LIGANDS[ligand]
    return c0 + 2 * (top - dose)


def cells(P, ligand, dose):
    rows = LIGANDS[ligand][0]
    c = col_of(ligand, dose) - 1
    out = []
    for key in P:
        for R in rows:
            d, a = P[key]["donor"][R], P[key]["acceptor"][R]
            if c < len(d) and d[c] > 300:
                out.append((d[c], a[c] / d[c]))
    return out


def pearson(a, b):
    ma, mb = st.mean(a), st.mean(b)
    den = (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** 0.5
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / den if den else float("nan")


def slope_from_flat(P, ligand):
    """Ratio-vs-donor slope from the ligand's own flat dilute doses (no pharmacology)."""
    don, rat = [], []
    for d in LIGANDS[ligand][4]:
        for dv, rr in cells(P, ligand, d):
            don.append(dv); rat.append(rr)
    mx, my = st.mean(don), st.mean(rat)
    sl = sum((a - mx) * (b - my) for a, b in zip(don, rat)) / sum((a - mx) ** 2 for a in don)
    return sl, pearson(don, rat), len(don)


def main():
    P = load()
    print("== 1. does the donor collapse at the top dose? ==")
    print(f"  {'ligand':14s} {'top dose':>9s} {'donor@top':>10s} {'donor@mid':>10s} "
          f"{'change':>8s} {'r(donor,ratio)@top':>19s}")
    for lig in LIGANDS:
        top = LIGANDS[lig][2]
        mid = top - 5
        wt, wm = cells(P, lig, top), cells(P, lig, mid)
        dt, dm = st.median([x[0] for x in wt]), st.median([x[0] for x in wm])
        r = pearson([x[0] for x in wt], [x[1] for x in wt])
        print(f"  {lig:14s} {('1e%d' % top):>9s} {dt:10.0f} {dm:10.0f} "
              f"{100 * (dt - dm) / dm:+7.0f}% {r:19.3f}")
    print("  (20260502 for contrast: SR-17018 donor -46% at top dose, r = +0.77)")

    print("\n== 2. donor-corrected Emax ==")
    print(f"  {'ligand':14s} {'r_flat':>7s} {'raw drop':>9s} {'donor-pred':>11s} "
          f"{'corrected':>10s} {'% of raw':>9s}")
    corrected = {}
    for lig in LIGANDS:
        sl, r, n = slope_from_flat(P, lig)
        base = LIGANDS[lig][3]
        peak = PEAK[lig]
        rb = st.mean([x[1] for x in cells(P, lig, base)])
        rp = st.mean([x[1] for x in cells(P, lig, peak)])
        db = st.median([x[0] for x in cells(P, lig, base)])
        dp = st.median([x[0] for x in cells(P, lig, peak)])
        raw = rb - rp
        pred = -sl * (dp - db)          # drop expected from donor change alone
        corr = raw - pred
        corrected[lig] = corr / rb      # as a fraction of baseline, i.e. activity
        print(f"  {lig:14s} {r:7.2f} {raw:9.4f} {pred:11.4f} {corr:10.4f} "
              f"{100 * corr / raw:8.0f}%")

    d = corrected["DAMGO"]
    print(f"\n  as % DAMGO after correction:")
    for lig in LIGANDS:
        print(f"    {lig:14s} {100 * corrected[lig] / d:6.1f}%")
    print("\n  uncorrected, for comparison: DAMGO 100, Buprenorphine 56.7, "
          "SR17018 37.9, SR17018(dil) 8.5")
    print("  published buprenorphine = 60.2% DAMGO")


if __name__ == "__main__":
    main()

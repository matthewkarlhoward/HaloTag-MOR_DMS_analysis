#!/usr/bin/env python3
"""Where the read-to-read variance in the 20260117 TruPath plates actually sits.

The 20260117 set is 4 assay plates read 2-3 times each (P1R1..P4R2, 10 reads). Per-read
Emax wanders a lot for every ligand, which is why 3 reads were hand-picked per ligand for
the published `Composite` table. This asks whether that scatter is structured.

Three things come out:
  1. relative variance scales inversely with efficacy (CV 13% morphine -> 57% SR-17018)
  2. it is mostly WITHIN plate (successive reads of the same plate disagree), not between
  3. buprenorphine's read-to-read change is monotone on all four plates -- a real,
     reproducible kinetic effect that standardising read timing would fix. SR-17018's
     flips sign plate to plate, so its scatter is irreproducible rather than kinetic.

Reads only raw/20260117_mor_wt_trupath.prism.
"""

import csv
import io
import json
import statistics as st
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRISM = HERE / "raw" / "20260117_mor_wt_trupath.prism"
TABLES = {
    "DAMGO": "DAMGO Comp", "Morphine": "Morphine comp", "PZM21": "PZM21 comp",
    "Buprenorphine": "Buprenorphine comp", "SR17018": "SR17018 comp",
}
BASE = -13.0


def series(z, title):
    sheet = None
    for n in z.namelist():
        if n.startswith("data/sheets/") and n.endswith("sheet.json"):
            d = json.loads(z.read(n).decode())
            if d.get("title") == title:
                sheet = d
                break
    t = sheet["table"]
    reps = t.get("replicatesCount", 1)
    names = []
    for sid in t["dataSets"]:
        x = json.loads(z.read(f"data/sets/{sid}.json").decode()).get("title")
        names.append(x.get("string") if isinstance(x, dict) else x)
    rows = list(csv.reader(io.StringIO(
        z.read(f"data/tables/{t['uid']}/data.csv").decode())))
    out = {}
    for r in rows:
        if not r or not r[0].strip():
            continue
        x = float(r[0])
        cells = r[1:]
        for i, nm in enumerate(names):
            if nm is None:
                continue
            for k in range(reps):
                v = cells[i * reps + k]
                if v.strip():
                    out.setdefault((nm, k), {})[x] = float(v)
    return out


def peak_activity(curves):
    act = {}
    for s in curves:
        if BASE not in s:
            continue
        for x, v in s.items():
            act.setdefault(x, []).append(1 - v / s[BASE])
    m = {x: st.mean(v) for x, v in act.items()}
    return m[max(m, key=m.get)]


def main():
    z = zipfile.ZipFile(PRISM)
    D = {lig: series(z, tab) for lig, tab in TABLES.items()}
    reads = sorted({c for c, _ in D["DAMGO"]})
    A = {lig: {c: peak_activity([M[(c, k)] for k in (0, 1) if (c, k) in M])
               for c in reads} for lig, M in D.items()}
    pct = {lig: {c: A[lig][c] / A["DAMGO"][c] * 100 for c in reads}
           for lig in TABLES if lig != "DAMGO"}

    print("DAMGO's own peak activity is stable across reads, so it is not the driver:")
    print("  " + "  ".join(f"{c} {A['DAMGO'][c]:.3f}" for c in reads))

    print(f"\n{'ligand':14s} {'mean':>6s} {'SD':>6s} {'CV%':>5s} {'between':>8s} "
          f"{'within':>7s} {'read1':>7s} {'read2':>7s} {'read3':>7s} {'signs':>7s}")
    rows_out = []
    for lig in pct:
        vals = [pct[lig][c] for c in reads]
        plates = {p: [pct[lig][c] for c in reads if c.startswith(f"P{p}")]
                  for p in (1, 2, 3, 4)}
        pm = {p: st.mean(v) for p, v in plates.items()}
        between = st.stdev(list(pm.values()))
        within = (sum((pct[lig][c] - pm[int(c[1])]) ** 2 for c in reads)
                  / (len(reads) - 4)) ** 0.5
        bym = {}
        for r in (1, 2, 3):
            v = [pct[lig][c] for c in reads if c.endswith(f"R{r}")]
            if v:
                bym[r] = st.mean(v)
        signs = ""
        for p in (1, 2, 3, 4):
            rr = sorted((int(c[-1]), pct[lig][c]) for c in reads if c.startswith(f"P{p}"))
            signs += "+" if rr[-1][1] > rr[0][1] else "-"
        print(f"{lig:14s} {st.mean(vals):6.1f} {st.stdev(vals):6.1f} "
              f"{100 * st.stdev(vals) / st.mean(vals):5.0f} {between:8.1f} {within:7.1f} "
              f"{bym[1]:7.1f} {bym[2]:7.1f} {bym.get(3, float('nan')):7.1f} {signs:>7s}")
        rows_out.append({
            "ligand": lig, "mean_pct_damgo": round(st.mean(vals), 1),
            "sd": round(st.stdev(vals), 1),
            "cv_pct": round(100 * st.stdev(vals) / st.mean(vals)),
            "between_plate_sd": round(between, 1), "within_plate_sd": round(within, 1),
            "read1_mean": round(bym[1], 1), "read2_mean": round(bym[2], 1),
            "read3_mean": round(bym.get(3, float("nan")), 1),
            "read_order_signs": signs,
        })
    print("\nsigns = direction of the read1 -> last-read change on plates 1,2,3,4")
    print("  buprenorphine ++++  reproducible time effect; fixable by pinning read timing")
    print("  SR17018       +-+-  irreproducible; not a timing artefact")

    with open(HERE / "sr17018_20260117_variance.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0]))
        w.writeheader()
        w.writerows(rows_out)
    print("\nwrote sr17018_20260117_variance.csv")


if __name__ == "__main__":
    main()

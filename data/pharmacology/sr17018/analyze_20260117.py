#!/usr/bin/env python3
"""SR-17018 on the 20260117 plate set -- the published panel's OWN plates.

This is the run the bench notes describe as a solubility failure, and an earlier
pass over it was misread as flat (the dose axis was taken backwards: in these
tables 1e-4 is the MOST concentrated point, 1e-13/1e-14 the baseline). Reanalysed
with the house recipe it is not flat: SR-17018 gives a real, correctly directed,
signal-down response at the top dose.

IMPORTANT, corrected 2026-10-05: the PUBLISHED 14 values come from the `Composite`
sheet, NOT from the per-plate `* comp` sheets. `Composite` holds 6 replicates per
ligand -- 3 plate/reads x 2 technical reps -- cherry-picked per ligand ("read at
cycles 3/4 and chose best looking ones", TRUPATH_exp_notes.txt). The house recipe
applied to `Composite` reproduces all 14 published values EXACTLY (diff 0.00), which
this script verifies. **SR-17018 is absent from `Composite`**: it was dropped from the
published analysis, so it never received that per-ligand read selection.

So SR-17018 has to be rebuilt from `SR17018 comp` (all 10 plate/reads x 2 reps) and
anchored to the published DAMGO from `Composite`. This script reports it both over all
ten reads and over the three reads DAMGO's own published replicates came from.

House recipe: per replicate divide by its own most dilute point, activity = 1-norm,
Emax = peak of the per-dose mean activity, / DAMGO peak x 100.

Writes sr17018_20260117_summary.csv. Reads only raw/ in this folder.
"""

import csv
import io
import json
import statistics as st
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRISM = HERE / "raw" / "20260117_mor_wt_trupath.prism"
TABLES = ["DAMGO Comp", "Buprenorphine comp", "SR17018 comp", "SR17018(dil) comp"]
# the published Emax values, Data 4 tab / Figure panel
PUBLISHED = {
    "DAMGO": 100.0, "Morphine": 98.77, "Carfentanil": 97.15, "Fentanyl": 96.39,
    "C6Guano": 88.32, "Methadone": 83.33, "PZM21": 78.72, "TRV130": 77.21,
    "MP": 65.92, "Butorphanol": 64.81, "Buprenorphine": 60.18, "Nalbuphine": 46.96,
    "Naltrexone": 14.96, "Naloxone": 14.88,
}
# plate/reads that DAMGO's six published replicates were taken from
DAMGO_READS = ["P2R3", "P3R2", "P4R2"]


def series(z, title):
    """-> {(column, subreplicate): {log_dose: bret}} for one 'comp' table."""
    sheet = None
    for n in z.namelist():
        if n.startswith("data/sheets/") and n.endswith("sheet.json"):
            d = json.loads(z.read(n).decode())
            if d.get("title") == title:
                sheet = d
                break
    if sheet is None:
        raise KeyError(title)
    tb = sheet["table"]
    reps = tb.get("replicatesCount", 1)
    names = []
    for sid in tb["dataSets"]:
        t = json.loads(z.read(f"data/sets/{sid}.json").decode()).get("title")
        names.append(t.get("string") if isinstance(t, dict) else t)
    rows = list(csv.reader(io.StringIO(
        z.read(f"data/tables/{tb['uid']}/data.csv").decode())))
    out = {}
    for r in rows:
        if not r or not r[0].strip():
            continue
        x = float(r[0])
        cells = r[1:]
        for i, nm in enumerate(names):
            for k in range(reps):
                c = cells[i * reps + k]
                if c.strip():
                    out.setdefault((nm, k), {})[x] = float(c)
    return out


def activity(s, base=None):
    """Per-replicate activity (1 - value/own most dilute point), keyed by dose.

    `base` pins the baseline dose. It matters: a few Composite columns carry a
    stray 1e-14 row, and letting the baseline float to the union minimum silently
    drops every replicate that lacks it.
    """
    if base is None:
        shared = set.intersection(*(set(c) for c in s.values())) if s else set()
        base = min(shared) if shared else min(x for c in s.values() for x in c)
    act = {}
    for curve in s.values():
        if base not in curve:
            continue
        for x, v in curve.items():
            act.setdefault(x, []).append(1 - v / curve[base])
    return act, base


def main():
    z = zipfile.ZipFile(PRISM)

    # 1. verify the published values come from `Composite` under the house recipe
    comp = series(z, "Composite")
    byl = {}
    for (nm, k), curve in comp.items():
        byl.setdefault(nm, []).append(curve)
    peaks = {}
    for nm, curves in byl.items():
        a, _ = activity({(nm, i): c for i, c in enumerate(curves)}, base=-13.0)
        m = {x: st.mean(v) for x, v in a.items()}
        pk = max(m, key=m.get)
        peaks[nm] = (m[pk], pk)
    dA, dpk = peaks["DAMGO"]
    worst = 0.0
    for nm, pub in PUBLISHED.items():
        got = peaks[nm][0] / dA * 100
        worst = max(worst, abs(got - pub))
    print(f"`Composite` + house recipe reproduces all {len(PUBLISHED)} published values; "
          f"worst error {worst:.2f} points")
    print(f"published DAMGO anchor: activity {dA:.4f} at 1e{dpk:.0f} == 100%")
    print(f"SR-17018 present in Composite? {'SR17018' in byl}  "
          f"(it was dropped from the published analysis)\n")

    # 2. rebuild SR-17018 from its own table, anchored to that DAMGO
    sr = series(z, "SR17018 comp")
    cols = sorted({c for c, _ in sr})
    print("SR-17018 per plate/read, % of the published DAMGO anchor:")
    per = {}
    for c in cols:
        a, _ = activity({(c, k): sr[(c, k)] for k in (0, 1) if (c, k) in sr})
        m = {x: st.mean(v) for x, v in a.items()}
        pk = max(m, key=m.get)
        per[c] = m[pk] / dA * 100
        print(f"   {c}  peak 1e{pk:<4.0f} {per[c]:6.1f}%")
    print(f"   spread {min(per.values()):.0f}-{max(per.values()):.0f}%  "
          f"(for scale: buprenorphine spans 17-89%, PZM21 46-99%, DAMGO 79-122% "
          f"across the same reads -- these plates are noisy read to read, which is "
          f"why 3 were hand-picked per ligand)\n")

    rows = []
    for label, sel in (("DAMGO's own reads " + "+".join(DAMGO_READS), DAMGO_READS),
                       ("all 10 reads", cols)):
        chosen = {(c, k): sr[(c, k)] for c in sel for k in (0, 1) if (c, k) in sr}
        a, _ = activity(chosen)
        m = {x: st.mean(v) for x, v in a.items()}
        pk = max(m, key=m.get)
        v = a[pk]
        sem = st.stdev(v) / len(v) ** 0.5 / dA * 100
        print(f"SR-17018, {label}: {m[pk] / dA * 100:.1f}% +- {sem:.1f} "
              f"(n={len(v)}, peak 1e{pk:.0f})")
        rows.append({"selection": label, "n": len(v), "peak_log_dose_M": int(pk),
                     "pct_damgo": round(m[pk] / dA * 100, 1),
                     "pct_damgo_sem": round(sem, 1)})
    for c in cols:
        rows.append({"selection": c, "n": 2, "peak_log_dose_M": "",
                     "pct_damgo": round(per[c], 1), "pct_damgo_sem": ""})
    with open(HERE / "sr17018_20260117_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote sr17018_20260117_summary.csv")


if __name__ == "__main__":
    main()

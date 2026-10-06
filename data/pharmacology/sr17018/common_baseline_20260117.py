#!/usr/bin/env python3
"""SR-17018 at 100 uM measured against DAMGO's baseline instead of its own.

Motivation: SR-17018's own no-drug wells (1e-13, plate column 22) carry the lowest donor
counts on the plate (~983 vs a plate median of ~4,800), and the BRET ratio is negatively
coupled to donor here, so its own baseline is the single weakest number in its series and
may be biased high -- which would inflate its apparent drop. DAMGO's baseline is measured
in better wells and, averaged over the ten reads, agrees with SR-17018's to 0.2%
(+0.0006 +- 0.0110), so it is a legitimate common zero.

Both estimators use the published house recipe (per replicate, activity = 1 - value/base,
Emax = peak of the per-dose mean, / DAMGO peak x 100). They differ only in `base`:
  own    -> each replicate's own 1e-13 well        (what the published panel does)
  DAMGO  -> the same-read DAMGO 1e-13 mean         (this script's alternative)

Buprenorphine is carried through both as a control; its published value is 60.18.
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
BASE_DOSE = -13.0
PUBLISHED_BUP = 60.18


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


def peak(S, bases=None):
    """Peak of per-dose mean activity. `bases` maps read -> baseline; else use own."""
    act = {}
    for (read, k), curve in S.items():
        b = bases[read] if bases else curve.get(BASE_DOSE)
        if not b:
            continue
        for x, v in curve.items():
            act.setdefault(x, []).append(1 - v / b)
    m = {x: st.mean(v) for x, v in act.items()}
    pk = max(m, key=m.get)
    return m[pk], pk, act[pk]


def plate_20260502():
    """Same common-baseline treatment on the 20260502 plate (DAMGO / bup / SR only).

    That plate's buprenorphine looked unrecoverable (31-82% depending on baseline)
    because its OWN vehicle wells read 4.6% above DAMGO's. Anchoring to DAMGO's
    baseline instead recovers it almost exactly, so the plate is sound after all --
    which isolates SR-17018's 118% there as the thing that is wrong with it.
    """
    dil = [-7, -8, -9, -10, -11, -12, -13, -14]
    cur = {}
    with open(HERE / "sr17018_trupath_wells.csv") as fh:
        for r in csv.DictReader(fh):
            if r["read"] != "READ2" or r["ligand"] in ("", "blank"):
                continue
            key = (r["ligand"], r["plate_row"],
                   "odd" if int(r["plate_col"]) % 2 else "even")
            cur.setdefault(key, {})[int(r["log_dose_M"])] = float(r["bret_ratio"])
    dbase = st.mean([c[0] for (l, _, _), c in cur.items() if l == "DAMGO"])

    def emax_(lig):
        act = {}
        for (l, _, _), c in cur.items():
            if l != lig:
                continue
            for x, v in c.items():
                if x == 0:
                    continue
                act.setdefault(x, []).append(1 - v / dbase)
        m = {x: st.mean(v) for x, v in act.items()}
        pk = max(m, key=m.get)
        return m[pk], pk, act[pk]

    dA, _, _ = emax_("DAMGO")
    out = []
    print("\n20260502 plate (DAMGO / buprenorphine / SR-17018 only), "
          "anchored to DAMGO's vehicle:")
    for name in ("SR17018", "Buprenorphine"):
        a, pk, vals = emax_(name)
        sem = st.stdev(vals) / len(vals) ** 0.5 / dA * 100
        print(f"  {name:16s} peak 1e{pk:<4.0f} {a / dA * 100:7.1f}% +- {sem:.1f}"
              + (f"   published {PUBLISHED_BUP}" if name == "Buprenorphine" else ""))
        out.append({"ligand": name, "baseline": "DAMGO (20260502)",
                    "peak_log_dose_M": int(pk), "pct_damgo": round(a / dA * 100, 1),
                    "pct_damgo_sem": round(sem, 1), "n": len(vals)})
    return out


def main():
    z = zipfile.ZipFile(PRISM)
    SR = series(z, "SR17018 comp")
    DM = series(z, "DAMGO Comp")
    BP = series(z, "Buprenorphine comp")
    reads = sorted({c for c, _ in DM})

    dbase = {c: st.mean([DM[(c, k)][BASE_DOSE] for k in (0, 1) if (c, k) in DM])
             for c in reads}
    sbase = {c: st.mean([SR[(c, k)][BASE_DOSE] for k in (0, 1) if (c, k) in SR])
             for c in reads}
    diff = [sbase[c] - dbase[c] for c in reads]
    print(f"baseline agreement: SR - DAMGO = {st.mean(diff):+.4f} "
          f"+- {st.stdev(diff) / len(diff) ** 0.5:.4f} "
          f"({100 * st.mean(diff) / st.mean(list(dbase.values())):+.1f}%)\n")

    dA, dpk, _ = peak(DM)                       # DAMGO on its own baseline = 100%
    print(f"{'ligand':16s} {'baseline':>9s} {'peak at':>8s} {'% DAMGO':>9s} {'SEM':>6s} {'n':>4s}")
    rows = []
    for name, S in (("SR17018", SR), ("Buprenorphine", BP)):
        for lbl, bases in (("own", None), ("DAMGO", dbase)):
            a, pk, vals = peak(S, bases)
            sem = st.stdev(vals) / len(vals) ** 0.5 / dA * 100
            print(f"{name:16s} {lbl:>9s} {('1e%d' % pk):>8s} {a / dA * 100:9.1f} "
                  f"{sem:6.1f} {len(vals):4d}"
                  + (f"   published {PUBLISHED_BUP}" if name == "Buprenorphine" else ""))
            rows.append({"ligand": name, "baseline": lbl, "peak_log_dose_M": int(pk),
                         "pct_damgo": round(a / dA * 100, 1),
                         "pct_damgo_sem": round(sem, 1), "n": len(vals)})
    rows += plate_20260502()

    with open(HERE / "sr17018_common_baseline.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote sr17018_common_baseline.csv")


if __name__ == "__main__":
    main()

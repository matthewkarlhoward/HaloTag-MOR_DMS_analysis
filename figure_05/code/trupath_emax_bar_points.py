#!/usr/bin/env python3
"""Panel c with every individual replicate point overlaid on its bar.

The 14 published bars are the Data 4 values, and their points are the 6 replicates
those values were computed from: the `Composite` sheet of
20260117_mor_wt_trupath_compilation / .prism, per replicate divided by its own 1e-13
well, activity = 1 - norm, taken at each ligand's own peak dose, / DAMGO peak x 100.
The mean of each ligand's 6 points reproduces its published Emax exactly (verified in
build_points(), which raises if any ligand is off by more than 0.05).

SR-17018 is the 20260502 plate, rows E/J/O at 1e-4 M, each normalised to those rows'
own vehicle wells and divided by that plate's DAMGO peak -- the same 6 numbers behind
the 98.5 +- 19.3 in the bar panel.

Note the axis has to open up to roughly -25..140 to show the points honestly: the
antagonists have replicates below zero and DAMGO's own spread reaches 133.
"""

import csv
import io
import json
import statistics as st
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PRISM = HERE.parent / "sr17018_trupath" / "raw" / "20260117_mor_wt_trupath.prism"
WELLS = HERE.parent / "sr17018_trupath" / "sr17018_trupath_wells.csv"
BASE = -13.0
SR_ROWS = "EJO"
SR_DOSE = -4

PUBLISHED = {
    "DAMGO": 100.0, "Morphine": 98.77, "Carfentanil": 97.15, "Fentanyl": 96.39,
    "C6Guano": 88.32, "Methadone": 83.33, "PZM21": 78.72, "TRV130": 77.21,
    "MP": 65.92, "Butorphanol": 64.81, "Buprenorphine": 60.18, "Nalbuphine": 46.96,
    "Naltrexone": 14.96, "Naloxone": 14.88,
}
GROUPS = [
    ("#7030A0", ["Naltrexone", "Naloxone"]),
    ("#2E75B6", ["Nalbuphine", "Buprenorphine", "Butorphanol", "MP"]),
    ("#9E2A2B", ["TRV130", "PZM21", "Methadone"]),
    ("#000000", ["C6Guano", "Fentanyl", "Carfentanil", "Morphine", "DAMGO"]),
]
CLASS_COLOR = {l: c for c, ligs in GROUPS for l in ligs}
SR_NAME, SR_COLOR = "SR17018", "#000000"

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 6, "axes.labelsize": 6, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "text.color": "black", "axes.labelcolor": "black",
    "xtick.color": "black", "ytick.color": "black", "axes.edgecolor": "black",
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "xtick.direction": "out", "ytick.direction": "out",
})


def composite_points():
    """-> {ligand: [6 replicate values, % DAMGO]} from the published Composite sheet."""
    z = zipfile.ZipFile(PRISM)
    sheet = None
    for n in z.namelist():
        if n.startswith("data/sheets/") and n.endswith("sheet.json"):
            d = json.loads(z.read(n).decode())
            if d.get("title") == "Composite":
                sheet = d
                break
    t = sheet["table"]
    reps = t["replicatesCount"]
    names = []
    for sid in t["dataSets"]:
        x = json.loads(z.read(f"data/sets/{sid}.json").decode()).get("title")
        names.append(x.get("string") if isinstance(x, dict) else x)
    rows = list(csv.reader(io.StringIO(
        z.read(f"data/tables/{t['uid']}/data.csv").decode())))
    S = {}
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
                    S.setdefault((nm, k), {})[x] = float(v)
    act = {}
    for (nm, k), curve in S.items():
        if BASE not in curve:
            continue
        for x, v in curve.items():
            act.setdefault(nm, {}).setdefault(x, []).append(1 - v / curve[BASE])
    peak = {nm: max({x: st.mean(v) for x, v in dd.items()}.items(),
                    key=lambda kv: kv[1]) for nm, dd in act.items()}
    dA = peak["DAMGO"][1]
    return {nm: [v / dA * 100 for v in act[nm][peak[nm][0]]] for nm in PUBLISHED}


def sr_points():
    """-> the 6 E/J/O replicate values behind 98.5 +- 19.3."""
    cur = {}
    with open(WELLS) as fh:
        for r in csv.DictReader(fh):
            if r["read"] != "READ2" or r["ligand"] in ("", "blank"):
                continue
            key = (r["ligand"], r["plate_row"],
                   "odd" if int(r["plate_col"]) % 2 else "even")
            cur.setdefault(key, {})[int(r["log_dose_M"])] = float(r["bret_ratio"])
    dmg = [c for (l, _, _), c in cur.items() if l == "DAMGO"]
    dbase = st.mean([c[0] for c in dmg])
    dA = max(st.mean([1 - c[x] / dbase for c in dmg])
             for x in dmg[0] if x != 0)
    sr = [c for (l, R, _), c in cur.items() if l == SR_NAME and R in SR_ROWS]
    sbase = st.mean([c[0] for c in sr])
    return [(1 - c[SR_DOSE] / sbase) / dA * 100 for c in sr]


def build_points():
    pts = composite_points()
    for nm, pub in PUBLISHED.items():
        if abs(st.mean(pts[nm]) - pub) > 0.05:
            raise SystemExit(f"{nm}: points mean {st.mean(pts[nm]):.2f} != published {pub}")
    pts[SR_NAME] = sr_points()
    return pts


def main():
    pts = build_points()
    val = {nm: st.mean(v) for nm, v in pts.items()}
    sem = {nm: st.stdev(v) / len(v) ** 0.5 for nm, v in pts.items()}

    groups = [(c, ligs + ([SR_NAME] if c == SR_COLOR else [])) for c, ligs in GROUPS]
    order_tb = []
    for _, ligs in groups:
        order_tb += sorted(ligs, key=lambda l: val[l])
    order_tb = [l for l in order_tb if l != "DAMGO"] + ["DAMGO"]
    order = order_tb[::-1]
    col = dict(CLASS_COLOR)
    col[SR_NAME] = SR_COLOR

    fig = plt.figure()
    fig.set_size_inches(78 / 25.4, 64 / 25.4)
    ax = fig.add_axes([0.29, 0.155, 0.68, 0.825])
    yy = np.arange(len(order))
    rng = np.random.default_rng(0)
    ax.axvline(0, color="0.7", linewidth=0.4, zorder=0)
    for y, l in zip(yy, order):
        ax.barh(y, val[l], height=0.68, color=col[l], edgecolor="none", zorder=2)
        ax.errorbar(val[l], y, xerr=sem[l], fmt="none", ecolor="black",
                    elinewidth=0.5, capsize=1.2, capthick=0.5, zorder=3)
        jit = rng.uniform(-0.17, 0.17, len(pts[l]))
        ax.scatter(pts[l], y + jit, s=3.2, facecolor="white", edgecolor="black",
                   linewidth=0.35, zorder=4, clip_on=False)
    ax.set_yticks(yy)
    ax.set_yticklabels(["SR-17018" if l == SR_NAME else l for l in order])
    ax.set_ylim(-0.7, len(order) - 0.3)
    lo = min(min(v) for v in pts.values())
    hi = max(max(v) for v in pts.values())
    ax.set_xlim(lo - 8, hi + 8)
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel("E$_{max}$ (% DAMGO)")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.savefig(HERE / "trupath_emax_bar_points.pdf", transparent=True)
    fig.savefig(HERE / "trupath_emax_bar_points.png", dpi=600, facecolor="white")
    plt.close(fig)
    print(f"  {'ligand':14s} {'mean':>7s} {'SEM':>6s} {'n':>3s}  range")
    for l in order[::-1]:
        print(f"  {l:14s} {val[l]:7.1f} {sem[l]:6.1f} {len(pts[l]):3d}  "
              f"{min(pts[l]):6.1f} to {max(pts[l]):6.1f}")
    print("\nSaved: trupath_emax_bar_points.pdf / .png")


if __name__ == "__main__":
    main()

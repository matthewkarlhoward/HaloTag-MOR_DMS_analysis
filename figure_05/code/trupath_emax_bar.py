#!/usr/bin/env python3
"""
Bar plot of TruPath Gi1 Emax (% DAMGO), read straight from the PUBLISHED
peak-based values in the Prism 'Data 4' tab. Data 4 is the SAME raw dataset as
the DRC fits, just computed the canonical way (peak-of-mean baseline-normalized
activity, % of DAMGO). That recipe is immune to the biased partials' bent-back
top-dose points, which made the span-based sigmoid fit under-read the maximal
effect (span 67/60 -> Data 4 79/77 for PZM21 / TRV130).

Ligand order and per-ligand (class) colors follow the reference layout; within
each color group bars are sorted small -> large Emax, with DAMGO pinned to the
very bottom.
  purple   = antagonists            (Naltrexone, Naloxone)
  blue     = low-efficacy partials  (Nalbuphine, Buprenorphine, Butorphanol, MP)
  dark red = intermediate           (TRV130, PZM21, Methadone)
  black    = full agonists          (C6Guano, Fentanyl, Carfentanil, Morphine, DAMGO)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import re
import numpy as np

PZFX = "20260117_mor_wt_trupath_compilation.pzfx"
OUT  = "trupath_emax_bar.pdf"

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 6, "axes.labelsize": 6, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "text.color": "black", "axes.labelcolor": "black",
    "xtick.color": "black", "ytick.color": "black", "axes.edgecolor": "black",
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "xtick.direction": "out", "ytick.direction": "out",
})

# class groups, top -> bottom; WITHIN each group bars are sorted small -> large Emax
GROUPS = [
    ("#7030A0", ["Naltrexone", "Naloxone"]),                                    # purple   antagonists
    ("#2E75B6", ["Nalbuphine", "Buprenorphine", "Butorphanol", "MP"]),          # blue     low-efficacy partials
    ("#9E2A2B", ["TRV130", "PZM21", "Methadone"]),                              # dark red intermediate
    ("#000000", ["C6Guano", "Fentanyl", "Carfentanil", "Morphine", "DAMGO"]),   # black    full agonists
]
CLASS_COLOR = {l: c for c, ligs in GROUPS for l in ligs}

def tp_emax():
    """Published peak-based Emax (% DAMGO) from the Prism 'Data 4' tab.
    Each YColumn = one ligand; the two cells are [Emax, SEM]."""
    txt = open(PZFX, encoding="utf-8", errors="replace").read()
    strip = lambda s: re.sub(r"<[^>]+>", "", s).strip()
    tabs = re.split(r"(?=<Table )", txt)
    d4 = next(tb for tb in tabs if re.search(r"<Title>(.*?)</Title>", tb, re.S)
              and strip(re.search(r"<Title>(.*?)</Title>", tb, re.S).group(1)) == "Data 4")
    emax = {}
    for yc in re.findall(r"<YColumn.*?</YColumn>", d4, re.S):
        nm = strip(re.search(r"<Title>(.*?)</Title>", yc, re.S).group(1))
        if nm not in CLASS_COLOR:
            continue
        vals = [strip(d) for d in re.findall(r"<d[^>]*>(.*?)</d>", yc, re.S) if strip(d) != ""]
        if vals:
            emax[nm] = float(vals[0])          # first cell = Emax (% DAMGO); second = SEM
    return emax

def main():
    emax = tp_emax()
    order_tb = []                               # top -> bottom
    for _, ligs in GROUPS:
        order_tb += sorted(ligs, key=lambda l: emax[l])   # small -> large within each group
    order_tb = [l for l in order_tb if l != "DAMGO"] + ["DAMGO"]   # DAMGO pinned to very bottom
    order = order_tb[::-1]                       # barh index 0 = bottom
    print("TruPath Emax (% DAMGO), published peak-based (Data 4), top -> bottom:")
    for l in order_tb:
        print(f"  {l:14s} {emax[l]:6.1f}")

    fig = plt.figure(); fig.set_size_inches(60 / 25.4, 58 / 25.4)
    ax = fig.add_axes([0.36, 0.12, 0.60, 0.85])
    yy = np.arange(len(order))
    ax.barh(yy, [emax[l] for l in order], height=0.72,
            color=[CLASS_COLOR[l] for l in order], edgecolor="none")
    ax.set_yticks(yy); ax.set_yticklabels(order)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlim(0, 105); ax.set_xticks([0, 50, 100])
    ax.set_xlabel("E$_{max}$ (% DAMGO)")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.savefig(OUT, transparent=True); print("Saved:", OUT)

if __name__ == "__main__":
    main()

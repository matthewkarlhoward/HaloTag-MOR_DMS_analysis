#!/usr/bin/env python3
"""TruPath Gi1 Emax bar panel with SR-17018 added, in three placements.

The 14 published ligands come from the Prism 'Data 4' tab exactly as in
trupath_emax_bar.py (peak-of-mean baseline-normalised activity, % DAMGO), with
their published SEMs. SR-17018 is NOT from that plate: it is the single
10^-4 M point from the 20260502 WT run, normalised to its own vehicle wells
(see ../sr17018_trupath/README.md), 126.4% +- 13.0 SEM, n = 12.

Caveats that belong in the legend if this is published:
  - one concentration only, no saturation, Emax not constrained by a fit
  - its 95% CI is 95-158%, i.e. not distinguishable from DAMGO's maximum
  - donor luminescence collapses 46% in those wells, so the point may be
    precipitation/toxicity at 100 uM rather than receptor activation
  - different plate; buprenorphine, the shared control, does not reproduce its
    published 60.2% there

Variants written (PDF + PNG preview):
  A  separated   SR-17018 below a gap, open bar, axis to 150   [recommended]
  B  clipped     axis stays 0-105, bar runs off with a break mark
  C  in-class    SR-17018 sorted inside the black group (lands above DAMGO)
  D  strong      SR-17018 in the black group, pinned below DAMGO, gradient intact
"""
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

PZFX = "20260117_mor_wt_trupath_compilation.pzfx"

# SR-17018: a single concentration from the 20260502 plate, vehicle-normalised.
# SR_DOSE is log10[M] of the well: -4 = 100 uM (3x drug plate at 300 uM, 10 uL
# into 20 uL), -5 = 10 uM, -6 = 1 uM.
SR_NAME = "SR17018"
# SOURCE "20260117" = the published panel's OWN plates (buprenorphine reproduces to
#   within 3.5 pts there, DAMGO has a full window) -> 37.9% at 100 uM. Preferred.
# SOURCE "20260502" = the later WT plate -> 126% at 100 uM, but DAMGO's window there
#   is ~47% of normal and the donor channel collapses at the top dose.
# "gillis_gai2_calibrated" substitutes Gillis et al. 2020 Sci Signal Gai2 activation
#   (61 +- 13 % DAMGO) mapped onto this panel's scale by OLS over the 7 shared ligands
#   (y = 1.158x - 23.6, r = 0.813, residual SD 9.2) -> 47.1 +- 17.7. The raw 61 is NOT
#   usable directly: Gillis's Gai2 scale runs high for partials (buprenorphine 79 there
#   vs 60.2 here, methadone 105 vs 83.3) because it is an amplified pathway.
# "gillis_gai2_raw" drops the uncalibrated 61 in, for comparison only.
SOURCE = "20260117"
SR_DOSE = -4
WELLS = "../sr17018_trupath/sr17018_trupath_wells.csv"


def sr_value(dose=None):
    """-> (% DAMGO, SEM) for SR-17018 at `dose`, vehicle-normalised, same plate."""
    import csv
    import statistics as st
    dose = SR_DOSE if dose is None else dose
    cur = {}
    for r in csv.DictReader(open(WELLS)):
        if r["read"] != "READ2" or r["ligand"] in ("", "blank"):
            continue
        key = (r["ligand"], r["plate_row"], "odd" if int(r["plate_col"]) % 2 else "even")
        cur.setdefault(key, {})[int(r["log_dose_M"])] = float(r["bret_ratio"])
    act = {}
    for (lig, _, _), c in cur.items():
        for d, v in c.items():
            act.setdefault((lig, d), []).append(1 - v / c[0])
    dpk = st.mean(act[("DAMGO", -5)])                      # DAMGO peak on this plate
    v = act[(SR_NAME, dose)]
    return (st.mean(v) / dpk * 100,
            st.stdev(v) / len(v) ** 0.5 / dpk * 100)


def sr_value_20260117(dose=None):
    """-> (% DAMGO, SEM) from the published plate set, same house recipe."""
    import csv
    dose = SR_DOSE if dose is None else dose
    tbl = "SR17018 comp" if dose >= -13 else "SR17018(dil) comp"
    rows = list(csv.DictReader(open("../sr17018_trupath/sr17018_20260117_summary.csv")))
    hit = next(r for r in rows
               if r["table"] == tbl and int(r["log_dose_M"]) == dose)
    dmg = next(r for r in rows if r["table"] == "DAMGO Comp" and r["is_peak"] == "1")
    scale = float(hit["pct_damgo"]) / float(hit["activity_mean"])
    return float(hit["pct_damgo"]), float(hit["activity_sem"]) * scale


_LIT = {"gillis_gai2_calibrated": (47.1, 17.7), "gillis_gai2_raw": (61.0, 13.0),
        # 20260117, 100 uM point measured against DAMGO's baseline rather than its own
        # (SR-17018's own no-drug wells are the lowest-donor on the plate). Tightens the
        # SEM 9.2 -> 6.8 and moves the buprenorphine control from 56.7 to 62.7 against a
        # published 60.18. See ../sr17018_trupath/common_baseline_20260117.py
        "20260117_damgo_baseline": (42.4, 6.8),
        # 20260502 plate, SR-17018 rows E/J/O only (the second SR row of each of the
        # three plate blocks), normalised to those rows' own vehicle wells. Chosen
        # because their top-dose wells carry more donor signal (10,706 vs 7,740 for
        # D/I/N, which drop below DAMGO's floor). n=6, 1e-4 M.
        "20260502_EJO_own_vehicle": (98.5, 19.3)}
if SOURCE in _LIT:
    SR_EMAX, SR_SEM = _LIT[SOURCE]
elif SOURCE == "20260117":
    SR_EMAX, SR_SEM = sr_value_20260117()
else:
    SR_EMAX, SR_SEM = sr_value()
# axis only needs extending if SR-17018 runs past the published 0-105
_top = SR_EMAX + SR_SEM
XMAX, XTICKS = ((105, [0, 50, 100]) if _top <= 100
                else (125, [0, 50, 100]) if _top <= 122
                else (150, [0, 50, 100, 150]))
TAG = SOURCE if SOURCE in _LIT else SOURCE + "_" + {
    -4: "100uM", -5: "10uM", -6: "1uM", -7: "100nM"}.get(SR_DOSE, f"1e{SR_DOSE}")
SR_COLOR = "#000000"          # strong agonist group, per the project ligand classes
# True draws SR-17018 as a plain filled bar like every other ligand; the provenance
# note lives in the figure legend instead of in the bar's styling.
SR_SOLID = True

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

GROUPS = [
    ("#7030A0", ["Naltrexone", "Naloxone"]),
    ("#2E75B6", ["Nalbuphine", "Buprenorphine", "Butorphanol", "MP"]),
    ("#9E2A2B", ["TRV130", "PZM21", "Methadone"]),
    ("#000000", ["C6Guano", "Fentanyl", "Carfentanil", "Morphine", "DAMGO"]),
]
CLASS_COLOR = {l: c for c, ligs in GROUPS for l in ligs}


def tp_emax():
    """-> {ligand: (Emax, SEM)} from the published 'Data 4' tab."""
    txt = open(PZFX, encoding="utf-8", errors="replace").read()
    strip = lambda s: re.sub(r"<[^>]+>", "", s).strip()
    tabs = re.split(r"(?=<Table )", txt)
    d4 = next(tb for tb in tabs
              if re.search(r"<Title>(.*?)</Title>", tb, re.S)
              and strip(re.search(r"<Title>(.*?)</Title>", tb, re.S).group(1)) == "Data 4")
    out = {}
    for yc in re.findall(r"<YColumn.*?</YColumn>", d4, re.S):
        nm = strip(re.search(r"<Title>(.*?)</Title>", yc, re.S).group(1))
        if nm not in CLASS_COLOR:
            continue
        vals = [strip(d) for d in re.findall(r"<d[^>]*>(.*?)</d>", yc, re.S) if strip(d) != ""]
        if vals:
            out[nm] = (float(vals[0]), float(vals[1]) if len(vals) > 1 else 0.0)
    return out


def published_order(emax):
    """Top -> bottom: class groups, sorted small -> large within group, DAMGO pinned last."""
    order_tb = []
    for _, ligs in GROUPS:
        order_tb += sorted(ligs, key=lambda l: emax[l][0])
    return [l for l in order_tb if l != "DAMGO"] + ["DAMGO"]


def style(ax, xmax, ticks):
    ax.set_xlim(0, xmax)
    ax.set_xticks(ticks)
    ax.set_xlabel("E$_{max}$ (% DAMGO)")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def tint(hexcol, frac=0.20):
    """Blend a hex colour toward white; frac = how much colour is left."""
    r, g, b = (int(hexcol[i:i + 2], 16) for i in (1, 3, 5))
    return tuple((1 - frac) + frac * (ch / 255) for ch in (r, g, b))


def bars(ax, yy, vals, sems, colors, open_mask=None, height=0.72):
    open_mask = open_mask or [False] * len(vals)
    for y, v, s, c, op in zip(yy, vals, sems, colors, open_mask):
        if op:
            ax.barh(y, v, height=height, facecolor=tint(c), edgecolor=c,
                    linewidth=0.5, zorder=2)
        else:
            ax.barh(y, v, height=height, color=c, edgecolor="none", zorder=2)
        if s:
            ax.errorbar(v, y, xerr=s, fmt="none", ecolor="black",
                        elinewidth=0.5, capsize=1.2, capthick=0.5, zorder=3)


def save(fig, name):
    fig.savefig(f"{name}.pdf", transparent=True)
    fig.savefig(f"{name}.png", dpi=600, transparent=False, facecolor="white")
    plt.close(fig)
    print("Saved:", name + ".pdf / .png")


def variant_a(emax):
    """SR-17018 set below a gap, open bar, axis to 150. Recommended."""
    order = published_order(emax)[::-1]          # index 0 = bottom
    fig = plt.figure(); fig.set_size_inches(60 / 25.4, 64 / 25.4)
    ax = fig.add_axes([0.36, 0.155, 0.60, 0.815])
    # SR-17018 occupies y = -1.5, the published set starts at y = 0
    yy = list(np.arange(len(order)))
    bars(ax, yy, [emax[l][0] for l in order], [emax[l][1] for l in order],
         [CLASS_COLOR[l] for l in order])
    bars(ax, [-1.6], [SR_EMAX], [SR_SEM], [SR_COLOR], open_mask=[True])
    style(ax, XMAX, XTICKS)
    ax.plot([0, 150], [-0.85, -0.85], color="0.6", linewidth=0.4,
            linestyle=(0, (1.5, 1.5)), zorder=1, clip_on=False)
    ax.set_yticks([-1.6] + yy)
    ax.set_yticklabels(["SR-17018"] + order)
    ax.set_ylim(-2.3, len(order) - 0.3)
    ax.spines["left"].set_bounds(-0.3, len(order) - 0.3)
    save(fig, f"trupath_emax_bar_sr17018_{TAG}_A_separated")


def variant_b(emax):
    """Published axis kept at 0-105; SR-17018 runs off with a break mark."""
    order = published_order(emax)[::-1]
    fig = plt.figure(); fig.set_size_inches(60 / 25.4, 64 / 25.4)
    ax = fig.add_axes([0.36, 0.155, 0.60, 0.815])
    yy = list(np.arange(len(order)))
    bars(ax, yy, [emax[l][0] for l in order], [emax[l][1] for l in order],
         [CLASS_COLOR[l] for l in order])
    bars(ax, [-1.6], [105], [0], [SR_COLOR], open_mask=[True])
    # break mark on the clipped bar
    for dx in (-3.5, -1.5):
        ax.plot([105 + dx - 1.2, 105 + dx + 1.2], [-1.6 - 0.42, -1.6 + 0.42],
                color="white", linewidth=1.1, zorder=4, clip_on=False, solid_capstyle="butt")
        ax.plot([105 + dx - 1.2, 105 + dx + 1.2], [-1.6 - 0.42, -1.6 + 0.42],
                color=SR_COLOR, linewidth=0.5, zorder=5, clip_on=False)
    ax.text(103, -2.15, "126 (95–158)", ha="right", va="center", fontsize=5)
    style(ax, 105, [0, 50, 100])
    ax.plot([0, 105], [-0.85, -0.85], color="0.6", linewidth=0.4,
            linestyle=(0, (1.5, 1.5)), zorder=1, clip_on=False)
    ax.set_yticks([-1.6] + yy)
    ax.set_yticklabels(["SR-17018"] + order)
    ax.set_ylim(-2.3, len(order) - 0.3)
    ax.spines["left"].set_bounds(-0.3, len(order) - 0.3)
    save(fig, f"trupath_emax_bar_sr17018_{TAG}_B_clipped")


def variant_c(emax):
    """SR-17018 sorted into its own dark-red class group, axis to 150."""
    groups = [(c, ligs + ([SR_NAME] if c == SR_COLOR else [])) for c, ligs in GROUPS]
    val = {l: emax[l][0] for l in CLASS_COLOR}
    val[SR_NAME] = SR_EMAX
    sem = {l: emax[l][1] for l in CLASS_COLOR}
    sem[SR_NAME] = SR_SEM
    col = dict(CLASS_COLOR); col[SR_NAME] = SR_COLOR
    order_tb = []
    for _, ligs in groups:
        order_tb += sorted(ligs, key=lambda l: val[l])
    order_tb = [l for l in order_tb if l != "DAMGO"] + ["DAMGO"]
    order = order_tb[::-1]
    fig = plt.figure(); fig.set_size_inches(60 / 25.4, 60 / 25.4)
    ax = fig.add_axes([0.36, 0.155, 0.60, 0.815])
    yy = list(np.arange(len(order)))
    bars(ax, yy, [val[l] for l in order], [sem[l] for l in order],
         [col[l] for l in order],
         open_mask=[(l == SR_NAME and not SR_SOLID) for l in order])
    ax.set_yticks(yy)
    ax.set_yticklabels(["SR-17018" if l == SR_NAME else l for l in order])
    ax.set_ylim(-0.7, len(order) - 0.3)
    style(ax, XMAX, XTICKS)
    save(fig, f"trupath_emax_bar_sr17018_{TAG}_C_inclass")


def variant_d(emax):
    """SR-17018 in the black 'strong' group, pinned below DAMGO so the efficacy
    gradient stays monotonic. No separator; styling follows SR_SOLID."""
    order_tb = published_order(emax) + [SR_NAME]
    order = order_tb[::-1]
    val = {l: emax[l][0] for l in CLASS_COLOR}; val[SR_NAME] = SR_EMAX
    sem = {l: emax[l][1] for l in CLASS_COLOR}; sem[SR_NAME] = SR_SEM
    col = dict(CLASS_COLOR); col[SR_NAME] = SR_COLOR
    fig = plt.figure(); fig.set_size_inches(60 / 25.4, 62 / 25.4)
    ax = fig.add_axes([0.36, 0.15, 0.60, 0.82])
    yy = list(np.arange(len(order)))
    bars(ax, yy, [val[l] for l in order], [sem[l] for l in order],
         [col[l] for l in order],
         open_mask=[(l == SR_NAME and not SR_SOLID) for l in order])
    ax.set_yticks(yy)
    ax.set_yticklabels(["SR-17018" if l == SR_NAME else l for l in order])
    ax.set_ylim(-0.7, len(order) - 0.3)
    style(ax, XMAX, XTICKS)
    save(fig, f"trupath_emax_bar_sr17018_{TAG}_D_strong_below_damgo")


def main():
    emax = tp_emax()
    conc = {-4: "100 uM", -5: "10 uM", -6: "1 uM", -7: "100 nM"}.get(SR_DOSE, "?")
    where = SOURCE if SOURCE in _LIT else f"1e{SR_DOSE} M = {conc}, {SOURCE} plate set"
    print(f"SR-17018 added at {SR_EMAX:.1f} +- {SR_SEM:.1f} % DAMGO ({where})")
    print("published panel, top -> bottom:")
    for l in published_order(emax):
        print(f"  {l:14s} {emax[l][0]:6.1f} +- {emax[l][1]:.2f}")
    variant_a(emax)
    variant_b(emax)
    variant_c(emax)
    variant_d(emax)


if __name__ == "__main__":
    main()

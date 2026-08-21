#!/usr/bin/env python3
"""
WT TruPath Gi1 BRET dose-response curves, normalized to MATCH the published
peak-based Emax bar (trupath_emax_bar.py / Prism 'Data 4' tab).

Metric = activity as % of the DAMGO peak — the SAME recipe as Data 4, which it
reproduces EXACTLY for all 14 ligands:
  * per replicate, divide by its own -13 (no-drug) baseline -> norm
  * activity = 1 - norm
  * point(dose) = per-dose mean activity / DAMGO_peak_activity * 100
  * Emax(lig)  = max per-dose mean activity / DAMGO_peak * 100   (= the bar value)

Each ligand's curve PLATEAUS EXACTLY at its Emax (= bar):
  y_curve(x) = Emax[lig] / (1 + 10^(le - x)) ,
with le (LogEC50) from a per-ligand Hill=1 fit of the raw BRET (reproduces the
published EC50s). Because the plateau is anchored to the peak Emax rather than a
sigmoid span, the biased partials whose signal declines at supra-max doses
(PZM21 78.7, TRV130 77.2) now top out at their reported efficacy; those
high-dose points simply fall below the plateau. Highest-conc MP point (-4)
dropped as an outlier; '*'-flagged points excluded (as in Data 4).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import re
import numpy as np
from scipy.optimize import curve_fit

PZFX = "20260117_mor_wt_trupath_compilation.pzfx"
OUT  = "trupath_drc_curves.pdf"

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 6, "axes.labelsize": 6, "axes.titlesize": 6,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "text.color": "black", "axes.labelcolor": "black",
    "xtick.color": "black", "ytick.color": "black", "axes.edgecolor": "black",
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "xtick.direction": "out", "ytick.direction": "out",
})

ligand_palette = {
    "DAMGO":"#1F77B4","Fentanyl":"#FF7F0E","Morphine":"#2CA02C","Naloxone":"#D62728",
    "Buprenorphine":"#9467BD","Methadone":"#8C564B","PZM21":"#E377C2","Carfentanil":"#17BECF",
    "Nalbuphine":"#BCBD22","Butorphanol":"#7F7F7F","TRV130":"#AEC7E8","MP":"#FFBB78",
    "Naltrexone":"#843C39","C6Guano":"#5254A3",
}
def color_for(l): return "black" if l == "DAMGO" else ligand_palette[l]
NAMES = ["DAMGO","PZM21","Methadone","Naloxone","Naltrexone","Carfentanil","Morphine",
         "Nalbuphine","Butorphanol","Buprenorphine","C6Guano","TRV130","MP","Fentanyl"]
EXCLUDE = {("MP", -4.0)}       # highest-conc MP point is an outlier

def model(x, top, bot, le): return bot + (top - bot) / (1 + 10.0**(le - x))
def up(x, emax, le):        return emax / (1 + 10.0**(le - x))   # anchored plateau=emax, starts at 0

def parse_raw():
    """Return {ligand: [replicate, ...]} where replicate = [(dose, value, flagged), ...]."""
    txt = open(PZFX, encoding="utf-8", errors="replace").read()
    strip = lambda s: re.sub(r"<[^>]+>", "", s).strip()
    tabs = re.split(r"(?=<Table )", txt)
    rt = next(tb for tb in tabs if re.search(r"<Title>(.*?)</Title>", tb, re.S)
              and strip(re.search(r"<Title>(.*?)</Title>", tb, re.S).group(1)) == "20260117 µOR WT Gi1 TRUPATH")
    def subs(cb):
        return [[(strip(d) if d != "" else "") for d in re.findall(r"<d[^>]*/>|<d[^>]*>(.*?)</d>", s, re.S)]
                for s in re.findall(r"<Subcolumn>(.*?)</Subcolumn>", cb, re.S)]
    Xr = subs(re.search(r"<XColumn.*?</XColumn>", rt, re.S).group(0))[0]
    data = {}
    for yc in re.findall(r"<YColumn.*?</YColumn>", rt, re.S):
        nm = strip(re.search(r"<Title>(.*?)</Title>", yc, re.S).group(1))
        if nm not in ligand_palette:
            continue
        reps = []
        for sub in subs(yc):
            rep = []
            for k, v in enumerate(sub):
                if k < len(Xr) and v != "" and Xr[k] != "":
                    fl = "*" in v
                    try: fv = float(v.replace("*", "").strip())
                    except ValueError: continue
                    x = float(Xr[k])
                    if (nm, x) in EXCLUDE: continue
                    rep.append((x, fv, fl))
            reps.append(rep)
        data[nm] = reps
    return data

def fit_le(reps, lig):
    pts = [(x, y) for rep in reps for (x, y, fl) in rep if not fl]
    x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts])
    p, _ = curve_fit(model, x, y, p0=[np.min(y), np.max(y), -8.0], maxfev=40000)
    return p[2]

def activity_by_dose(reps):
    """Per-dose mean & SEM of activity (1 - value/own_-13_baseline), replicate-wise."""
    per = {}
    for rep in reps:
        base = next((y for (x, y, fl) in rep if x == -13.0 and not fl), None)
        if base in (None, 0):
            continue
        for (x, y, fl) in rep:
            if fl:
                continue
            per.setdefault(x, []).append(1.0 - y / base)
    doses = sorted(per)
    mean = {d: float(np.mean(per[d])) for d in doses}
    sem  = {d: (float(np.std(per[d], ddof=1) / np.sqrt(len(per[d]))) if len(per[d]) > 1 else 0.0) for d in doses}
    return mean, sem

def main():
    data = parse_raw()
    le  = {l: fit_le(data[l], l) for l in NAMES}
    act = {l: activity_by_dose(data[l]) for l in NAMES}
    dpk = max(act["DAMGO"][0].values())                     # DAMGO peak mean activity
    emax = {l: max(act[l][0].values()) / dpk * 100.0 for l in NAMES}   # = Data 4 / bar values
    order = ["DAMGO"] + [l for l in sorted(NAMES, key=lambda l: -emax[l]) if l != "DAMGO"]
    print("TruPath Emax (% DAMGO), peak-based (matches bar / Data 4):")
    for l in order:
        print(f"  {l:14s} {emax[l]:6.1f}   (LogEC50 {le[l]:7.2f})")

    fig = plt.figure(); fig.set_size_inches(78 / 25.4, 50 / 25.4)
    ax = fig.add_axes([0.115, 0.19, 0.55, 0.68])
    bar_hi, bar_lo = [], []
    for lig in order:
        col = color_for(lig); z = 10 if lig == "DAMGO" else 2
        mean, sem = act[lig]
        doses = [d for d in sorted(mean) if d > -13]
        mx = doses
        my = [mean[d] / dpk * 100.0 for d in doses]
        se = [sem[d] / dpk * 100.0 for d in doses]
        xlo = max(-12.5, min(doses) - 0.5); xhi = min(-3.5, max(doses) + 0.5)
        xx = np.linspace(xlo, xhi, 300)
        ax.plot(xx, up(xx, emax[lig], le[lig]), color=col, lw=1.0, zorder=z)
        ax.errorbar(mx, my, yerr=se, fmt="o", ms=2.2, mfc=col, mec="none",
                    ecolor=col, elinewidth=0.5, capsize=1, capthick=0.5, zorder=z + 1)
        bar_hi += [a + b for a, b in zip(my, se)]; bar_lo += [a - b for a, b in zip(my, se)]

    ax.axhline(0, color="grey", lw=0.4, ls=(0, (3, 3)), zorder=0)
    ax.set_xlim(-12.5, -3.5); ax.set_ylim(min(-10, min(bar_lo) - 3), max(bar_hi) + 4)
    ax.set_xticks(range(-12, -3, 1))
    ax.set_xticklabels([str(v) if v % 2 == 0 else "" for v in range(-12, -3, 1)])
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Log[agonist] (M)"); ax.set_ylabel("Activity (% DAMGO)")
    ax.set_title("TRUPATH GI1 BRET assay")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    handles = [Line2D([0], [0], marker="o", linestyle="none", markersize=2.6,
                      markerfacecolor=color_for(l), markeredgecolor="none", label=l) for l in order]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
              handlelength=1.0, handletextpad=0.4, labelspacing=0.35, borderpad=0.0, frameon=False)
    fig.savefig(OUT, transparent=True); print("Saved:", OUT)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Panel c at exactly 40 x 40 mm, 6 pt Helvetica, margins trimmed to the text.

Same content as trupath_emax_bar_sr17018_20260502_EJO_own_vehicle_C_inclass: the 14
published Data 4 values plus SR-17018 at 98.5 +- 19.3 (20260502 rows E/J/O, own
vehicle), sorted into the black strong-agonist group.

Fitting 15 rows into 40 mm leaves ~2.2 mm per bar, so the y labels are set solid with
no leading to spare. Everything that can be reclaimed is: tick marks shortened to
1 pt, label pads to 1, no title, no top/right spines, and the axes rectangle is
computed from the measured width of the longest y label rather than guessed.

Per the exact-mm note: force Agg, set_size_inches, and save WITHOUT bbox_inches so
the MediaBox lands on 40 mm exactly (verified at the end).
"""

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PZFX = HERE / "20260117_mor_wt_trupath_compilation.pzfx"
OUT = HERE / "trupath_emax_bar_40mm"
MM = 40.0
SR_NAME, SR_EMAX, SR_SEM = "SR17018", 98.5, 19.3

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 6, "axes.labelsize": 6, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "text.color": "black", "axes.labelcolor": "black",
    "xtick.color": "black", "ytick.color": "black", "axes.edgecolor": "black",
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 1.0, "ytick.major.size": 1.0,
    "xtick.major.pad": 1.0, "ytick.major.pad": 1.0,
    "axes.labelpad": 1.0,
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


def label_width_mm(labels, fig):
    """Widest rendered y label, in mm, at the current rcParams."""
    r = fig.canvas.get_renderer()
    w = 0.0
    for s in labels:
        t = fig.text(0, 0, s, fontsize=mpl.rcParams["ytick.labelsize"])
        w = max(w, t.get_window_extent(renderer=r).width)
        t.remove()
    return w / fig.dpi * 25.4


def main():
    emax = tp_emax()
    val = {l: v[0] for l, v in emax.items()}
    sem = {l: v[1] for l, v in emax.items()}
    val[SR_NAME], sem[SR_NAME] = SR_EMAX, SR_SEM
    col = dict(CLASS_COLOR)
    col[SR_NAME] = "#000000"

    groups = [(c, ligs + ([SR_NAME] if c == "#000000" else [])) for c, ligs in GROUPS]
    order_tb = []
    for _, ligs in groups:
        order_tb += sorted(ligs, key=lambda l: val[l])
    order_tb = [l for l in order_tb if l != "DAMGO"] + ["DAMGO"]
    order = order_tb[::-1]
    labels = ["SR-17018" if l == SR_NAME else l for l in order]

    fig = plt.figure()
    fig.set_size_inches(MM / 25.4, MM / 25.4)
    fig.canvas.draw()

    lab_mm = label_width_mm(labels, fig)
    left = (lab_mm + 1.0) / MM            # labels + tick + pad
    bottom = 5.3 / MM                     # x ticks + axis label
    right_pad = 0.4 / MM                  # last tick is 100, well inside xlim 120
    top_pad = 0.3 / MM
    ax = fig.add_axes([left, bottom, 1 - left - right_pad, 1 - bottom - top_pad])

    yy = np.arange(len(order))
    for y, l in zip(yy, order):
        ax.barh(y, val[l], height=0.75, color=col[l], edgecolor="none", zorder=2)
        ax.errorbar(val[l], y, xerr=sem[l], fmt="none", ecolor="black",
                    elinewidth=0.4, capsize=0.8, capthick=0.4, zorder=3)
    ax.set_yticks(yy)
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.65, len(order) - 0.35)
    ax.set_xlim(0, 120)   # SR-17018's error bar tops out at 117.8
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel("E$_{max}$ (% DAMGO)")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    fig.savefig(f"{OUT}.pdf", transparent=True)
    fig.savefig(f"{OUT}.png", dpi=1200, facecolor="white")
    plt.close(fig)

    data = open(f"{OUT}.pdf", "rb").read()
    box = re.search(rb"/MediaBox\s*\[([^\]]*)\]", data).group(1).decode().split()
    w, h = float(box[2]) - float(box[0]), float(box[3]) - float(box[1])
    print(f"widest y label {lab_mm:.2f} mm -> left margin {left * MM:.2f} mm")
    print(f"row pitch {(1 - bottom - top_pad) * MM / len(order):.2f} mm per bar")
    print(f"MediaBox {w:.3f} x {h:.3f} pt = {w / 72 * 25.4:.3f} x {h / 72 * 25.4:.3f} mm")
    fonts = sorted(set(f.decode() for f in re.findall(rb"/BaseFont\s*/([A-Za-z0-9+-]+)", data)))
    print("fonts:", fonts)


if __name__ == "__main__":
    main()

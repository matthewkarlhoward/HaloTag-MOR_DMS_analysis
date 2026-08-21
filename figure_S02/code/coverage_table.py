#!/usr/bin/env python3
"""
Library coverage summary table for the baseline sample (b1 + b2).

Columns: Total | Synonymous | Missense | INDEL (= D + I)
Rows:    Designed Variants / Observed Variants
Each cell: raw count and (% of designed in that category)

Output: lib_generation/coverage_table.pdf (60 mm wide x 15 mm tall)
"""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

MM = 1.0 / 25.4

mpl.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 6,
    "text.color": "black",
    "axes.edgecolor": "black",
})

HERE = Path(__file__).resolve().parent
B1 = pd.read_csv(HERE / "baseline_1.csv")
B2 = pd.read_csv(HERE / "baseline_2.csv")
MV = pd.read_csv(HERE / "mor_variants.csv")

b = B1[["mutation_type", "mutation"]].copy()
b["count"] = B1["count"].values + B2["count"].values
obs = b["count"] > 0

# Exclude impossible synonymous variants from the designed set:
# Met (ATG) and Trp (TGG) have only one codon, so a "synonymous" change cannot
# exist and these 20 entries were not designed as oligos.
impossible_syn = (MV["mutation_type"] == "S") & MV["mutation"].isin(["M", "W"])
MV_designed = MV[~impossible_syn]

def category(mt_series, mt):
    if mt == "INDEL":
        return mt_series.isin(["D", "I"])
    return mt_series == mt

cats = [("Total", None), ("Synonymous", "S"), ("Missense", "M"), ("INDEL", "INDEL")]
designed_n = []
observed_n = []
for label, key in cats:
    if key is None:
        d = len(MV_designed)
        o = int(obs.sum())
    else:
        d = int(category(MV_designed["mutation_type"], key).sum())
        o = int((category(b["mutation_type"], key) & obs).sum())
    designed_n.append(d)
    observed_n.append(o)

print("Designed:", dict(zip([c[0] for c in cats], designed_n)))
print("Observed:", dict(zip([c[0] for c in cats], observed_n)))

# Layout (mm)
W, H = 60, 15
label_w = 19
data_w = (W - label_w) / 4
header_h = 4
row_h = (H - header_h) / 2

fig = plt.figure(figsize=(W * MM, H * MM))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis("off")

LW = 0.4
# vertical lines (5 of them, framing 4 data columns)
for i in range(5):
    x = label_w + i * data_w
    ax.plot([x, x], [0, H], color="black", lw=LW, solid_capstyle="butt")
# horizontal lines (top, under header, bottom)
for y in [0, H - header_h, H]:
    ax.plot([label_w, W], [y, y], color="black", lw=LW, solid_capstyle="butt")

def col_xc(i):
    return label_w + (i + 0.5) * data_w

y_header = H - header_h / 2
for i, (label, _) in enumerate(cats):
    ax.text(col_xc(i), y_header, label, ha="center", va="center", fontsize=6)

row_centers_y = [H - header_h - row_h / 2, H - header_h - 1.5 * row_h]
row_labels = ["Designed Variants", "Observed Variants"]
for i, lbl in enumerate(row_labels):
    ax.text(label_w - 1, row_centers_y[i], lbl, ha="right", va="center", fontsize=6)

def fmt(n, pct):
    if pct >= 99.99:
        return f"{n}", "(100%)"
    return f"{n}", f"({pct:.1f}%)"

for r, vals in enumerate([designed_n, observed_n]):
    yc = row_centers_y[r]
    for c, v in enumerate(vals):
        pct = 100.0 if r == 0 else (v / designed_n[c] * 100)
        n_s, p_s = fmt(v, pct)
        xc = col_xc(c)
        ax.text(xc, yc + 1.0, n_s, ha="center", va="center", fontsize=6)
        ax.text(xc, yc - 1.0, p_s, ha="center", va="center", fontsize=6)

out = HERE / "coverage_table.pdf"
fig.savefig(out, dpi=600)
print(f"Wrote {out}")

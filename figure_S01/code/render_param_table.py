#!/usr/bin/env python3
"""
Render param_table.csv as a house-style (Helvetica 6pt, black, booktabs rules)
vector PDF table for placing into Illustrator. Run build_param_table.py first.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42, "mathtext.default": "regular",
    "font.size": 6, "text.color": "black",
})
MINUS = "−"; PM = "±"
def num(v, d):                      # signed number with typographic minus
    if v == "" or v is None: return None
    s = ("%%.%df" % d) % float(v)
    return s.replace("-", MINUS)
def cell(v, s, d):
    if v == "" or v is None: return MINUS
    return f"{num(v,d)} {PM} {num(s,d)}"

rows = list(csv.DictReader(open("param_table.csv")))

# column centers (axis fraction) and alignment
COLS = [
    ("lig",  0.005, "left"),
    ("cE",   0.235, "center"),
    ("cM",   0.400, "center"),
    ("cN",   0.505, "center"),
    ("tE",   0.660, "center"),
    ("tM",   0.825, "center"),
    ("tN",   0.930, "center"),
]
cx = {k: (x, a) for k, x, a in COLS}

fig = plt.figure(); fig.set_size_inches(150 / 25.4, 62 / 25.4)
ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
def T(k, y, s, weight="normal", size=6):
    x, a = cx[k]; ax.text(x, y, s, ha=a, va="center", fontsize=size, fontweight=weight)
def rule(x0, x1, y, lw=0.5):
    ax.plot([x0, x1], [y, y], color="black", lw=lw, solid_capstyle="butt")

y = 0.955; step = 0.0595
# group headers
ax.text((cx["cE"][0] + cx["cN"][0]) / 2, y, "Pooled cAMP", ha="center", va="center", fontsize=6)
ax.text((cx["tE"][0] + cx["tN"][0]) / 2, y, "TRUPATH Gi1", ha="center", va="center", fontsize=6)
rule(0.155, 0.545, y - 0.028); rule(0.575, 0.965, y - 0.028)      # group underlines
# column headers
y -= step
for k, s in [("lig","Ligand"), ("cE","LogEC$_{50}$"), ("cM","E$_{max}$, % DAMGO"), ("cN","n"),
             ("tE","LogEC$_{50}$"), ("tM","E$_{max}$, % DAMGO"), ("tN","n")]:
    T(k, y, s)
rule(0.0, 0.985, 0.985, lw=0.7)                                   # top rule
rule(0.0, 0.985, y - step/2)                                      # header rule
# data rows
for r in rows:
    y -= step
    T("lig", y, r["ligand"])
    T("cE", y, cell(r["camp_logec50"], r["camp_logec50_sem"], 2))
    T("cM", y, cell(r["camp_emax"], r["camp_emax_sem"], 1))
    T("cN", y, r["camp_n"] if r["camp_n"] else MINUS)
    T("tE", y, cell(r["tp_logec50"], r["tp_logec50_sem"], 2))
    T("tM", y, cell(r["tp_emax"], r["tp_emax_sem"], 1))
    T("tN", y, r["tp_n"])
rule(0.0, 0.985, y - step/2, lw=0.7)                              # bottom rule

fig.savefig("param_table.pdf", transparent=True)
print("Saved: param_table.pdf")

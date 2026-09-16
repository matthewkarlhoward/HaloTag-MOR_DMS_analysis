#!/usr/bin/env python3
"""Venn of EC50 vs Emax divergent determinants at the 2 SD cutoff.
   EC50-only 29 | Emax-only 49 | both 27  (EC50 total 56, Emax total 76)."""
from pathlib import Path
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
MM = 1/25.4
OUT = Path("/Users/mkh/GitHub/mor_dms_analysis/plots/ec50_emax_overlap")
BLUE, RED, PURPLE = "#1f4fd8", "#d42b2b", "#7b3fa0"
EC_ONLY, EM_ONLY, BOTH = 29, 49, 27

mpl.rcParams.update({"font.family":"Helvetica","font.size":6,"pdf.fonttype":42})
fig, ax = plt.subplots(figsize=(70*MM, 55*MM))

try:
    from matplotlib_venn import venn2, venn2_circles
    v = venn2(subsets=(EC_ONLY, EM_ONLY, BOTH), set_labels=("EC50", "Emax"), ax=ax)
    for pid, col in [("10", BLUE), ("01", RED), ("11", PURPLE)]:
        p = v.get_patch_by_id(pid)
        if p: p.set_color(col); p.set_alpha(0.85); p.set_edgecolor("white"); p.set_linewidth(0.8)
    for lid in ("10", "01", "11"):
        lbl = v.get_label_by_id(lid)
        if lbl: lbl.set_color("white"); lbl.set_fontsize(8); lbl.set_fontweight("bold")
    for sl in v.set_labels or []:
        if sl: sl.set_fontsize(7)
    venn2_circles(subsets=(EC_ONLY, EM_ONLY, BOTH), ax=ax, lw=0.5, color="0.4")
    mode = "matplotlib_venn"
except Exception as e:
    # fallback: two overlapping circles by hand
    from matplotlib.patches import Circle
    ax.add_patch(Circle((-0.42,0), 0.62, fc=BLUE, ec="0.4", lw=0.5, alpha=0.55))
    ax.add_patch(Circle(( 0.42,0), 0.72, fc=RED,  ec="0.4", lw=0.5, alpha=0.55))
    ax.text(-0.62,0, str(EC_ONLY), ha="center", va="center", color="white", fontsize=8, fontweight="bold")
    ax.text( 0.62,0, str(EM_ONLY), ha="center", va="center", color="white", fontsize=8, fontweight="bold")
    ax.text( 0.02,0, str(BOTH),   ha="center", va="center", color="white", fontsize=8, fontweight="bold")
    ax.text(-0.6,0.78,"EC50", ha="center", fontsize=7); ax.text(0.6,0.86,"Emax", ha="center", fontsize=7)
    ax.set_xlim(-1.3,1.4); ax.set_ylim(-0.95,1.05); ax.set_aspect("equal")
    mode = "manual fallback"

ax.set_title("Morphine determinants: EC50 (potency) vs Emax (efficacy), >2 SD", fontsize=6)
ax.axis("off")
fig.tight_layout()
fig.savefig(OUT/"ec50_emax_venn_2SD.pdf", transparent=True)
fig.savefig(OUT/"ec50_emax_venn_2SD.png", dpi=300, transparent=False)
print("rendered via:", mode)
print(f"EC50 total {EC_ONLY+BOTH} | Emax total {EM_ONLY+BOTH} | shared {BOTH}")

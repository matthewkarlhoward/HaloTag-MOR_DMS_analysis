#!/usr/bin/env python3
"""
Per-variant beeswarm of EC50 and Emax at a single position (default V83),
Morphine only. Two panels side by side.

One point = one high-R^2 sigmoid missense variant at the position (same set as
the V83 dose-response showcase: fit_ok, missense, sigmoid, r2 >= R2_MIN).
  - left  : EC50 (log[M]), absolute fitted value  — more negative = more potent
  - right : Emax on the Activity scale (syn Emin -> 0, syn Emax -> 1), so 1 is
            WT-like efficacy and 0 is baseline — higher = more efficacious
A thin grey dashed line marks the synonymous mean in each panel (EC50 syn mean;
Emax = 1).

Points are colored by EC50 (magma), with the SAME color in both panels, so a
variant's potency can be tracked against its efficacy.

Reads:  curve_refitting/refit_3param_robust_morphine.csv
Outputs: curve_refitting/position_beeswarm_morphine_<WTpos>.{pdf,png}
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT  = ROOT / "curve_refitting"

POSITION = 83
R2_MIN   = 0.85

plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 6,
    "axes.labelsize": 6, "axes.titlesize": 6,
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 5,
    "text.color": "black", "axes.edgecolor": "black",
    "axes.labelcolor": "black", "axes.titlecolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "lines.linewidth": 0.5, "patch.linewidth": 0.5,
})

def beeswarm_positions(values, width=0.6, col_gap=0.09):
    """Greedy beeswarm with FIXED column spacing. Place y-sorted points at the
    smallest |x| that keeps every marker >= one diameter d (in y-units) from its
    placed neighbours (true 2-D distance), then map each unit of d to a fixed
    axis gap (col_gap) so dense clusters stay tight and contiguous. Two failure
    modes this avoids: an x-only exclusion test (step < exclusion radius) leaves
    a gap after the center column -> looks like two columns; normalizing offsets
    to fill width/2 over-spreads a few rings -> also looks like separate columns.
    """
    values = np.asarray(values, dtype=float)
    n = len(values)
    if n == 0:
        return np.array([])
    vrange = values.max() - values.min()
    d = vrange / 15 if vrange > 0 else 1.0
    step = d * 0.5
    order = np.argsort(values, kind="mergesort")
    xo = np.zeros(n)
    placed = []
    for idx in order:
        y = values[idx]
        near = [(px, py) for px, py in placed if abs(py - y) < d]
        k = 0
        while True:
            cands = (0.0,) if k == 0 else (k * step, -k * step)
            hit = next((cx for cx in cands
                        if all((cx - px) ** 2 + (y - py) ** 2 >= d * d
                               for px, py in near)), None)
            if hit is not None:
                xo[idx] = hit
                placed.append((hit, y))
                break
            k += 1
    half = width / 2
    return np.clip(xo / d * col_gap, -half, half)

# ── Load + synonymous reference ──────────────────────────────────────────────
fit = pd.read_csv(OUT / "refit_3param_robust_morphine.csv")
syn = fit[(fit.fit_ok) & (fit.type == "synonymous")
          & (fit.new_curve_type == "sigmoid")]
syn_emin = float(syn.fitted_emin.mean())
syn_emax = float(syn.fitted_emax.mean())
syn_ec50 = float(syn.fitted_ec50_logM.mean())

sel = fit[(fit.fit_ok) & (fit.position == POSITION)
          & (fit.type == "missense") & (fit.new_curve_type == "sigmoid")
          & (fit.r2 >= R2_MIN)].copy()
wt = sel.wildtype.iloc[0]
n = len(sel)

ec50 = sel.fitted_ec50_logM.to_numpy()
# Emax on the Activity scale: syn Emin -> 0, syn Emax -> 1.
emax_act = (sel.fitted_emax.to_numpy() - syn_emin) / (syn_emax - syn_emin)
muts = sel.mutation.to_numpy()

# Color every variant by EC50 (magma), shared across both panels.
cnorm = (ec50 - ec50.min()) / (ec50.max() - ec50.min()) if n > 1 else np.zeros(n)
colors = plt.get_cmap("magma")(cnorm * 0.88)

PANELS = [
    (ec50,     syn_ec50, "EC50 (log[M])"),
    (emax_act, 1.0,      "Emax (Activity)"),
]

MM = 1 / 25.4
fig, axes = plt.subplots(1, 2, figsize=(50 * MM, 40 * MM),
                         constrained_layout=True)

for ax, (yv, synref, ylab) in zip(axes, PANELS):
    ax.axhline(synref, color="grey", lw=0.5, ls="--", zorder=1)
    xs = beeswarm_positions(yv, width=0.6)
    for x, y, c, m in zip(xs, yv, colors, muts):
        ax.text(x, y, m, color=c, ha="center", va="center",
                fontsize=5, fontweight="bold", zorder=3)
    # Text artists are ignored by autoscale, so set y-limits explicitly.
    lo, hi = float(np.min(yv)), float(np.max(yv))
    pad = 0.10 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlim(-0.5, 0.5)
    ax.set_xticks([])
    ax.set_ylabel(ylab, labelpad=1)
    ax.tick_params(axis="y", pad=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)

fig.suptitle(f"{wt}{POSITION}  (n={n})", fontsize=5.5)

tag = f"{wt}{POSITION}"
fig.savefig(OUT / f"position_beeswarm_morphine_{tag}.pdf", dpi=600)
fig.savefig(OUT / f"position_beeswarm_morphine_{tag}.png", dpi=600)
plt.close(fig)
print(f"Saved: position_beeswarm_morphine_{tag}.pdf  (n={n}; "
      f"EC50 {ec50.min():.2f}..{ec50.max():.2f}, "
      f"Emax(Activity) {emax_act.min():.2f}..{emax_act.max():.2f})")

#!/usr/bin/env python3
"""
Combined single-position figure (Morphine, default V83), 80 x 40 mm:
  [ dose-response curves + legend ]  [ EC50 beeswarm ]  [ Emax beeswarm ]

Left: per-dose Activity points + robust 3-parameter Hill fit for every high-R^2
sigmoid missense variant at the position, over the synonymous 95% band + mean.
Variants are colored by EC50 (magma, value-based). The legend lists each variant
by its single-letter mutant code only (color = its EC50 color).

Right: the SAME variants as beeswarms of EC50 (log[M]) and Emax (Activity scale,
syn Emin -> 0, syn Emax -> 1). Each variant is drawn as its single-letter code in
the same color it has in the curves, so potency/efficacy track the curve panel.

Reads:
  dms_scores/composite_dms_scores.csv               (raw per-dose rescaled scores)
  curve_refitting/refit_3param_robust_morphine.csv  (fitted params, class, R^2)
Outputs:
  curve_refitting/position_combo_morphine_<WTpos>.{pdf,png}
"""
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT  = ROOT / "curve_refitting"

POSITION = 83
R2_MIN   = 0.85

plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 6,
    "mathtext.default": "regular",
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

def hill3(x, emin, emax, log_ec50):
    return emin + (emax - emin) / (1 + 10**(log_ec50 - x))

def beeswarm_positions(values, width=0.6, col_gap=0.09):
    """Greedy beeswarm with FIXED column spacing (see position_beeswarm_morphine)."""
    values = np.asarray(values, dtype=float)
    n = len(values)
    if n == 0:
        return np.array([])
    vrange = values.max() - values.min()
    d = vrange / 13 if vrange > 0 else 1.0
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

# ── Load data ────────────────────────────────────────────────────────────────
dms = pd.read_csv(ROOT / "dms_scores/composite_dms_scores.csv")
pat = re.compile(r"MOR_Morphine_(-\d+)M_rescaled_score$")
dose_pairs = sorted(((c, int(pat.search(c).group(1)) / 10.0)
                     for c in dms.columns if pat.search(c)),
                    key=lambda t: t[1])
dose_cols = [c for c, _ in dose_pairs]
doses     = np.array([d for _, d in dose_pairs])
dms_idx   = dms.set_index("hgvs")

fit = pd.read_csv(OUT / "refit_3param_robust_morphine.csv")
syn = fit[(fit.fit_ok) & (fit.type == "synonymous")
          & (fit.new_curve_type == "sigmoid")]
syn_emin = float(syn.fitted_emin.mean())
syn_emax = float(syn.fitted_emax.mean())
syn_ec50 = float(syn.fitted_ec50_logM.mean())

def norm(y):
    return (y - syn_emin) / (syn_emax - syn_emin)

sel = fit[(fit.fit_ok) & (fit.position == POSITION)
          & (fit.type == "missense") & (fit.new_curve_type == "sigmoid")
          & (fit.r2 >= R2_MIN)].sort_values("fitted_ec50_logM").reset_index(drop=True)
wt = sel.wildtype.iloc[0]
n  = len(sel)

ec50     = sel.fitted_ec50_logM.to_numpy()
emax_act = (sel.fitted_emax.to_numpy() - syn_emin) / (syn_emax - syn_emin)
muts     = sel.mutation.to_numpy()

# One magma color per variant, by EC50 value, reused in curves + both beeswarms.
cnorm  = (ec50 - ec50.min()) / (ec50.max() - ec50.min()) if n > 1 else np.zeros(n)
colors = plt.get_cmap("magma")(cnorm * 0.88)

# Dense x-grid for drawing the fitted Hill curves.
xs = np.linspace(-12, -4, 200)

# ── Figure ───────────────────────────────────────────────────────────────────
MM = 1 / 25.4
FIGW, FIGH = 80, 40
fig = plt.figure(figsize=(FIGW * MM, FIGH * MM))

# Shared vertical span: panels fill the height and the legend matches it.
Y0, H = 0.165, 0.795
# Horizontal split of the FILE width: a curve+legend block, then two even swarms.
mmf  = lambda v: v / FIGW            # millimetres -> figure-width fraction
B1R  = 0.56                          # right edge of the curve+legend block
LEGW = mmf(5.5)                      # legend strip at the right of block 1
LAB  = mmf(5.5)                      # y-label + tick allowance left of each swarm
SW   = mmf(10.5)                     # each swarm's plot width; laid out from the right
leg_x0 = B1R - LEGW
ax_c   = fig.add_axes([mmf(8.5), Y0, leg_x0 - mmf(8.5) - mmf(0.5), H])
ax_leg = fig.add_axes([leg_x0, Y0, LEGW, H]); ax_leg.axis("off")
ax_e   = fig.add_axes([1 - mmf(0.3) - 2 * SW - LAB - mmf(2.0), Y0, SW, H])
ax_m   = fig.add_axes([1 - mmf(0.3) - SW, Y0, SW, H])

# Curves (colored by EC50; syn mean for reference, not in legend).
ax_c.plot(xs, norm(hill3(xs, syn_emin, syn_emax, syn_ec50)),
          color="black", lw=0.7, ls="--", zorder=2)
handles = [Line2D([0], [0], color="black", lw=0.7, ls="--", label="syn")]
for color, r in zip(colors, sel.itertuples(index=False)):
    yraw = dms_idx.loc[r.hgvs, dose_cols].to_numpy(dtype=float)
    ax_c.scatter(doses, norm(yraw), s=4, color=color, edgecolors="none",
                 alpha=0.9, zorder=4)
    ax_c.plot(xs, norm(hill3(xs, r.fitted_emin, r.fitted_emax, r.fitted_ec50_logM)),
              color=color, lw=0.8, zorder=3)
    handles.append(Line2D([0], [0], color=color, lw=1.0, label=r.mutation))
ax_c.set_xlabel("log[morphine]", labelpad=1)
ax_c.set_ylabel("Activity", labelpad=1)
ax_c.set_xlim(-11.8, -4.2)
# Tick mark at every molar decade; label alternate ones so they don't collide.
_xt = list(range(-11, -4))
ax_c.set_xticks(_xt)
ax_c.set_xticklabels([str(t) if t % 2 else "" for t in _xt])
ax_c.tick_params(axis="both", pad=1)
ax_c.spines["top"].set_visible(False)
ax_c.spines["right"].set_visible(False)
for s in ("left", "bottom"):
    ax_c.spines[s].set_linewidth(0.5)

# Legend: single-letter mutant code only, in each variant's color.
ax_leg.legend(handles=handles, loc="center left", bbox_to_anchor=(0.0, 0.5),
              frameon=False, handlelength=0.9, handletextpad=0.3,
              labelspacing=0.12, ncol=1, borderaxespad=0)

# Beeswarms: same variants, same colors, drawn as single-letter codes.
for ax, yv, synref, ylab in [(ax_e, ec50, syn_ec50, r"EC$_{50}$"),
                             (ax_m, emax_act, 1.0, r"E$_{max}$")]:
    ax.axhline(synref, color="black", lw=0.5, ls="--", zorder=1)
    bx = beeswarm_positions(yv, width=1.8, col_gap=0.28)
    for x, y, c, m in zip(bx, yv, colors, muts):
        ax.text(x, y, m, color=c, ha="center", va="center",
                fontsize=5.5, fontweight="bold", zorder=3,
                path_effects=[pe.withStroke(linewidth=0.5, foreground="0.2")])
    lo, hi = float(np.min(yv)), float(np.max(yv))
    pad = 0.10 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlim(-0.9, 0.9)
    ax.set_xticks([])
    ax.set_ylabel(ylab, labelpad=1)
    ax.tick_params(axis="y", pad=1)
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)

tag = f"{wt}{POSITION}"
fig.savefig(OUT / f"position_combo_morphine_{tag}.pdf", dpi=600)
fig.savefig(OUT / f"position_combo_morphine_{tag}.png", dpi=600)
plt.close(fig)
print(f"Saved: position_combo_morphine_{tag}.pdf  (n={n}; "
      f"EC50 {ec50.min():.2f}..{ec50.max():.2f}, "
      f"Emax(Activity) {emax_act.min():.2f}..{emax_act.max():.2f})")

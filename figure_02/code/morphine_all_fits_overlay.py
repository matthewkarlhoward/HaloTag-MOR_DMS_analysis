#!/usr/bin/env python3
"""
All robust 3-parameter refit curves for morphine, overlaid (40 x 40 mm).

Reads refit_3param_robust_morphine.csv and draws the fitted Hill curve for
every variant that converged (fit_ok = True). Responses are normalized to the
synonymous (WT-like) reference: mean sigmoid-syn Emin -> 0, mean sigmoid-syn
Emax -> 1. Pathological fits are dropped (params pinned at the +-2 bound, or
inverted curves dropping below -1); super-responders are kept and clipped by
the y-axis limit [-0.4, 1.5]. Thin, low-alpha black traces so the
density/envelope of fits is visible. The synonymous (WT-like) average curve
(mean sigmoid-syn Emin/Emax/EC50) is overlaid on top in red as the reference.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT  = ROOT / "curve_refitting"

plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 6,
    "axes.labelsize": 6, "axes.titlesize": 6,
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "text.color": "black", "axes.edgecolor": "black",
    "axes.labelcolor": "black", "axes.titlecolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "lines.linewidth": 0.5, "patch.linewidth": 0.5,
})

def hill(x, emin, emax, ec50_logM):
    return emin + (emax - emin) / (1.0 + 10.0**(ec50_logM - x))

df = pd.read_csv(OUT / "refit_3param_robust_morphine.csv")
df = df[df.fit_ok == True]

# Synonymous (WT-like) reference for normalization: sigmoid syn only.
syn = df[(df.type == "synonymous") & (df.new_curve_type == "sigmoid")]
syn_emin = float(syn.fitted_emin.mean())
syn_emax = float(syn.fitted_emax.mean())
syn_ec50 = float(syn.fitted_ec50_logM.mean())

xs = np.linspace(-13, -4, 200)

MM = 1 / 25.4
fig, ax = plt.subplots(figsize=(40 * MM, 40 * MM))

n_plotted = 0
for r in df.itertuples(index=False):
    # syn emax is more negative than syn emin, so the negative denominator
    # also flips the curves positive-going (no explicit *-1 needed).
    ys = (hill(xs, r.fitted_emin, r.fitted_emax, r.fitted_ec50_logM)
          - syn_emin) / (syn_emax - syn_emin)
    # Drop pathological fits only: params pinned at the +-2 bound, or inverted
    # curves dropping below -1. Super-responders (>1) are kept; the axis clips them.
    if r.hit_param_bound or ys.min() < -1:
        continue
    ax.plot(xs, ys, color="black", alpha=0.03, linewidth=0.2)
    n_plotted += 1

# Synonymous (WT-like) average curve, drawn last so it sits on top, in red.
ys_syn = (hill(xs, syn_emin, syn_emax, syn_ec50) - syn_emin) / (syn_emax - syn_emin)
ax.plot(xs, ys_syn, color="red", linewidth=0.9, zorder=5)

ax.set_xlabel("log[morphine]", labelpad=1)
ax.set_ylabel("Activity", labelpad=1)
ax.set_ylim(-0.4, 1.5)
ax.tick_params(axis="both", pad=1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.5)

fig.subplots_adjust(left=0.21, bottom=0.16, right=0.97, top=0.98)
out_pdf = OUT / "morphine_all_fits_overlay.pdf"
fig.savefig(out_pdf, dpi=600)
print(f"Saved: {out_pdf}  (n = {n_plotted} of {len(df)} curves; "
      f"syn_emin={syn_emin:.3f}->0, syn_emax={syn_emax:.3f}->1)")
plt.close()

#!/usr/bin/env python3
"""
Variant-level scatter of Morphine EC50 (x, log[M]) vs Emax (y, Activity scale),
from the robust 3-parameter refit (refit_3param_robust_morphine.csv).

One point = one high-quality sigmoid missense variant (fit_ok, missense,
new_curve_type == "sigmoid"). Emax is on the Activity scale: the synonymous
(WT-like) mean Emin -> 0 and mean Emax -> 1, so 1 is WT-like efficacy and 0 is
baseline. Points are colored by local 2-D density (gaussian KDE) and drawn
dense-on-top, so the crowded core is distinguishable from the sparse tails.

Grey dotted reference lines mark the synonymous mean EC50 (vertical) and the
WT-like Emax of 1.0 (horizontal).

QC: missense are filtered to EC50 >= -10.5 (at/above the 2nd-lowest tested dose,
so potency is anchored by data, not extrapolated off the low-dose edge) and
Emax >= 0 (drops sub-baseline noise fits). Same cuts as
morphine_variant_ec50_vs_emax_TM.py and ..._byregion.py.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from pathlib import Path

HERE = Path(__file__).resolve().parent
EC50_FLOOR = -10.5          # 2nd-lowest tested dose; EC50 below this is unanchored
EMAX_FLOOR = 0.0            # Activity scale; fits below the syn baseline are noise

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

# ── Load + synonymous (WT-like) reference ────────────────────────────────────
df = pd.read_csv(HERE / "refit_3param_robust_morphine.csv")
df = df[df.fit_ok & (df.new_curve_type == "sigmoid")]

syn = df[df.type == "synonymous"]
syn_emin = float(syn.fitted_emin.mean())
syn_emax = float(syn.fitted_emax.mean())
syn_ec50 = float(syn.fitted_ec50_logM.mean())

mis = df[df.type == "missense"].copy()
# Emax on the Activity scale: syn Emin -> 0, syn Emax -> 1.
mis["emax_act"] = (mis.fitted_emax - syn_emin) / (syn_emax - syn_emin)
# QC: anchored potency + non-negative efficacy (see module docstring).
n_all = len(mis)
mis = mis[(mis.fitted_ec50_logM >= EC50_FLOOR) & (mis.emax_act >= EMAX_FLOOR)]
ec50 = mis.fitted_ec50_logM.to_numpy()
emax = mis.emax_act.to_numpy()

# ── Local 2-D density (KDE); sort so dense points draw on top ─────────────────
xy = np.vstack([ec50, emax])
z = gaussian_kde(xy)(xy)
order = z.argsort()
ec50, emax, z = ec50[order], emax[order], z[order]

MM = 1 / 25.4
fig, ax = plt.subplots(figsize=(40 * MM, 40 * MM))

ax.axvline(syn_ec50, color="grey", lw=0.5, ls=":", zorder=1)
ax.axhline(1.0, color="grey", lw=0.5, ls=":", zorder=1)

ax.scatter(ec50, emax, c=z, s=1, cmap="magma", edgecolors="none",
           alpha=0.85, zorder=2)

ax.set_xlabel("EC50 (log[M])", labelpad=1)
ax.set_ylabel("Emax (Activity)", labelpad=1)
ax.tick_params(axis="both", pad=1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.5)

fig.subplots_adjust(left=0.205, bottom=0.175, right=0.97, top=0.98)
out_pdf = HERE / "morphine_variant_ec50_vs_emax.pdf"
fig.savefig(out_pdf, dpi=600)
fig.savefig(out_pdf.with_suffix(".png"), dpi=600)
print(f"Saved: {out_pdf}  (n = {len(ec50)} kept, {n_all - len(ec50)} removed; "
      f"syn_ec50 = {syn_ec50:.2f}, syn Emin/Emax = {syn_emin:.3f}/{syn_emax:.3f})")
plt.close()

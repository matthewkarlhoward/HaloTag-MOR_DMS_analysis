#!/usr/bin/env python3
"""
Supplemental operational-model figure (morphine), assembled as a 3-column grid.
Top-left cell is left BLANK for manually-added text / equations.

Layout (3 cols x 4 rows):
  row0: [ BLANK ]                 Potency concordance      Efficacy concordance
  row1: Basal concordance         Surface A (Emax+surface) Surface B (decomposition)
  row2: Efficacy loss top-down    Efficacy loss side       [efficacy colorbar]
  row3: Potency loss top-down     Potency loss side        [potency colorbar]

Panels reproduce 06_hill_concordance.py, 12_surface_incorporation.py, and
13_structure_maps.py. Output: operational_model/plots/fig_supp_operational.*
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, spearmanr, pearsonr

ROOT = Path(__file__).resolve().parents[1]
NET = ROOT / "network_tools"
OUT = ROOT / "operational_model" / "plots"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 6, "mathtext.default": "regular",
    "text.color": "black", "axes.labelcolor": "black", "axes.edgecolor": "black",
    "axes.titlesize": 6, "axes.labelsize": 6, "xtick.labelsize": 5,
    "ytick.labelsize": 5, "xtick.color": "black", "ytick.color": "black",
    "axes.linewidth": 0.5, "lines.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
})
MM = 1 / 25.4

# ══ data ══════════════════════════════════════════════════════════════════════
op = pd.read_csv(ROOT / "operational_model" / "morphine_operational_fits.csv")
op["hill_activation"] = op.hill_emin - op.hill_emax
rob = pd.read_csv(ROOT / "curve_refitting" / "refit_3param_robust_morphine.csv")
con = op.merge(rob[["hgvs", "fitted_ec50_logM", "fitted_emin", "new_curve_type"]],
               on="hgvs", how="inner")
con = con[(con.new_curve_type == "sigmoid") & (con.curve_type == "sigmoid")
          & (con.rmse < 0.08) & ~con.saturating.astype(bool)]
clean = op[(op.curve_type == "sigmoid") & (op.rmse < 0.12)
           & ~op.saturating.astype(bool) & op.type.eq("missense")].copy()

# per-position operational loss (for the structure maps)
syn = op[(op.type == "synonymous") & op.tau.notna() & ~op.saturating.astype(bool)]
EMAX_WT, EC50_WT = float(np.nanmedian(syn.emax_obs_fit)), float(np.nanmedian(syn.logEC50))
ppos = clean.groupby("position").agg(emax=("emax_obs_fit", "mean"),
                                     ec50=("logEC50", "mean")).reset_index()
ppos["eff_loss"] = ((EMAX_WT - ppos.emax) / EMAX_WT).clip(lower=0)   # fraction of WT
ppos["pot_loss"] = (ppos.ec50 - EC50_WT).clip(lower=0)
VAL = {"eff_loss": dict(zip(ppos.position, ppos.eff_loss)),
       "pot_loss": dict(zip(ppos.position, ppos.pot_loss))}
EFF_MAX = float(np.nanquantile(ppos.eff_loss, 0.97))
POT_MAX = float(np.nanquantile(ppos.pot_loss, 0.97))

# ══ structure projection (PCA frame + your side-view camera) ══════════════════
ca = pd.read_csv(NET / "pdb_8efq_chain_R_ca.csv"); ca["position"] = ca.position.astype(int)
X = ca[["x", "y", "z"]].to_numpy(float); X -= X.mean(0)
_, _, Vt = np.linalg.svd(X, full_matrices=False); P = X @ Vt.T
idx = {int(p): i for i, p in enumerate(ca.position)}
sgn = 1.0
if np.mean([P[idx[p], 0] for p in (73, 76, 200, 210) if p in idx]) < \
   np.mean([P[idx[p], 0] for p in (167, 262, 344, 348) if p in idx]):
    P[:, 0] *= -1; sgn = -1.0
coord = {int(p): (P[i, 1], P[i, 2], P[i, 0]) for i, p in enumerate(ca.position)}
ALL = np.array([coord[p] for p in ca.position])
A_L2P = np.array([Vt[1], Vt[2], sgn * Vt[0]])
CAM = np.array([0.17894, -0.015396, 0.98374, 178.08, 0.22556, 0.97389, -0.025787,
                -3.0794, -0.95765, 0.22651, 0.17774, 34.126]).reshape(3, 4)
MOD = np.array([0.96975, -0.18802, 0.15566, -135.48, -0.073502, 0.38316, 0.92075,
                -186.46, -0.23276, -0.90434, 0.35775, 95.405]).reshape(3, 4)
_N = CAM[:, :3].T @ MOD[:, :3]
SIDE = (A_L2P @ _N[0], A_L2P @ _N[1])
TOPD = (np.eye(3)[0], np.eye(3)[1])

# ══ panel helpers ═════════════════════════════════════════════════════════════
def concordance(ax, x, y, xlab, ylab, title, lims):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y); x, y = x[ok], y[ok]
    ax.plot(lims, lims, color="grey", ls=":", lw=0.5, zorder=1)
    z = gaussian_kde(np.vstack([x, y]))(np.vstack([x, y])); o = z.argsort()
    ax.scatter(x[o], y[o], c=z[o], s=1.2, cmap="magma", lw=0, alpha=0.85, zorder=2)
    ax.set_xlim(lims); ax.set_ylim(lims); ax.set_box_aspect(1)
    ax.set_xlabel(xlab, labelpad=1); ax.set_ylabel(ylab, labelpad=1)
    ax.set_title(f"{title}\nρ={spearmanr(x, y).statistic:.2f}, "
                 f"r={pearsonr(x, y).statistic:.2f}, n={len(x)}", pad=2)
    ax.tick_params(pad=1)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def struct(ax, metric, cmap, vmax, view, title=None):
    rx, uy = view; B = np.column_stack([rx, uy]); A2 = ALL @ B
    lo, hi = A2.min(0), A2.max(0); pad = 0.06 * (hi - lo); lo -= pad; hi += pad
    span = max(hi - lo); ctr = (lo + hi) / 2
    ax.set_xlim(ctr[0] - span / 2, ctr[0] + span / 2)
    ax.set_ylim(ctr[1] - span / 2, ctr[1] + span / 2); ax.set_aspect("equal")
    ax.scatter(A2[:, 0], A2[:, 1], s=1.6, color="#dcdcdc", edgecolors="none", zorder=0)
    mp = [p for p in coord if p in VAL[metric]]
    xs = np.array([coord[p] for p in mp]) @ B
    cs = np.array([VAL[metric][p] for p in mp]); o = np.argsort(cs)
    sc = ax.scatter(xs[o, 0], xs[o, 1], c=cs[o], s=7, cmap=cmap, vmin=0, vmax=vmax,
                    edgecolors="#555", linewidths=0.15, zorder=3)
    ax.set_xticks([-20, 0, 20]); ax.set_yticks([-20, 0, 20]); ax.tick_params(pad=1)
    ax.set_xlabel(f"{'x'} (Å)", labelpad=1)
    ax.set_ylabel(("y" if view is TOPD else "z") + " (Å)", labelpad=1)
    if title:
        ax.set_title(title, pad=2)
    return sc


# ══ assemble 3-col grid ═══════════════════════════════════════════════════════
fig = plt.figure(figsize=(180 * MM, 215 * MM))
gs = fig.add_gridspec(4, 3, hspace=0.55, wspace=0.5,
                      left=0.08, right=0.95, top=0.96, bottom=0.05)

# row0: basal | potency | efficacy(observed)
concordance(fig.add_subplot(gs[0, 0]), con.fitted_emin, con.B,
            "Hill E$_{min}$", "operational baseline", "Basal", (-0.5, 0.5))
concordance(fig.add_subplot(gs[0, 1]), con.fitted_ec50_logM, con.logEC50,
            "Hill logEC$_{50}$ (M)", "operational logEC$_{50}$ (M)", "Potency", (-11, -4))
concordance(fig.add_subplot(gs[0, 2]), con.hill_activation, con.emax_obs_fit,
            "Hill E$_{max}$", "operational E$_{max}$", "Efficacy (observed)", (0, 0.75))

# row1: surface A | surface B | reclassification bar
axA = fig.add_subplot(gs[1, 0])
a = clean.dropna(subset=["hill_activation", "emax_obs_fit", "expr_effect"])
sA = axA.scatter(a.hill_activation, a.emax_obs_fit, c=a.expr_effect, cmap="viridis",
                 vmin=np.nanquantile(a.expr_effect, .02), vmax=np.nanquantile(a.expr_effect, .98),
                 s=2.5, lw=0, alpha=0.8)
axA.plot((0, .75), (0, .75), color="grey", ls=":", lw=0.5)
axA.set_xlim(0, .75); axA.set_ylim(0, .75); axA.set_box_aspect(1)
axA.set_xlabel("Hill E$_{max}$", labelpad=1); axA.set_ylabel("operational E$_{max}$", labelpad=1)
axA.set_title("Efficacy + independent\nsurface measurement", pad=2)
for s in ("top", "right"): axA.spines[s].set_visible(False)
cA = fig.colorbar(sA, ax=axA, fraction=0.05, pad=0.03); cA.set_label("surface score", fontsize=5)
cA.ax.tick_params(labelsize=4.5, width=0.4, length=1.5)

axB = fig.add_subplot(gs[1, 1])
b = clean.dropna(subset=["log10_S", "dlog10_tau", "log10_rho"])   # for the bar
dER = clean.dropna(subset=["emax_obs_fit", "log10_rho", "log10_S"]).copy()
dER["emax_wt1"] = dER.emax_obs_fit / EMAX_WT
o = dER.log10_S.abs().argsort()
sB = axB.scatter(dER.emax_wt1.values[o], dER.log10_rho.values[o], c=dER.log10_S.values[o],
                 cmap="PuOr", vmin=-0.4, vmax=0.4, s=2.5, lw=0, alpha=0.85)
axB.axvline(1.0, color="grey", ls=":", lw=0.5); axB.axhline(0.0, color="grey", ls=":", lw=0.5)
axB.text(0.30, 0.18, "low E$_{max}$, ρ≈0\nexpression-limited", fontsize=4.3, color="#8c510a", va="bottom")
axB.text(0.30, -0.78, "low E$_{max}$, ρ<0\ncoupling loss", fontsize=4.3, color="#542788", va="bottom")
axB.set_xlim(0.25, 1.45); axB.set_ylim(-1.5, 1.5); axB.set_box_aspect(1)
axB.set_xlabel("E$_{max}$ (efficacy, WT=1)", labelpad=1)
axB.set_ylabel("log$_{10}$ ρ (coupling)", labelpad=1)
axB.set_title("Efficacy vs coupling,\ncolored by expression", pad=2)
for s in ("top", "right"): axB.spines[s].set_visible(False)
cB = fig.colorbar(sB, ax=axB, fraction=0.05, pad=0.03); cB.set_label("log$_{10}$ S", fontsize=5)
cB.ax.tick_params(labelsize=4.5, width=0.4, length=1.5)

# reclassification of efficacy-loss variants: expression-limited vs coupling loss
axC = fig.add_subplot(gs[1, 2])
lofb = b[b.dlog10_tau < -0.10]
expr_lim = (lofb.log10_rho > -0.05).mean() * 100
coupling = (lofb.log10_rho < -0.15).mean() * 100
inter = 100 - expr_lim - coupling
vals = [expr_lim, inter, coupling]
cols = ["#2c7fb8", "#bdbdbd", "#b2182b"]
bars = axC.bar(range(3), vals, color=cols, edgecolor="black", linewidth=0.4, width=0.72)
for x, v in zip(range(3), vals):
    axC.text(x, v + 1.5, f"{v:.0f}%", ha="center", va="bottom", fontsize=5)
axC.set_xticks(range(3))
axC.set_xticklabels(["expression-\nlimited\n(ρ≈1)", "inter-\nmediate", "coupling\nloss\n(ρ<1)"],
                    fontsize=4.5)
axC.set_ylabel("% of efficacy-loss variants", fontsize=5.5, labelpad=1)
axC.set_ylim(0, max(vals) * 1.25)
axC.set_title(f"Reclassification of\nefficacy loss (n={len(lofb)})", pad=2)
axC.set_box_aspect(1)
axC.tick_params(pad=1)
for s in ("top", "right"):
    axC.spines[s].set_visible(False)

# row2: efficacy loss top-down | side | colorbar
scE = struct(fig.add_subplot(gs[2, 0]), "eff_loss", "Reds", EFF_MAX, TOPD,
             "Efficacy loss · top-down")
struct(fig.add_subplot(gs[2, 1]), "eff_loss", "Reds", EFF_MAX, SIDE, "Efficacy loss · side")
ceax = fig.add_subplot(gs[2, 2]); ceax.axis("off")
ce = fig.colorbar(scE, ax=ceax, fraction=0.28, pad=0.02, aspect=12)
ce.set_label("fraction of WT E$_{max}$ lost", fontsize=5.5)
ce.ax.tick_params(labelsize=5, width=0.4, length=1.5)

# row3: potency loss top-down | side | colorbar
scP = struct(fig.add_subplot(gs[3, 0]), "pot_loss", "Blues", POT_MAX, TOPD,
             "Potency loss · top-down")
struct(fig.add_subplot(gs[3, 1]), "pot_loss", "Blues", POT_MAX, SIDE, "Potency loss · side")
cpax = fig.add_subplot(gs[3, 2]); cpax.axis("off")
cp = fig.colorbar(scP, ax=cpax, fraction=0.28, pad=0.02, aspect=12)
cp.set_label("rightward EC$_{50}$ shift  (ΔlogEC$_{50}$)", fontsize=5.5)
cp.ax.tick_params(labelsize=5, width=0.4, length=1.5)

fig.savefig(OUT / "fig_supp_operational.pdf")
fig.savefig(OUT / "fig_supp_operational.png", dpi=300)
print("saved", OUT / "fig_supp_operational.pdf")

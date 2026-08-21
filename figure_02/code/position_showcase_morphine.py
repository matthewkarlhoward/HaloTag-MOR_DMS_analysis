#!/usr/bin/env python3
"""
Single-position dose-response showcases for Morphine: raw per-dose data points
plus the robust 3-parameter Hill fit for every high-R^2 sigmoid missense
variant at a chosen position. One figure per position.

Two POSITIONS groups (see below):
  - High-coverage, EC50+Emax dual-variation candidates, auto-labeled from
    structural_master (SSE, GPCRdb generic number, ligand shell) plus coverage.
  - Curated structural / physicochem stories with hand-written role labels
    (L114, I244, N334, orthosteric D149/Y150/M153, ...).
Set a tuple's role to None to auto-derive the label; otherwise the string is used.

Each curve is colored by an ordering variable (`color_by`): mutant-residue
hydrophobicity ('kd'), volume ('vol'), charge ('chg'), or fitted potency
('ec50'). The legend lists each mutation with that variable's value, and the
property->effect Pearson r is annotated on the axes.

A grey band marks the synonymous 2.5-97.5 percentile envelope (the spread of
WT-like fits); variant curves leaving it diverge from synonymous behavior.

Activity scale: points and curves are normalized to the synonymous (WT-like)
reference (mean sigmoid-syn Emin -> 0, Emax -> 1), matching
morphine_all_fits_overlay.py. Set NORMALIZE = False for raw rescaled scores.

Reads:
  dms_scores/composite_dms_scores.csv               (raw per-dose rescaled scores)
  curve_refitting/refit_3param_robust_morphine.csv  (fitted params, class, R^2)
  structures/processed/structural_master.csv        (ligand distance / shell)
Outputs:
  curve_refitting/position_showcase_morphine_<WTpos>.{pdf,png}
"""
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT  = ROOT / "curve_refitting"

R2_MIN    = 0.85
NORMALIZE = True
# Positions where flat (non-sigmoid) missense variants are ALSO drawn, in grey,
# so loss-of-dose-response variants are visible alongside the sigmoid fits.
SHOW_FLAT = set()
# Per-position tuning: colormap override (default viridis) and suppression of
# the property->effect correlation annotation.
CMAP    = {83: "magma", 334: "magma"}
NO_CORR = {334}

# position, color_by (kd|vol|chg|ec50), role label (None -> auto-derive)
POSITIONS = [
    # High-coverage, EC50+Emax dual-variation candidates (auto-labeled, sorted
    # by coverage in the original ranking; review these to pick the exemplar).
    (73,  "ec50", None),
    (122, "ec50", None),
    (110, "ec50", None),
    (100, "ec50", None),
    (85,  "ec50", None),
    (128, "ec50", None),
    (83,  "ec50", None),
    (99,  "ec50", None),
    (185, "ec50", None),
    (202, "ec50", None),
    (90,  "ec50", None),
    (329, "ec50", None),
    (119, "ec50", None),
    (147, "ec50", None),
    # Curated structural / physicochem stories (hand-written labels).
    (114, "ec50", "TM2 (2x48); tunes EC50 + Emax"),
    (244, "kd",   "TM5 (5x48), distal; hydrophobicity -> EC50"),
    (334, "kd",   "NPxxY Asn (7x49); hydrophobicity -> Emax"),
    (148, "ec50", "TM3 (3x31), distal"),
    (149, "ec50", "orthosteric Asp (D3.32) salt bridge"),
    (150, "ec50", "orthosteric Tyr (Y3.33), closest contact"),
    (153, "ec50", "orthosteric (M3.36); tunes EC50 + Emax"),
]

# ── Amino-acid physicochemistry ──────────────────────────────────────────────
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2,
          I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9,
          Y=-1.3, V=4.2)
VOL = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4, G=60.1,
           H=153.2, I=166.7, L=166.7, K=168.6, M=162.9, F=189.9, P=112.7,
           S=89.0, T=116.1, W=227.8, Y=193.6, V=140.0)
CHG = dict(D=-1.0, E=-1.0, K=1.0, R=1.0, H=0.1)
PROP = {"kd": KD, "vol": VOL, "chg": CHG}
PROP_LABEL = {"kd": "hydrophobicity (KD)", "vol": "volume (A^3)", "chg": "charge"}

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

def hill3(x, emin, emax, log_ec50):
    return emin + (emax - emin) / (1 + 10**(log_ec50 - x))

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

try:
    sm = pd.read_csv(ROOT / "structures/processed/structural_master.csv").set_index(
        "residue_number_human")
    dist_lookup   = sm["dist_Morphine_min"].to_dict()
    gpcrdb_lookup = sm["GPCRdb"].to_dict()
    sse_lookup    = sm["SSE"].to_dict()
    shell_lookup  = sm["shell_Morphine"].to_dict()
except Exception:
    dist_lookup = gpcrdb_lookup = sse_lookup = shell_lookup = {}

# Coverage: number of distinct missense substitutions assayed at each position.
n_meas = fit[fit.type == "missense"].groupby("position")["mutation"].nunique().to_dict()

def auto_role(position):
    """SSE + GPCRdb generic number (+ ligand shell) from structural_master."""
    sse   = sse_lookup.get(position, "")
    gp    = gpcrdb_lookup.get(position, "")
    shell = shell_lookup.get(position, "none")
    parts = [f"{sse} ({gp})" if (sse and gp) else (sse or "")]
    if isinstance(shell, str) and shell.endswith("_shell"):
        parts.append(shell.replace("_", " "))
    return "; ".join(p for p in parts if p) or "—"

def norm(y):
    return (y - syn_emin) / (syn_emax - syn_emin) if NORMALIZE else y

xs = np.linspace(-12, -4, 200)
MM = 1 / 25.4

# Synonymous 95% band: 2.5-97.5 percentile envelope across syn sigmoid fits
# (on the Activity scale). Variants leaving this band diverge from WT-like behavior.
_syn_curves = np.array([
    norm(hill3(xs, r.fitted_emin, r.fitted_emax, r.fitted_ec50_logM))
    for r in syn.itertuples(index=False)
])
syn_lo = np.percentile(_syn_curves, 2.5, axis=0)
syn_hi = np.percentile(_syn_curves, 97.5, axis=0)

def plot_position(position, color_by, role):
    sel = fit[(fit.fit_ok) & (fit.position == position)
              & (fit.type == "missense") & (fit.new_curve_type == "sigmoid")
              & (fit.r2 >= R2_MIN)].copy()
    if sel.empty:
        print(f"  position {position}: no qualifying variants, skipped")
        return
    wt = sel.wildtype.iloc[0]
    if role is None:
        role = auto_role(position)

    # Ordering / coloring variable
    if color_by in PROP:
        sel["cvar"] = sel.mutation.map(PROP[color_by])
        if color_by == "chg":
            sel["cvar"] = sel["cvar"].fillna(0.0)
        sel = sel.dropna(subset=["cvar"])
    else:  # potency
        sel["cvar"] = sel.fitted_ec50_logM
    sel = sel.sort_values("cvar")

    # Property -> effect correlations for the annotation
    ann = ""
    if color_by in PROP and position not in NO_CORR and sel["cvar"].nunique() >= 4:
        r_ec, _ = pearsonr(sel["cvar"], sel.fitted_ec50_logM)
        r_em, _ = pearsonr(sel["cvar"], sel.fitted_emax)
        ann = (f"{PROP_LABEL[color_by]} vs\nEC50: r={r_ec:+.2f}\n"
               f"Emax: r={r_em:+.2f}")

    cmap_name = CMAP.get(position, "viridis")
    cmax = 0.88 if cmap_name == "magma" else 1.0
    colors = plt.get_cmap(cmap_name)(np.linspace(0, cmax, len(sel)))

    fig, ax = plt.subplots(figsize=(85 * MM, 55 * MM))
    ax.fill_between(xs, syn_lo, syn_hi, color="grey", alpha=0.25, lw=0,
                    zorder=1, label="syn 95%")
    ax.plot(xs, norm(hill3(xs, syn_emin, syn_emax, syn_ec50)),
            color="black", lw=0.8, ls="--", zorder=2, label="syn mean")

    for color, r in zip(colors, sel.itertuples(index=False)):
        yraw = dms_idx.loc[r.hgvs, dose_cols].to_numpy(dtype=float)
        ax.scatter(doses, norm(yraw), s=7, color=color, edgecolors="none",
                   alpha=0.9, zorder=4)
        if color_by in PROP:
            lab = f"{wt}{position}{r.mutation}  {color_by}={r.cvar:.1f}"
        else:
            lab = f"{wt}{position}{r.mutation}  (R²={r.r2:.2f})"
        ax.plot(xs, norm(hill3(xs, r.fitted_emin, r.fitted_emax, r.fitted_ec50_logM)),
                color=color, lw=1.0, zorder=3, label=lab)

    # Flat (loss-of-dose-response) variants: grey, so the LOF set is visible
    # without competing with the viridis sigmoid fits. 'x' markers + dotted fit.
    n_flat = 0
    if position in SHOW_FLAT:
        flat = fit[(fit.fit_ok) & (fit.position == position)
                   & (fit.type == "missense")
                   & (fit.new_curve_type == "flat")].sort_values("fitted_emax")
        n_flat = len(flat)
        for gc, r in zip(plt.cm.Greys(np.linspace(0.45, 0.85, max(n_flat, 1))),
                         flat.itertuples(index=False)):
            yraw = dms_idx.loc[r.hgvs, dose_cols].to_numpy(dtype=float)
            ax.scatter(doses, norm(yraw), s=5, color=gc, marker="x",
                       linewidths=0.6, alpha=0.7, zorder=2)
            ax.plot(xs, norm(hill3(xs, r.fitted_emin, r.fitted_emax, r.fitted_ec50_logM)),
                    color=gc, lw=0.8, ls=":", zorder=2,
                    label=f"{wt}{position}{r.mutation}  (flat)")

    if ann:
        ax.text(0.03, 0.97, ann, transform=ax.transAxes, va="top", ha="left",
                fontsize=5, color="black")

    dist = dist_lookup.get(position)
    dtxt = f"  |  {dist:.1f} Å to morphine" if isinstance(dist, (int, float)) and np.isfinite(dist) else ""
    nmeas = n_meas.get(position)
    flat_txt = f", {n_flat} flat" if n_flat else ""
    ntxt  = f"  (n={len(sel)}/{nmeas} sig{flat_txt})" if nmeas else f"  (n={len(sel)})"
    ax.set_xlabel("log[Morphine] (M)", labelpad=1)
    ax.set_ylabel("Activity" if NORMALIZE else "rescaled score", labelpad=1)
    ax.set_title(f"{wt}{position} — {role}{dtxt}{ntxt}",
                 color="black", pad=2, fontsize=5.5)
    ax.tick_params(axis="both", pad=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.5)

    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
              handlelength=1.2, handletextpad=0.4, labelspacing=0.22, borderaxespad=0)
    fig.subplots_adjust(left=0.11, bottom=0.15, right=0.72, top=0.90)

    tag = f"{wt}{position}"
    fig.savefig(OUT / f"position_showcase_morphine_{tag}.pdf", dpi=600)
    fig.savefig(OUT / f"position_showcase_morphine_{tag}.png", dpi=600)
    plt.close(fig)
    print(f"  {tag}: n={len(sel)}, color_by={color_by}, "
          f"EC50 {sel.fitted_ec50_logM.min():.2f}..{sel.fitted_ec50_logM.max():.2f}, "
          f"Emax {sel.fitted_emax.min():.2f}..{sel.fitted_emax.max():.2f}"
          + (f", {ann.split(chr(10))[1]} {ann.split(chr(10))[2]}" if ann else ""))

print("Generating per-position showcases:")
for position, color_by, role in POSITIONS:
    plot_position(position, color_by, role)

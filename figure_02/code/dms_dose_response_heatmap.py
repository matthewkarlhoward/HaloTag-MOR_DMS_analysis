#!/usr/bin/env python3
"""
DMS dose-response heatmap summary figure.

Top panel: Morphine / Fentanyl / DAMGO per-dose heatmaps (variant x log[M]),
           plus No Ligand (FSK) and Cell Surface Expression in the DAMGO col.
Bottom panel: per-drug Emin / EC50 / Emax heatmaps.

Outputs:
  plots/heatmaps/dms_dose_response_heatmap.pdf         (combined)
  plots/heatmaps/dms_dose_response_heatmap_top.pdf     (top strip only)
  plots/heatmaps/dms_dose_response_heatmap_bottom.pdf  (bottom strip only)
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import (LinearSegmentedColormap, TwoSlopeNorm,
                               PowerNorm, Normalize)

mpl.rcParams.update({
    "font.family":  "Helvetica",
    "font.weight":  "normal",
    "font.size":    6,
    "axes.titlesize":  6,
    "axes.labelsize":  6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "axes.titleweight": "normal",
    "axes.labelweight": "normal",
    "axes.linewidth":   0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "mathtext.default": "regular",
})

ROOT    = Path(__file__).resolve().parents[2]
CSV     = ROOT / "dms_scores" / "composite_dms_scores.csv"
OUT_DIR = ROOT / "plots" / "heatmaps"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Palettes to render one PDF for each. Surface Expression always uses the
# red/white/blue divergent scale regardless of this choice.
## Keep full magma (darkest end intact) but reshape the ramp so it reaches
## the yellow/bright end faster. A power-warped linspace pulls more of the
## colormap into the upper (bright) portion without dropping any colors.
_magma = plt.get_cmap("magma")
_t = np.linspace(0.0, 1.0, 256) ** 0.55  # gamma < 1 -> faster climb to yellow
CMAP = LinearSegmentedColormap.from_list("magma_fastyellow", _magma(_t))
norm_factory = lambda lim: Normalize(vmin=-lim, vmax=lim)

# --------------------------------------------------------------------------
# Load and prep
# --------------------------------------------------------------------------
df = pd.read_csv(CSV)

indel_order = ["insGSG", "insGS", "insG", "del2", "del1", "del"]
aa_order    = list("PRKHDEFWYNQCSTILVMGA")
row_order   = indel_order + aa_order
n_rows      = len(row_order)

pos_min   = 1
pos_max   = int(df["position"].max())
positions = np.arange(pos_min, pos_max + 1)

def get_mat(col):
    """Pivot to (row_order x positions) matrix, NaN-padded."""
    sub = df[df["mutation"].isin(row_order)]
    p = sub.pivot_table(index="mutation", columns="position",
                        values=col, aggfunc="first")
    return p.reindex(index=row_order, columns=positions).values

# Dose columns, ordered top -> bottom (highest log[M] first)
mor_doses = [
    ("MOR_Morphine_-45M_rescaled_score",  -4.5),
    ("MOR_Morphine_-55M_rescaled_score",  -5.5),
    ("MOR_Morphine_-65M_rescaled_score",  -6.5),
    ("MOR_Morphine_-75M_rescaled_score",  -7.5),
    ("MOR_Morphine_-85M_rescaled_score",  -8.5),
    ("MOR_Morphine_-95M_rescaled_score",  -9.5),
    ("MOR_Morphine_-105M_rescaled_score", -10.5),
    ("MOR_Morphine_-115M_rescaled_score", -11.5),
]
fen_doses = [(f"MOR_Fentanyl_{v}M_rescaled_score", float(v))
             for v in range(-5, -13, -1)]
# DAMGO: dms_1=10uM (-5), dms_2=1uM (-6), ..., dms_5=0.001uM (-9).
# dms_6 is "no ligand" and excluded from the dose-response panel.
dam_doses = [
    ("DAMGO_mor_dms_1_rescaled_score", -5.0),
    ("DAMGO_mor_dms_2_rescaled_score", -6.0),
    ("DAMGO_mor_dms_3_rescaled_score", -7.0),
    ("DAMGO_mor_dms_4_rescaled_score", -8.0),
    ("DAMGO_mor_dms_5_rescaled_score", -9.0),
]

mor_mats = [(get_mat(c), v) for c, v in mor_doses]
fen_mats = [(get_mat(c), v) for c, v in fen_doses]
dam_mats = [(get_mat(c), v) for c, v in dam_doses]

fsk_mat  = get_mat("FSK_effect")
surf_mat = get_mat("Surface_effect")

def ec50_idx_to_logM(drug, idx):
    """Convert fit EC50 from dose-index space back to log[M].
    Assumes the sigmoid fits used x = 1..N in order of ascending concentration.
      Morphine: idx 1..8 -> -11.5..-4.5
      Fentanyl: idx 1..8 -> -12  ..-5
      DAMGO:    idx 1..5 -> -9..-5      (idx 1 = 0.001 uM, idx 5 = 10 uM);
                idx 6    = no-ligand -> NaN
    """
    idx = np.asarray(idx, dtype=float)
    if drug == "Morphine":
        return idx - 12.5
    if drug == "Fentanyl":
        return idx - 13.0
    if drug == "DAMGO":
        out = idx - 10.0                      # 1->-9, 2->-8, ..., 5->-5
        out = np.where(idx >= 5.5, np.nan, out)
        return out
    raise ValueError(drug)

fit = {
    drug: {
        "Emin": get_mat(f"{drug}_emin"),
        "EC50": ec50_idx_to_logM(drug, get_mat(f"{drug}_ec50")),
        "Emax": get_mat(f"{drug}_emax"),
    }
    for drug in ("Morphine", "Fentanyl", "DAMGO")
}

# --------------------------------------------------------------------------
# Color limits
# --------------------------------------------------------------------------
def q98(*arrs):
    v = np.concatenate([a.ravel() for a in arrs])
    return float(np.nanquantile(np.abs(v), 0.98))

def qrange(*arrs, lo=0.02, hi=0.98):
    v = np.concatenate([a.ravel() for a in arrs])
    return (float(np.nanquantile(v, lo)), float(np.nanquantile(v, hi)))

def _round_up(x, step=0.1):
    return float(np.ceil(x / step) * step)
def _round_outer(lo, hi, step=0.5):
    return (float(np.floor(lo / step) * step),
            float(np.ceil(hi / step) * step))

dose_lim = _round_up(q98(*[m for m, _ in mor_mats + fen_mats + dam_mats]), 0.1)
fsk_lim  = _round_up(q98(fsk_mat), 0.1)
surf_lim = _round_up(q98(surf_mat), 0.1)
lim_emin = _round_up(q98(*[fit[d]["Emin"] for d in fit]), 0.1)
lim_emax = _round_up(q98(*[fit[d]["Emax"] for d in fit]), 0.1)
# EC50 is in log[M] (asymmetric, not centered on 0) -> use actual 2-98 pct range
ec50_vmin, ec50_vmax = _round_outer(*qrange(*[fit[d]["EC50"] for d in fit]),
                                    step=0.5)

magma = plt.get_cmap("magma")
rwb   = LinearSegmentedColormap.from_list("rwb", ["#b2182b", "white", "#2166ac"])

# --------------------------------------------------------------------------
# Figure layout
# --------------------------------------------------------------------------
MM = 1 / 25.4  # mm -> inches

def draw(ax, m, cmap, lim, norm_factory, *, show_yrow_labels=False,
         show_xticks=False, show_xlabel=False):
    # Use pcolormesh so the heatmap renders as vector graphics in the PDF
    # (imshow would rasterize at fixed dpi and blur when zoomed).
    x_edges = np.arange(pos_min - 0.5, pos_max + 1.5)
    y_edges = np.arange(-0.5, n_rows + 0.5)
    ax.pcolormesh(
        x_edges, y_edges, m,
        cmap=cmap, norm=norm_factory(lim),
        shading="flat", linewidth=0, rasterized=False,
    )
    ax.set_xlim(pos_min - 0.5, pos_max + 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    if show_yrow_labels:
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels(row_order)
        ax.tick_params(axis="y", length=0, pad=1)
    else:
        ax.set_yticks([])
    if show_xticks:
        ax.set_xticks([pos_min, 100, 200, 300, pos_max])
        ax.tick_params(axis="x", length=2, pad=1)
        if show_xlabel:
            ax.set_xlabel("Position", labelpad=1)
    else:
        ax.set_xticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.5)

def put_right(ax, text):
    ax.text(1.018, 0.5, text, transform=ax.transAxes,
            ha="left", va="center")

def put_left(ax, text):
    ax.text(-0.012, 0.5, text, transform=ax.transAxes,
            ha="right", va="center")

# --------------------------------------------------------------------------
# Legend helpers (used by all three figure builders)
# --------------------------------------------------------------------------
def _fmt_sym(lim):
    return [-lim, 0.0, lim]

def _fmt_lin(lo, hi):
    return [lo, (lo + hi) / 2, hi]

def _fmt_tick(t, strip_sign):
    v = abs(t) if strip_sign else t
    if abs(v) < 1e-9:
        return "0"
    return f"{v:.1f}"

def _render_legend(fig, entries, cb_w, cb_h):
    """Each entry is a tuple:
       (left, bottom, cmap, norm, title, ticks, strip_sign, (l_anchor, r_anchor))
    """
    for left, bottom, cm, nm, title, ticks, strip_sign, anchors in entries:
        cax = fig.add_axes([left, bottom, cb_w, cb_h])
        cb = mpl.colorbar.ColorbarBase(
            cax, cmap=cm, norm=nm, orientation="horizontal",
        )
        cb.set_ticks(ticks)
        cb.ax.set_xticklabels([_fmt_tick(t, strip_sign) for t in ticks])
        cb.ax.tick_params(labelsize=6, length=2, pad=1, width=0.5)
        cb.outline.set_linewidth(0.5)
        cax.set_title(title, fontsize=6, pad=4)
        l_anchor, r_anchor = anchors
        cax.text(-0.03, 0.5, l_anchor, transform=cax.transAxes,
                 fontsize=6, ha="right", va="center")
        cax.text(1.03, 0.5, r_anchor, transform=cax.transAxes,
                 fontsize=6, ha="left", va="center")


def _render_top_block(fig, gs, cmap, norm_factory, row_offset=0):
    """Render the 8-row dose heatmap block into rows [row_offset, row_offset+8) of gs."""
    top_first_row_axes = [None, None, None]
    for r in range(8):
        is_bottom = (r == 7)

        m, lab = mor_mats[r]
        ax = fig.add_subplot(gs[row_offset + r, 0])
        draw(ax, m, cmap, dose_lim, norm_factory,
             show_xticks=is_bottom, show_xlabel=is_bottom)
        put_right(ax, f"{abs(lab):g}")
        if r == 0:
            top_first_row_axes[0] = ax

        m, lab = fen_mats[r]
        ax = fig.add_subplot(gs[row_offset + r, 1])
        draw(ax, m, cmap, dose_lim, norm_factory,
             show_xticks=is_bottom, show_xlabel=is_bottom)
        put_right(ax, f"{abs(lab):g}")
        if r == 0:
            top_first_row_axes[1] = ax

        ax = fig.add_subplot(gs[row_offset + r, 2])
        if r < 5:
            m, lab = dam_mats[r]
            draw(ax, m, cmap, dose_lim, norm_factory)
            put_right(ax, f"{abs(lab):g}")
        elif r == 5:
            draw(ax, fsk_mat, cmap, fsk_lim, norm_factory)
            ax.set_title("No Ligand", pad=1)
            bb = ax.get_position()
            ax.set_position([bb.x0, bb.y0 - bb.height / 2,
                             bb.width, bb.height])
        elif r == 6:
            ax.axis("off")
        else:
            draw(ax, surf_mat, rwb, surf_lim, norm_factory,
                 show_xticks=True, show_xlabel=True)
            ax.set_title("Cell Surface Expression", pad=1)
        if r == 0:
            top_first_row_axes[2] = ax

    for ax, title in zip(top_first_row_axes,
                         ["Morphine", "Fentanyl", "DAMGO"]):
        ax.set_title(title, pad=4)


def _render_bottom_block(fig, gs, cmap, norm_factory, row_offset=0):
    """Render the 3-row Emin/EC50/Emax block into rows [row_offset, row_offset+3) of gs."""
    ec50_norm = lambda _lim: Normalize(vmin=ec50_vmin, vmax=ec50_vmax)
    metrics = [
        ("Emin", lim_emin,             norm_factory),
        ("EC50", max(abs(ec50_vmin),
                     abs(ec50_vmax)),  ec50_norm),
        ("Emax", lim_emax,             norm_factory),
    ]
    drugs = ["Morphine", "Fentanyl", "DAMGO"]
    for i, (metric, mlim, nf) in enumerate(metrics):
        for j, drug in enumerate(drugs):
            ax = fig.add_subplot(gs[row_offset + i, j])
            draw(ax, fit[drug][metric], cmap, mlim, nf,
                 show_xticks=(i == 2), show_xlabel=(i == 2))
            if i == 0:
                ax.set_title(drug, pad=4)


def build_figure(cmap, norm_factory):
    """Combined figure: top + bottom blocks with 3-colorbar legend."""
    fig = plt.figure(figsize=(155 * MM, 95 * MM))
    gs = fig.add_gridspec(
        nrows=12, ncols=3,
        height_ratios=[1] * 8 + [1.8] + [1] * 3,
        hspace=0.06, wspace=0.12,
        left=0.012, right=0.965, top=0.968, bottom=0.205,
    )

    _render_top_block(fig, gs, cmap, norm_factory, row_offset=0)
    _render_bottom_block(fig, gs, cmap, norm_factory, row_offset=9)

    cb_w, cb_h = 0.16, 0.012
    col_x = [0.078, 0.420, 0.762]
    y_legend = 0.060
    entries = [
        (col_x[0], y_legend, cmap, norm_factory(dose_lim),
         "Signaling score",
         _fmt_sym(dose_lim), False,
         (r"$\uparrow$signaling", r"$\downarrow$signaling")),
        (col_x[1], y_legend, cmap, Normalize(vmin=ec50_vmin, vmax=ec50_vmax),
         r"$EC_{50}$ (log[ligand], M)",
         _fmt_lin(ec50_vmin, ec50_vmax), True,
         (r"$\uparrow$potency", r"$\downarrow$potency")),
        (col_x[2], y_legend, rwb, norm_factory(surf_lim),
         "Cell surface expression",
         _fmt_sym(surf_lim), False,
         (r"$\downarrow$surface", r"$\uparrow$surface")),
    ]
    _render_legend(fig, entries, cb_w, cb_h)
    return fig


def build_top_figure(cmap, norm_factory):
    """Top strip only: dose heatmaps + No Ligand + Cell Surface Expression.
       Legend: Signaling score + Cell surface expression."""
    fig = plt.figure(figsize=(155 * MM, 70 * MM))
    gs = fig.add_gridspec(
        nrows=8, ncols=3,
        hspace=0.06, wspace=0.12,
        left=0.012, right=0.965, top=0.955, bottom=0.305,
    )

    _render_top_block(fig, gs, cmap, norm_factory, row_offset=0)

    cb_w, cb_h = 0.16, 0.016
    col_x = [0.22, 0.62]
    y_legend = 0.075
    entries = [
        (col_x[0], y_legend, cmap, norm_factory(dose_lim),
         "Signaling score",
         _fmt_sym(dose_lim), False,
         (r"$\uparrow$signaling", r"$\downarrow$signaling")),
        (col_x[1], y_legend, rwb, norm_factory(surf_lim),
         "Cell surface expression",
         _fmt_sym(surf_lim), False,
         (r"$\downarrow$surface", r"$\uparrow$surface")),
    ]
    _render_legend(fig, entries, cb_w, cb_h)
    return fig


def build_bottom_figure(cmap, norm_factory):
    """Bottom strip only: Emin / EC50 / Emax heatmaps.
       Legend: Signaling score + EC50."""
    fig = plt.figure(figsize=(155 * MM, 42 * MM))
    gs = fig.add_gridspec(
        nrows=3, ncols=3,
        hspace=0.06, wspace=0.12,
        left=0.012, right=0.965, top=0.86, bottom=0.50,
    )

    _render_bottom_block(fig, gs, cmap, norm_factory, row_offset=0)

    cb_w, cb_h = 0.16, 0.026
    col_x = [0.22, 0.62]
    y_legend = 0.135
    entries = [
        (col_x[0], y_legend, cmap, norm_factory(lim_emax),
         "Signaling score",
         _fmt_sym(lim_emax), False,
         (r"$\uparrow$signaling", r"$\downarrow$signaling")),
        (col_x[1], y_legend, cmap, Normalize(vmin=ec50_vmin, vmax=ec50_vmax),
         r"$EC_{50}$ (log[ligand], M)",
         _fmt_lin(ec50_vmin, ec50_vmax), True,
         (r"$\uparrow$potency", r"$\downarrow$potency")),
    ]
    _render_legend(fig, entries, cb_w, cb_h)
    return fig


outputs = [
    (build_figure,        "dms_dose_response_heatmap.pdf"),
    (build_top_figure,    "dms_dose_response_heatmap_top.pdf"),
    (build_bottom_figure, "dms_dose_response_heatmap_bottom.pdf"),
]
for builder, fname in outputs:
    fig = builder(CMAP, norm_factory)
    out = OUT_DIR / fname
    fig.savefig(out, dpi=600)
    plt.close(fig)
    print(f"Saved: {out}")

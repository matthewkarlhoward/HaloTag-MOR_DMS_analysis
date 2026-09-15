"""House style + shared plotting helpers.

House style: black text, Helvetica, 6 pt, 0.5 pt lines.
Force the Agg backend so figsize is honoured exactly (the macosx backend snaps
figure size to 2 dp).
"""
from __future__ import annotations

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

MM = 1 / 25.4

RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 6,
    "axes.labelsize": 6,
    "axes.titlesize": 6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "text.color": "black",
    "axes.labelcolor": "black",
    "axes.edgecolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2,
    "ytick.major.size": 2,
    "lines.linewidth": 0.5,
    "patch.linewidth": 0.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.transparent": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
}
mpl.rcParams.update(RC)

# ligand efficacy classes, colours from plots/radars/radars.Rmd
CLASS_COLORS = {
    "strong": "#231f20",
    "intermediate": "#9b1c1d",
    "weak": "#3779b9",
    "antagonist": "#6f08a3",
    "no_ligand": "#8a8a8a",
}

# microswitch element colours (Phase 2/3 sign-test panels)
ELEMENT_COLORS = {
    "DRY": "#e08214",
    "CWxP": "#1b7837",
    "NPxxY": "#c51b7d",
    "PIF": "#4d4d4d",
    "Na_site": "#0571b0",
    "TM6_lock": "#ca0020",
    "orthosteric": "#7b3294",
    "none": "#cccccc",
}


def figure(w_mm: float, h_mm: float):
    """Figure with an exact mm canvas (see reference_matplotlib_exact_mm_figsize)."""
    fig = plt.figure()
    fig.set_size_inches(w_mm * MM, h_mm * MM)
    return fig


def save(fig, path, dpi=600):
    fig.savefig(path, dpi=dpi, facecolor="white")
    print(f"  wrote {path}")
    plt.close(fig)

#!/usr/bin/env python3
"""
Dose-response QC composite. 2 rows x 3 cols, aligned columns per drug.

Row 1: Score per dose (strip plot)
Row 2: SD per dose (strip plot)
Cols: Fentanyl / Morphine / DAMGO (DAMGO no-lig on the left)

Output: DMS_QC/04_dose_response_composite.pdf
"""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

mpl.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 6,
    "text.color": "black",
    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "axes.titlesize": 6,
    "axes.labelsize": 6,
    "xtick.color": "black", "ytick.color": "black",
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
})

HERE = Path(__file__).resolve().parent
CSV = HERE.parent / "dms_scores" / "composite_dms_scores.csv"
df = pd.read_csv(CSV, low_memory=False)

# (prefix, dose_x, label, is_no_lig)
FENTANYL = [
    (f"MOR_Fentanyl_{v}M", float(v), f"{v}", False)
    for v in ["-12", "-11", "-10", "-9", "-8", "-7", "-6", "-5"]
]
MORPHINE = [
    (f"MOR_Morphine_{v}M", float(v) / 10, f"{float(v)/10:g}", False)
    for v in ["-115", "-105", "-95", "-85", "-75", "-65", "-55", "-45"]
]
DAMGO_MAP = {5: -9, 4: -8, 3: -7, 2: -6, 1: -5}
DAMGO = (
    [("DAMGO_mor_dms_6", None, "no lig.", True)]
    + [(f"DAMGO_mor_dms_{i}", DAMGO_MAP[i], f"{DAMGO_MAP[i]}", False) for i in [5, 4, 3, 2, 1]]
)

SERIES = [("Morphine", MORPHINE),
          ("Fentanyl", FENTANYL),
          ("DAMGO",    DAMGO)]

GREY = "#888888"
GREY_DARK = "#222222"

def fetch(prefix, col):
    a = df[f"{prefix}_{col}"].to_numpy()
    return a[~np.isnan(a)]

def _violin_jitter(d, max_width, rng):
    """Return x-offsets shaped like a violin: each point's jitter range scaled
    by local KDE density so the cloud is fat where the distribution is dense."""
    if len(d) < 2 or np.allclose(d.std(), 0):
        return np.zeros(len(d))
    kde = gaussian_kde(d)
    dens = kde(d)
    norm = dens / dens.max()
    return rng.uniform(-1, 1, size=len(d)) * norm * max_width

SUBSAMPLE_N = 800  # per condition; keeps KDE shape but avoids opaque blobs

def strip_panel(ax, data_list, positions, no_lig_idx=None,
                ylim=None, ylabel=None, jitter=0.38, rng_seed=0,
                point_size=0.7, point_alpha=0.55):
    rng = np.random.default_rng(rng_seed)
    for i, (x, d) in enumerate(zip(positions, data_list)):
        # subsample for visibility while preserving distribution shape
        if len(d) > SUBSAMPLE_N:
            d_plot = rng.choice(d, size=SUBSAMPLE_N, replace=False)
        else:
            d_plot = d
        c = GREY_DARK if (no_lig_idx is not None and i == no_lig_idx) else GREY
        xs = x + _violin_jitter(d_plot, jitter, rng)
        ax.scatter(xs, d_plot, s=point_size, c=c, alpha=point_alpha,
                   linewidths=0, rasterized=True)
    meds = [np.median(d) for d in data_list]
    ax.plot(positions, meds, color="black", lw=0.8, marker="o", ms=2.4, zorder=4)
    for x, d in zip(positions, data_list):
        q25, q75 = np.percentile(d, [25, 75])
        ax.plot([x, x], [q25, q75], color="black", lw=0.7, zorder=3)
    ax.axhline(0, color="black", lw=0.3, ls="--", alpha=0.4)
    if ylim is not None:
        ax.set_ylim(*ylim)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

# ---- Layout: 160 x 50 mm final size ----
MM = 1.0 / 25.4
fig = plt.figure(figsize=(160 * MM, 50 * MM))
gs = fig.add_gridspec(
    2, 3,
    width_ratios=[len(s[1]) for s in SERIES],
    height_ratios=[1, 1],
    wspace=0.10, hspace=0.30,
    left=0.07, right=0.99, top=0.90, bottom=0.18,
)

SCORE_LIM = (-1.25, 0.75)
SCORE_TICKS = [-1.0, -0.5, 0.0, 0.5]
SD_LIM = (0.0, 0.4)

for c, (drug, series) in enumerate(SERIES):
    scores, sds, positions, ticklabels, no_lig_idx = [], [], [], [], None
    for i, (prefix, _, label, is_nl) in enumerate(series):
        scores.append(fetch(prefix, "rescaled_score"))
        sds.append(fetch(prefix, "rescaled_sd"))
        positions.append(i); ticklabels.append(label)
        if is_nl:
            no_lig_idx = i
    xlim = (-0.6, max(positions) + 0.6)

    # Row 1: scores
    ax_top = fig.add_subplot(gs[0, c])
    strip_panel(ax_top, scores, positions, no_lig_idx=no_lig_idx,
                ylim=SCORE_LIM, ylabel=("Score" if c == 0 else None))
    ax_top.set_yticks(SCORE_TICKS)
    if c > 0:
        ax_top.set_yticklabels([])
    ax_top.set_xticks(positions); ax_top.set_xticklabels([])
    ax_top.set_xlim(*xlim)
    ax_top.set_title(drug, loc="center", pad=4)

    # Row 2: SDs (aligned x-axis to row 1)
    ax_bot = fig.add_subplot(gs[1, c], sharex=ax_top)
    strip_panel(ax_bot, sds, positions, no_lig_idx=no_lig_idx,
                ylim=SD_LIM, ylabel=("SD" if c == 0 else None))
    if c > 0:
        ax_bot.set_yticklabels([])
    ax_bot.set_xticks(positions); ax_bot.set_xticklabels(ticklabels)
    ax_bot.set_xlim(*xlim)
    ax_bot.set_xlabel("log[ligand] (M)")

out = HERE / "04_dose_response_composite.pdf"
fig.savefig(out, dpi=600)  # no bbox_inches="tight" so the page stays 80x50 mm
print(f"Wrote {out}")

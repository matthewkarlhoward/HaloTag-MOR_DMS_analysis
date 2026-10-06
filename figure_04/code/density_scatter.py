#!/usr/bin/env python3
"""
Shared point style for the variant-level Morphine vs Fentanyl scatters
(panels b and d of the chemotype figure).

Three interchangeable looks, picked with the SCATTER_STYLE environment variable:

  magma  density-coloured (gaussian KDE, magma, dense drawn on top) -- the
         original panel
  grey   the same KDE density on a light-grey -> black ramp, so the cloud reads
         as density without carrying a second colour meaning
  alpha  plain black points at low opacity; density shows as accumulated ink

`grey` and `alpha` both leave colour free for the overlays (green Q126 in panel
b, blue/red bias elsewhere), which is why they were asked for. `alpha` is the
honest one for sparse regions -- a single variant stays a single faint dot
rather than being pulled onto a colour ramp -- while `grey` keeps the dense core
legible when the cloud saturates.
"""
import os

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde

# Greys truncated at the light end: the full colormap starts at white, which
# makes the sparse outliers invisible on a white page.
GREY_CMAP = LinearSegmentedColormap.from_list(
    "density_greys", ["#d0d0d0", "#8a8a8a", "#3a3a3a", "#000000"])

STYLE = os.environ.get("SCATTER_STYLE", "alpha")
# Ink per point for the "alpha" style: low enough that the core still grades,
# high enough that a lone variant in the sparse tail stays visible at 40 mm.
ALPHA = float(os.environ.get("SCATTER_ALPHA", 0.18))
SUFFIX = {"magma": "", "grey": "_grey", "alpha": "_alpha"}


def style_suffix(style=None):
    """Filename suffix for a style ('' for the original magma panel)."""
    return SUFFIX[style or STYLE]


def draw_density(ax, x, y, style=None, s=1.0, zorder=2, alpha=None):
    """Draw the background cloud of one scatter panel."""
    style = style or STYLE
    x, y = np.asarray(x, float), np.asarray(y, float)
    if style == "alpha":
        return ax.scatter(x, y, s=s, c="black", alpha=ALPHA if alpha is None else alpha,
                          edgecolors="none", zorder=zorder)
    z = gaussian_kde(np.vstack([x, y]))(np.vstack([x, y]))
    o = z.argsort()                      # dense points last = on top
    cmap = "magma" if style == "magma" else GREY_CMAP
    return ax.scatter(x[o], y[o], c=z[o], s=s, cmap=cmap,
                      edgecolors="none", alpha=0.85, zorder=zorder)

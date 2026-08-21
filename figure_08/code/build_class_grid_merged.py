#!/usr/bin/env python3
"""
Stitch the 10 per-class PNGs into a 2-row grid, for the *merged* grouping
where the old strong + full classes are pooled into "Strong Agonists":
  row 0 = LoF  (top)
  row 1 = GoF  (bottom)
  columns ordered: no_ligand -> antagonist -> weak -> intermediate -> strong_agonists

Prereq: run chimerax_scripts/fig_class_grid_merged_darkgrey.cxc in ChimeraX so
the PNGs exist in network_tools/chimerax_renders/combined_series_merged/.

Output: network_tools/chimerax_renders/class_grid_merged.png
"""
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

ROOT = Path(__file__).resolve().parent
PNG_DIR = ROOT / "chimerax_renders" / "combined_series_merged"
OUT = ROOT / "chimerax_renders" / "class_grid_merged.png"

CLASS_ORDER = ["no_ligand", "antagonist", "weak", "intermediate", "strong_agonists"]
DIRECTIONS = [("lof", "LoF"), ("gof", "GoF")]


def main():
    if not PNG_DIR.exists():
        sys.exit(f"No PNG directory at {PNG_DIR}. Run the slideshow in ChimeraX first.")

    n_rows = len(DIRECTIONS)
    n_cols = len(CLASS_ORDER)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * 2.4, n_rows * 2.4 + 0.4),
        gridspec_kw=dict(wspace=0.02, hspace=0.05),
        squeeze=False,
    )

    for j, cls in enumerate(CLASS_ORDER):
        axes[0, j].set_title(cls.replace("_", " ").title(), fontsize=12, pad=4)

    for i, (dir_key, dir_label) in enumerate(DIRECTIONS):
        for j, cls in enumerate(CLASS_ORDER):
            ax = axes[i, j]
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            png = PNG_DIR / f"class_{dir_key}_{cls}_darkgrey.png"
            if png.exists():
                ax.imshow(mpimg.imread(png))
            else:
                ax.text(0.5, 0.5, "missing", ha="center", va="center",
                        fontsize=8, color="red", transform=ax.transAxes)
            if j == 0:
                ax.set_ylabel(dir_label, fontsize=14, rotation=0,
                              labelpad=30, ha="right", va="center")

    plt.subplots_adjust(left=0.05, right=0.995, top=0.93, bottom=0.02)
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()

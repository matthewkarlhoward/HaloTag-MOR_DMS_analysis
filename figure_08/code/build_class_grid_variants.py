#!/usr/bin/env python3
"""
Stitch the per-class PNGs of each display variant (nolabel / scaled /
scaled_labeled) into the same 2-row grid as build_class_grid_merged.py:
  row 0 = LoF, row 1 = GoF
  cols  = no_ligand -> antagonist -> weak -> intermediate -> strong_agonists

Prereq: render chimerax_scripts/fig_class_variants_all.cxc in ChimeraX first.
Output: chimerax_renders/class_grid_merged_<variant>.png
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

ROOT = Path(__file__).resolve().parent
RENDER_ROOT = ROOT / "chimerax_renders"
VARIANTS = ["nolabel", "scaled", "scaled_labeled"]
CLASS_ORDER = ["no_ligand", "antagonist", "weak", "intermediate", "strong_agonists"]
DIRECTIONS = [("lof", "LoF"), ("gof", "GoF")]


def build(variant):
    png_dir = RENDER_ROOT / f"combined_series_merged_{variant}"
    out = RENDER_ROOT / f"class_grid_merged_{variant}.png"
    n_rows, n_cols = len(DIRECTIONS), len(CLASS_ORDER)
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
            png = png_dir / f"class_{dir_key}_{cls}_darkgrey.png"
            if png.exists():
                ax.imshow(mpimg.imread(png))
            else:
                ax.text(0.5, 0.5, "missing", ha="center", va="center",
                        fontsize=8, color="red", transform=ax.transAxes)
            if j == 0:
                ax.set_ylabel(dir_label, fontsize=14, rotation=0,
                              labelpad=30, ha="right", va="center")
    plt.subplots_adjust(left=0.05, right=0.995, top=0.93, bottom=0.02)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out.relative_to(ROOT)}")


if __name__ == "__main__":
    for v in VARIANTS:
        build(v)

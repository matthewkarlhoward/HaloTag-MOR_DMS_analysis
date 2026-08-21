#!/usr/bin/env python3
"""
build_class_variant_assets.py
----------------------------------------------------------------------------
Derive display variants of the merged per-class ChimeraX network panels that
build_merged_class_assets.py already produced in chimerax_scripts/classes_merged/.

Three variant sets are written (nothing existing is overwritten — each lands
in its own *_<variant> sibling dir, and every variant .cxc still references the
SAME shared .defattr/.pb/.cif assets by absolute path):

    nolabel         identical panels, but the per-residue number labels removed
    scaled          nolabel + node spheres sized by the LoF/GoF disruption score
    scaled_labeled  scaled, with the residue-number labels kept

Sphere sizing (the "scaled" sets):
  * radius is mapped from the per-position class count (the same value already
    loaded as the {direction}_{cls}_count residue attribute) onto [RMIN, RMAX]
  * AREA-proportional:  radius proportional to sqrt(count)  (range-normalized)
  * PER-DIRECTION scope: LoF panels share one count range, GoF panels another,
    so each row fills the full size range but a LoF sphere and a GoF sphere of
    equal radius do NOT encode the same count (LoF tops out ~17.7, GoF ~8).

A shared view treatment is applied to every variant panel (constants below):
the 2D title is dropped, the structure is rotated TURN_Y deg about the vertical
screen axis, recentered+fit on the receptor (FIT_SPEC) then scaled to fill the
frame while leaving margin (ZOOM), and saved on a transparent background. Only
the 3D residue-number labels distinguish nolabel from scaled_labeled.

Outputs:
    chimerax_scripts/classes_merged_<variant>/class_{dir}_{cls}_darkgrey.cxc
    chimerax_scripts/fig_class_grid_merged_darkgrey_<variant>.cxc   (per-variant render)
    chimerax_scripts/fig_class_variants_all.cxc                     (all 3, one session)
"""

import math
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent                 # network_tools/
ATTR_DIR = ROOT / "chimerax_attributes"
CXC_DIR = ROOT / "chimerax_scripts"
SRC_DIR = CXC_DIR / "classes_merged"                   # canonical (labeled, fixed-size)
RENDER_ROOT = ROOT / "chimerax_renders"

# ── sphere-size mapping (scaled variants) ────────────────────────────────────
RMIN, RMAX = 0.55, 2.40        # Angstrom atomRadius (canonical fixed size is 1.1)
RENDER_W, RENDER_H, SUPERSAMPLE = 1600, 1200, 3

# ── shared view treatment applied to every variant panel ─────────────────────
SHOW_TITLE = False     # drop the 2D panel title (the `2dlabels …` line)
FIT_SPEC = "#1/R"      # recenter+fit on the receptor only. Plain `view` fits the bounding
                       #   sphere of ALL models incl. the hidden G-protein, shrinking the
                       #   receptor to ~half size; fitting #1/R centers chain R reliably.
ZOOM = 1.6             # scale vs the receptor fit (1.0 = chain R just fits). Fills most of the
                       #   frame but leaves margin so the top/bottom never clip (1.10 clipped).
TURN_Y = -5            # deg about the vertical screen axis; negative = clockwise from top
TRANSPARENT_BG = True  # save on a transparent ("clear") background instead of white

DIRECTIONS = ["lof", "gof"]
# (labels?, scale?) for each variant
VARIANTS = {
    "nolabel":        dict(labels=False, scale=False),
    "scaled":         dict(labels=False, scale=True),
    "scaled_labeled": dict(labels=True,  scale=True),
}

# line matchers against the canonical .cxc
RE_SEL   = re.compile(r"^sel #2/R:([\d,]+)@CA", re.M)
RE_SIZE  = re.compile(r"^size sel atomRadius")
RE_LABEL = re.compile(r"^label #2/R:.*residues")
RE_2DLAB = re.compile(r"^2dlabels\b")
RE_ZOOM  = re.compile(r"^zoom\b")


def parse_defattr(path):
    d = {}
    for line in open(path):
        m = re.match(r"\s*/R:(\d+)\s+([-\d.]+)", line)
        if m:
            d[int(m.group(1))] = float(m.group(2))
    return d


def parse_nodes(text):
    m = RE_SEL.search(text)
    return [int(x) for x in m.group(1).split(",")] if m else []


def panel_name(cxc_path):
    """class_lof_strong_agonists_darkgrey -> ('lof', 'strong_agonists')."""
    m = re.match(r"class_(lof|gof)_(.+)_darkgrey", cxc_path.stem)
    return m.group(1), m.group(2)


def radius_for(count, cmin, cmax):
    """Range-normalized, area-proportional: radius proportional to sqrt(count)."""
    if cmax <= cmin:
        return round((RMIN + RMAX) / 2, 3)
    t = (math.sqrt(count) - math.sqrt(cmin)) / (math.sqrt(cmax) - math.sqrt(cmin))
    t = max(0.0, min(1.0, t))
    return round(RMIN + (RMAX - RMIN) * t, 3)


def size_lines(nodes, counts, cmin, cmax):
    """Grouped `size` commands: one line per distinct rounded radius."""
    groups = defaultdict(list)
    for p in nodes:
        groups[radius_for(counts.get(p, 0.0), cmin, cmax)].append(p)
    out = []
    for r in sorted(groups):
        spec = ",".join(str(p) for p in sorted(groups[r]))
        out.append(f"size #2/R:{spec}@CA atomRadius {r:g}")
    return out


# ── 1. per-direction count range over all node positions ─────────────────────
src_files = sorted(SRC_DIR.glob("class_*_darkgrey.cxc"))
dir_counts = {d: [] for d in DIRECTIONS}
panels = []   # (path, direction, cls, nodes, counts)
for cxc in src_files:
    direction, cls = panel_name(cxc)
    text = cxc.read_text()
    nodes = parse_nodes(text)
    counts = parse_defattr(ATTR_DIR / f"{direction}_{cls}_count.defattr")
    panels.append((cxc, direction, cls, nodes, counts))
    dir_counts[direction] += [counts.get(p, 0.0) for p in nodes]

dir_range = {d: (min(v), max(v)) for d, v in dir_counts.items() if v}
print("Per-direction node-count range (anchors the sphere size scale):")
for d in DIRECTIONS:
    if d in dir_range:
        lo, hi = dir_range[d]
        print(f"  {d}: count [{lo:.3f}, {hi:.3f}] -> radius [{RMIN}, {RMAX}] A")


# ── 2. write the variant .cxc files ──────────────────────────────────────────
for variant, opt in VARIANTS.items():
    out_dir = CXC_DIR / f"classes_merged_{variant}"
    out_dir.mkdir(parents=True, exist_ok=True)
    for cxc, direction, cls, nodes, counts in panels:
        cmin, cmax = dir_range.get(direction, (0.0, 0.0))
        out_lines = []
        for line in cxc.read_text().splitlines():
            if not opt["labels"] and RE_LABEL.match(line):
                continue                                   # drop residue-number labels
            if not SHOW_TITLE and RE_2DLAB.match(line):
                continue                                   # drop the 2D panel title
            if opt["scale"] and RE_SIZE.match(line) and nodes:
                out_lines += size_lines(nodes, counts, cmin, cmax)
                continue                                   # replace fixed size with scaled
            if RE_ZOOM.match(line):                        # came right after `view matrix camera`
                if TURN_Y:
                    out_lines.append(f"turn y {TURN_Y}")   # rotate about the vertical screen axis
                out_lines.append(f"view {FIT_SPEC}")       # recenter+fit the receptor
                out_lines.append(f"zoom {ZOOM}")           # then scale up to fill, leaving margin
                continue
            out_lines.append(line)
        (out_dir / cxc.name).write_text("\n".join(out_lines) + "\n")
    print(f"wrote {len(panels)} panels -> {out_dir.relative_to(ROOT)}")


# ── 3. render slideshows (open each panel, save PNG) ─────────────────────────
def render_block(variant):
    out_render = RENDER_ROOT / f"combined_series_merged_{variant}"
    out_render.mkdir(parents=True, exist_ok=True)
    scr_dir = CXC_DIR / f"classes_merged_{variant}"
    lines = []
    # keep the canonical column/row ordering used by the grid stitcher
    order = ["no_ligand", "antagonist", "weak", "intermediate", "strong_agonists"]
    for direction in DIRECTIONS:
        for cls in order:
            panel = scr_dir / f"class_{direction}_{cls}_darkgrey.cxc"
            if not panel.exists():
                continue
            png = out_render / f"class_{direction}_{cls}_darkgrey.png"
            save_opts = (f"width {RENDER_W} height {RENDER_H} supersample {SUPERSAMPLE}"
                         + (" transparentBackground true" if TRANSPARENT_BG else ""))
            lines.append(f"open {panel}")
            lines.append("wait 20")
            lines.append(f"save {png} {save_opts}")
    return lines


for variant in VARIANTS:
    body = [f"# Render the '{variant}' merged-class variant panels",
            f"# (auto-generated by {Path(__file__).name})",
            "# windowsize matches the save aspect so `zoom`/2dlabels frame the",
            "# same way every launch (ChimeraX otherwise reuses remembered window",
            "# geometry, which clips the 2D title and shifts the framing).",
            "",
            f"windowsize {RENDER_W} {RENDER_H}",
            ""]
    body += render_block(variant)
    (CXC_DIR / f"fig_class_grid_merged_darkgrey_{variant}.cxc").write_text(
        "\n".join(body) + "\n")

# one combined session that renders all three variants then quits (one window)
combined = ["# Render ALL merged-class display variants in one ChimeraX session,",
            "# then exit.  (auto-generated by build_class_variant_assets.py)",
            "# windowsize matches the save aspect (see per-variant scripts for why).",
            "",
            f"windowsize {RENDER_W} {RENDER_H}",
            ""]
for variant in VARIANTS:
    combined.append(f"# ---- {variant} ----")
    combined += render_block(variant)
    combined.append("")
combined.append("exit")
(CXC_DIR / "fig_class_variants_all.cxc").write_text("\n".join(combined) + "\n")

print("\nRender scripts written:")
for variant in VARIANTS:
    print(f"  chimerax_scripts/fig_class_grid_merged_darkgrey_{variant}.cxc")
print("  chimerax_scripts/fig_class_variants_all.cxc  (all three + exit)")

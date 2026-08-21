#!/usr/bin/env python3
"""
Per-class node/edge summaries for the functional-residue networks (Fig. panels A/B).

Node sets and edges are read DIRECTLY from the rendered ChimeraX assets in
chimerax_scripts/classes_merged/ (nodes: the `sel #2/R:...@CA` line) and
chimerax_pseudobonds/classes_merged/ (edges: the Ca-Ca pseudobond pairs), so every
count and coordinate here matches the published panels exactly.

Outputs (network_tools/class_network_plots/):
  class_node_counts_bar.pdf     # nodes per class, LoF vs GoF
  class_edge_counts_bar.pdf     # edges per class, LoF vs GoF
  class_node_density_lof.pdf    # 4 classes x 3 projections, node KDE   (LoF)
  class_node_density_gof.pdf    #   "                                    (GoF)
  class_edge_density_lof.pdf    # 4 classes x 3 projections, edge-midpoint KDE (LoF)
  class_edge_density_gof.pdf    #   "                                    (GoF)

Coordinates: 8EFQ chain-R Ca, PCA-aligned to the receptor principal axes so the
projections are interpretable — z = long axis (membrane normal), x/y = membrane
plane. Panels: x-y = top-down (membrane plane); x-z and y-z = side views.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

ROOT = Path(__file__).resolve().parent
CXC = ROOT / "chimerax_scripts" / "classes_merged"
PB = ROOT / "chimerax_pseudobonds" / "classes_merged"
CA_CSV = ROOT / "pdb_8efq_chain_R_ca.csv"
OUT = ROOT / "class_network_plots"
OUT.mkdir(exist_ok=True)

CLASSES = ["antagonist", "weak", "intermediate", "strong_agonists"]
TITLES = {"antagonist": "Antagonist", "weak": "Weak", "intermediate": "Intermediate",
          "strong_agonists": "Strong"}
COLORS = {"antagonist": "#6f08a3", "weak": "#3779b9",
          "intermediate": "#9b1c1d", "strong_agonists": "#231f20"}
DIRECTIONS = ["lof", "gof"]
DIR_LABEL = {"lof": "LoF", "gof": "GoF"}

# ── house style ───────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 6, "mathtext.default": "regular",
    "text.color": "black", "axes.labelcolor": "black", "axes.edgecolor": "black",
    "axes.titlesize": 6, "axes.labelsize": 6, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "xtick.color": "black", "ytick.color": "black", "legend.fontsize": 6,
    "axes.linewidth": 0.5, "lines.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
})
MM = 1 / 25.4


# ── parse rendered assets ─────────────────────────────────────────────────────
def get_nodes(direction, cls):
    t = (CXC / f"class_{direction}_{cls}_darkgrey.cxc").read_text()
    m = re.search(r"sel #2/R:([\d,]+)@CA", t)
    return sorted(int(x) for x in m.group(1).split(",")) if m else []


def get_edges(direction, cls):
    p = PB / f"{direction}_{cls}.pb"
    if not p.exists():
        return []
    out = []
    for ln in p.read_text().splitlines():
        m = re.findall(r"#2/R:(\d+)@CA", ln)
        if len(m) == 2:
            out.append(tuple(sorted(map(int, m))))
    return out


# ── PCA-aligned coordinates: z = long axis (normal), x/y = membrane plane ──────
ca = pd.read_csv(CA_CSV)
ca["position"] = ca["position"].astype(int)
XYZ_RAW = ca[["x", "y", "z"]].to_numpy(float)
CENTROID = XYZ_RAW.mean(0)                             # 8EFQ Ca centroid (file frame)
XYZ = XYZ_RAW - CENTROID
_, _, Vt = np.linalg.svd(XYZ, full_matrices=False)   # PC1 = largest extent = long axis
P = XYZ @ Vt.T                                         # columns: PC1, PC2, PC3
# pin PC1 sign so extracellular is UP (matches the figure's side-view convention):
# EC markers (N-term, ECL2) must sit at higher PC1 than IC markers (DRY R167, H8).
_idx = {int(p): i for i, p in enumerate(ca["position"])}
_ec = [p for p in (73, 76, 200, 210) if p in _idx]
_ic = [p for p in (167, 262, 344, 348) if p in _idx]
PC1_SIGN = 1.0
if np.mean([P[_idx[p], 0] for p in _ec]) < np.mean([P[_idx[p], 0] for p in _ic]):
    P[:, 0] *= -1
    PC1_SIGN = -1.0
# assign so x,y = membrane plane (PC2,PC3), z = long axis (PC1, extracellular up)
coord = {int(pos): (P[i, 1], P[i, 2], P[i, 0]) for i, pos in enumerate(ca["position"])}
# plot-axis (coord index) -> its unit vector in the 8EFQ file frame
#   0 -> x = PC2 (Vt[1]),  1 -> y = PC3 (Vt[2]),  2 -> z = PC1 (sign-pinned)
AXIS_VEC = {0: Vt[1], 1: Vt[2], 2: PC1_SIGN * Vt[0]}
ALL = np.array([coord[p] for p in ca["position"]])    # for backbone context

# TM centroids (for labels on the top-down view)
sse = pd.read_csv(ROOT / "GPCRdb_OPRM1_table.csv")[["pos", "SSE"]].dropna()
sse["pos"] = sse["pos"].astype(int)
TM_CENTROID = {}
for tm in [f"TM{k}" for k in range(1, 8)]:
    pts = np.array([coord[p] for p in sse.loc[sse.SSE == tm, "pos"] if p in coord])
    if len(pts):
        TM_CENTROID[tm] = pts.mean(0)

# ── Projections ─────────────────────────────────────────────────────────────
# A view is defined by two orthonormal 3-vectors in the coord/PCA frame:
# `rx` (screen right) and `uy` (screen up). Projected coord = (c·rx, c·uy).
#   top-down = PCA membrane plane (x=PC2, y=PC3);
#   side view = derived from a ChimeraX camera (matrices below) so the panel
#   reproduces exactly what that camera shows — line it up with a cartoon snap.
A_LOCAL_TO_PCA = np.array([Vt[1], Vt[2], PC1_SIGN * Vt[0]])   # 8EFQ frame -> PCA


def cam_axes(cam12, mod12):
    """(screen-right, screen-up) unit vectors in the PCA frame for a ChimeraX
    `view matrix camera` + `view matrix models #1` pair (each 12 numbers)."""
    Rc = np.asarray(cam12, float).reshape(3, 4)[:, :3]
    Rm = np.asarray(mod12, float).reshape(3, 4)[:, :3]
    N = Rc.T @ Rm                                    # 8EFQ-local -> screen
    return A_LOCAL_TO_PCA @ N[0], A_LOCAL_TO_PCA @ N[1]


# ChimeraX `view matrix` output defining the SIDE-view camera (edit to re-aim).
SIDE_CAM = [0.17894, -0.015396, 0.98374, 178.08,
            0.22556, 0.97389, -0.025787, -3.0794,
            -0.95765, 0.22651, 0.17774, 34.126]
SIDE_MODEL = [0.96975, -0.18802, 0.15566, -135.48,
              -0.073502, 0.38316, 0.92075, -186.46,
              -0.23276, -0.90434, 0.35775, 95.405]
SIDE_RX, SIDE_UY = cam_axes(SIDE_CAM, SIDE_MODEL)

_I = np.eye(3)
PROJ = [
    dict(rx=_I[0], uy=_I[1], la="x", lb="y",
         sub="top-down (membrane plane)", tm=True),
    dict(rx=SIDE_RX, uy=SIDE_UY, la="x", lb="z",
         sub="side view (ChimeraX camera)", tm=False),
]


def _view_span(proj):
    """Side length (Å) of a view's square data window (matches draw_view)."""
    A2 = ALL @ np.column_stack([proj["rx"], proj["uy"]])
    return max(np.ptp(A2[:, 0]) * 1.12, np.ptp(A2[:, 1]) * 1.12)


SPAN_REF = _view_span(PROJ[0])   # top-down span; node dots scale to this on-page size


def chimerax_matrix_topdown():
    """`view matrix models` (3x4) reproducing the PCA top-down projection."""
    R = np.array([AXIS_VEC[0], AXIS_VEC[1], AXIS_VEC[2]], float)
    if np.linalg.det(R) < 0:
        R[2] = -R[2]
    M = np.column_stack([R, -R @ CENTROID])
    return ",".join(f"{v:.6f}" for v in M.flatten())


def write_chimerax_views(fname="chimerax_view_matrices.cxc"):
    lines = ["# ChimeraX camera orientations matching the density-panel projections.",
             "# Open 8EFQ (receptor = model #1), run ONE block, `view` to fit, snapshot.", "",
             "# --- x-y top-down (generated from the PCA frame) ---",
             f"view matrix models #1,{chimerax_matrix_topdown()}",
             "view", "",
             "# --- x-z side view (the camera you supplied) ---",
             "view matrix camera " + ",".join(f"{v:g}" for v in SIDE_CAM),
             "view matrix models #1," + ",".join(f"{v:g}" for v in SIDE_MODEL), ""]
    (OUT / fname).write_text("\n".join(lines))
    print("saved", fname)

# collect node/edge data once
NODES = {(d, c): get_nodes(d, c) for d in DIRECTIONS for c in CLASSES}
EDGES = {(d, c): get_edges(d, c) for d in DIRECTIONS for c in CLASSES}
NODESET = {k: set(v) for k, v in NODES.items()}
NODE_R = 1.1   # node-dot radius in Angstroms (data units; aspect is equal)


# ════════════════════════════════════════════════════════════════════════════
# 1. Bar plots — nodes / edges per class, LoF vs GoF
# ════════════════════════════════════════════════════════════════════════════
def draw_bars(ax, kind, data_map, ylabel, title=True):
    from matplotlib.patches import Patch
    x = np.arange(len(CLASSES))
    w = 0.38
    for j, d in enumerate(DIRECTIONS):
        vals = [len(data_map[(d, c)]) for c in CLASSES]
        cols = [COLORS[c] for c in CLASSES]
        off = (j - 0.5) * w
        bars = ax.bar(x + off, vals, w, color=cols,
                      alpha=1.0 if d == "lof" else 0.45,
                      edgecolor="black", linewidth=0.4,
                      hatch=None if d == "lof" else "////")
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + max(vals) * 0.02, str(v),
                    ha="center", va="bottom", fontsize=5)
    ax.set_xticks(x)
    ax.set_xticklabels([TITLES[c] for c in CLASSES], rotation=20, ha="right")
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(f"{kind} per efficacy class", fontweight="bold", loc="left",
                     fontsize=6.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    handles = [Patch(facecolor="#888888", edgecolor="black", linewidth=0.4, label="LoF"),
               Patch(facecolor="#888888", edgecolor="black", linewidth=0.4,
                     alpha=0.45, hatch="////", label="GoF")]
    ax.legend(handles=handles, frameon=False, loc="upper left", fontsize=5.5)


def bar_plot(kind, data_map, ylabel, fname):
    fig, ax = plt.subplots(figsize=(75 * MM, 45 * MM))
    draw_bars(ax, kind, data_map, ylabel)
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=600)
    fig.savefig(OUT / fname.replace(".pdf", ".png"), dpi=300)
    plt.close()
    print("saved", fname)


bar_plot("Nodes", NODES, "n nodes", "class_node_counts_bar.pdf")
bar_plot("Edges", EDGES, "n edges", "class_edge_counts_bar.pdf")


def align_blocks(fig, blocks, gap=0.004, thick=0.26):
    """Seat each (top, right) marginal flush against its equal-aspect main box.

    With aspect='equal' the main axes shrinks inside its gridspec cell, so the
    marginals must be repositioned relative to the *drawn* main box: the top
    strip spans the main width just above it, the right strip spans the main
    height just to its right, each `thick` of the main box's size."""
    fig.canvas.draw()
    for ax_m, ax_t, ax_r in blocks:
        p = ax_m.get_position()
        ax_t.set_position([p.x0, p.y1 + gap, p.width, p.height * thick])
        ax_r.set_position([p.x1 + gap, p.y0, p.width * thick, p.height])


# ════════════════════════════════════════════════════════════════════════════
# 2. Density grids — 4 classes (rows) x 3 projections (cols), per direction
# ════════════════════════════════════════════════════════════════════════════
def kde_fill(ax, pts, color, ext):
    """Filled Gaussian KDE if enough non-degenerate points, else scatter only."""
    if len(pts) >= 5:
        try:
            kde = gaussian_kde(pts.T)
            xmin, xmax, ymin, ymax = ext
            xx, yy = np.mgrid[xmin:xmax:120j, ymin:ymax:120j]
            zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
            cmap = mpl.colors.LinearSegmentedColormap.from_list(
                "c", ["white", color])
            ax.contourf(xx, yy, zz, levels=8, cmap=cmap, alpha=0.9, zorder=1)
        except np.linalg.LinAlgError:
            pass
    ax.scatter(pts[:, 0], pts[:, 1], s=4, color=color,
               edgecolors="white", linewidths=0.2, zorder=3)


def density_grid(direction, source, fname, kind):
    fig, axes = plt.subplots(len(CLASSES), len(PROJ),
                             figsize=(44 * len(PROJ) * MM, 150 * MM),
                             squeeze=False)
    for r, cls in enumerate(CLASSES):
        # gather points for this class/direction
        if kind == "node":
            positions = source[(direction, cls)]
            xyz = np.array([coord[p] for p in positions]) if positions else np.empty((0, 3))
        else:  # edge midpoints
            mids = [np.mean([coord[a], coord[b]], axis=0)
                    for a, b in source[(direction, cls)]
                    if a in coord and b in coord]
            xyz = np.array(mids) if mids else np.empty((0, 3))
        for cix, pr in enumerate(PROJ):
            B = np.column_stack([pr["rx"], pr["uy"]])
            ALL2 = ALL @ B
            ext = (ALL2[:, 0].min(), ALL2[:, 0].max(),
                   ALL2[:, 1].min(), ALL2[:, 1].max())
            ax = axes[r, cix]
            ax.scatter(ALL2[:, 0], ALL2[:, 1], s=2.0, color="#b8b8b8",
                       edgecolors="none", zorder=0)
            if len(xyz):
                kde_fill(ax, xyz @ B, COLORS[cls], ext)
            ax.set_aspect("equal")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_linewidth(0.5)
            if r == 0:
                ax.set_title(f"{pr['la']}–{pr['lb']}\n{pr['sub']}", fontsize=5.5)
            if cix == 0:
                n = len(xyz)
                ax.set_ylabel(f"{TITLES[cls]}\n(n={n})", fontsize=6,
                              color=COLORS[cls], fontweight="bold")
    fig.suptitle(f"{kind.capitalize()} density by efficacy class — "
                 f"{DIR_LABEL[direction]}", fontweight="bold", fontsize=7, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(OUT / fname, dpi=600)
    fig.savefig(OUT / fname.replace(".pdf", ".png"), dpi=300)
    plt.close()
    print("saved", fname)


for d in DIRECTIONS:
    density_grid(d, NODES, f"class_node_density_{d}.pdf", "node")
    density_grid(d, EDGES, f"class_edge_density_{d}.pdf", "edge")


# ════════════════════════════════════════════════════════════════════════════
# 3. Overlay — all classes on one set of 3 projections (legible: HPD contours)
#    For each class we draw the smallest region containing FRAC of the density
#    mass as a colored outline + faint fill, so 4 classes coexist legibly.
# ════════════════════════════════════════════════════════════════════════════
GRAD_AMAX = 0.62   # peak opacity of each class's density layer
GRAD_GAMMA = 1.15  # >1 fades the low-density tails faster (cleaner periphery)


def hpd_contour(ax, pts, color, ext):
    """Continuous per-class density as an alpha-blended color gradient.

    Each class paints its hue with per-pixel alpha proportional to its
    (normalized) density, so four classes overlay legibly — overlaps blend
    instead of occluding, and empty regions stay transparent."""
    if len(pts) < 5:
        ax.scatter(pts[:, 0], pts[:, 1], s=6, color=color,
                   edgecolors="white", linewidths=0.3, zorder=4)
        return
    try:
        kde = gaussian_kde(pts.T)
    except np.linalg.LinAlgError:
        return
    xmin, xmax, ymin, ymax = ext
    xx, yy = np.mgrid[xmin:xmax:220j, ymin:ymax:220j]
    zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
    zz = zz / zz.max()
    rgb = mpl.colors.to_rgb(color)
    img = np.empty(zz.shape + (4,))
    img[..., :3] = rgb
    img[..., 3] = np.clip(zz ** GRAD_GAMMA, 0, 1) * GRAD_AMAX
    ax.imshow(img.transpose(1, 0, 2), origin="lower",
              extent=[xmin, xmax, ymin, ymax], interpolation="bilinear",
              aspect="auto", zorder=2)


def draw_node(ax, x, y, classes, r=NODE_R):
    """Full disc for a single-class node; equal pie slices for a multi-class node."""
    from matplotlib.patches import Wedge, Circle
    if len(classes) == 1:
        ax.add_patch(Circle((x, y), r, facecolor=COLORS[classes[0]],
                            edgecolor="white", lw=0.3, zorder=4))
        return
    n = len(classes)
    for k, cls in enumerate(classes):
        ax.add_patch(Wedge((x, y), r, 90 + k * 360 / n, 90 + (k + 1) * 360 / n,
                           facecolor=COLORS[cls], edgecolor="white", lw=0.3,
                           zorder=4))


def _pts_for(direction, source, kind, cls):
    if kind == "node":
        pos = source[(direction, cls)]
        return np.array([coord[p] for p in pos]) if pos else np.empty((0, 3))
    mids = [np.mean([coord[a], coord[b]], axis=0)
            for a, b in source[(direction, cls)] if a in coord and b in coord]
    return np.array(mids) if mids else np.empty((0, 3))


def marginal(ax, vals, color, lo, hi, orient):
    """1D KDE ridge (line only, no shading) for one class along one axis."""
    if len(vals) < 3 or np.ptp(vals) < 1e-6:
        return
    try:
        k = gaussian_kde(vals)
    except np.linalg.LinAlgError:
        return
    xs = np.linspace(lo, hi, 200)
    d = k(xs)
    if orient == "x":
        ax.plot(xs, d, color=color, lw=0.9)
    else:
        ax.plot(d, xs, color=color, lw=0.9)


def draw_view(fig, cell, direction, proj, kind, source, row_title=True):
    """Draw one projection as a jointplot (main + top/right line marginals) into
    a GridSpec `cell`. Returns (ax_m, ax_t, ax_r) for post-hoc alignment."""
    from matplotlib import patheffects
    from matplotlib.ticker import MultipleLocator
    rx, uy = proj["rx"], proj["uy"]
    la, lb, sub = proj["la"], proj["lb"], proj["sub"]
    B = np.column_stack([rx, uy])                 # coord3 @ B -> (screen x, y)
    ALL2 = ALL @ B
    lo_x, hi_x = ALL2[:, 0].min(), ALL2[:, 0].max()
    lo_y, hi_y = ALL2[:, 1].min(), ALL2[:, 1].max()
    mx, my = 0.06 * (hi_x - lo_x), 0.06 * (hi_y - lo_y)
    lo_x, hi_x, lo_y, hi_y = lo_x - mx, hi_x + mx, lo_y - my, hi_y + my
    # square data window: equal Angstrom span on both axes, receptor centered
    span = max(hi_x - lo_x, hi_y - lo_y)
    cx, cy = (lo_x + hi_x) / 2, (lo_y + hi_y) / 2
    lo_x, hi_x, lo_y, hi_y = cx - span / 2, cx + span / 2, cy - span / 2, cy + span / 2
    ext = (ALL2[:, 0].min(), ALL2[:, 0].max(), ALL2[:, 1].min(), ALL2[:, 1].max())
    inner = cell.subgridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4],
                             wspace=0.04, hspace=0.04)
    ax_m = fig.add_subplot(inner[1, 0])
    ax_t = fig.add_subplot(inner[0, 0], sharex=ax_m)
    ax_r = fig.add_subplot(inner[1, 1], sharey=ax_m)

    node_r = NODE_R * span / SPAN_REF        # equalize on-page dot size across views
    ax_m.scatter(ALL2[:, 0], ALL2[:, 1], s=3, color="#cbcbcb",
                 edgecolors="none", zorder=1)
    if kind == "node":
        for p in coord:
            hits = [c for c in CLASSES if p in NODESET[(direction, c)]]
            if hits:
                c3 = np.asarray(coord[p])
                draw_node(ax_m, c3 @ rx, c3 @ uy, hits, r=node_r)
    else:
        for cls in CLASSES:
            xyz = _pts_for(direction, source, kind, cls)
            if len(xyz):
                hpd_contour(ax_m, xyz @ B, COLORS[cls], ext)
    if proj.get("tm"):
        for tm, cen in TM_CENTROID.items():
            ax_m.text(cen @ rx, cen @ uy, tm, fontsize=5, fontweight="bold",
                      ha="center", va="center", color="#222222", zorder=6,
                      path_effects=[patheffects.withStroke(
                          linewidth=1.6, foreground="white")])
    ax_m.set_xlim(lo_x, hi_x); ax_m.set_ylim(lo_y, hi_y)
    ax_m.set_aspect("equal")
    # Angstrom ticks (PCA-frame coords are in A from the receptor centroid):
    # labeled majors every 10 A, unlabeled minors every 5 A.
    for axis in (ax_m.xaxis, ax_m.yaxis):
        axis.set_major_locator(MultipleLocator(10))
        axis.set_minor_locator(MultipleLocator(5))
    ax_m.tick_params(which="major", labelsize=5, length=2.4, width=0.5, pad=1)
    ax_m.tick_params(which="minor", length=1.2, width=0.4)
    ax_m.set_xlabel(f"{la} (Å)", fontsize=5.5)
    ax_m.set_ylabel(f"{lb} (Å)", fontsize=5.5)

    for cls in CLASSES:
        xyz = _pts_for(direction, source, kind, cls)
        if not len(xyz):
            continue
        p2 = xyz @ B
        marginal(ax_t, p2[:, 0], COLORS[cls], lo_x, hi_x, "x")
        marginal(ax_r, p2[:, 1], COLORS[cls], lo_y, hi_y, "y")
    for a in (ax_t, ax_r):
        a.tick_params(which="both", bottom=False, top=False, left=False,
                      right=False, labelbottom=False, labelleft=False)
        for sp in a.spines.values():
            sp.set_visible(False)
    if row_title:
        ax_t.set_title(f"{DIR_LABEL[direction]} · {la}–{lb} ({sub})",
                       fontsize=6, pad=2, loc="left")
    return ax_m, ax_t, ax_r


def joint_overlay(direction, source, fname, kind):
    """Two views (x-y top, x-z side) side by side, each a jointplot."""
    from matplotlib.lines import Line2D
    fig = plt.figure(figsize=(155 * MM, 88 * MM))
    outer = fig.add_gridspec(1, len(PROJ), wspace=0.32)
    blocks = [draw_view(fig, outer[v], direction, proj, kind, source)
              for v, proj in enumerate(PROJ)]
    align_blocks(fig, blocks)
    handles = [Line2D([0], [0], color=COLORS[c], lw=1.6, label=TITLES[c])
               for c in CLASSES]
    fig.legend(handles=handles, frameon=False, ncol=len(CLASSES),
               loc="lower center", fontsize=6, bbox_to_anchor=(0.5, -0.04))
    main_desc = ("nodes colored by class, pie-split if shared"
                 if kind == "node" else "per-class density gradient")
    fig.suptitle(f"{kind.capitalize()} overlay — {DIR_LABEL[direction]} "
                 f"(main: {main_desc}; margins: per-class KDE)",
                 fontweight="bold", fontsize=7, y=1.0)
    fig.savefig(OUT / fname, dpi=600, bbox_inches="tight")
    fig.savefig(OUT / fname.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("saved", fname)


for d in DIRECTIONS:
    joint_overlay(d, NODES, f"class_node_density_overlay_{d}.pdf", "node")
    joint_overlay(d, EDGES, f"class_edge_density_overlay_{d}.pdf", "edge")


# ════════════════════════════════════════════════════════════════════════════
# 4. Composite panel: bars (top) | LoF node overlay (a row per view) |
#    GoF node overlay (a row per view)
# ════════════════════════════════════════════════════════════════════════════
def composite_panel(fname="class_network_composite.pdf"):
    """Fixed 165 x 178 mm sheet. Row 1: node+edge bars (50/50). Rows 2-3: LoF /
    GoF node overlays split into 3 columns with the LEFT column blank for a
    ChimeraX cartoon (matrices in chimerax_view_matrices.cxc); views in cols 2-3.
    Saved WITHOUT tight bbox so the output is exactly the requested size."""
    from matplotlib.lines import Line2D
    fig = plt.figure(figsize=(165 * MM, 158 * MM))
    master = fig.add_gridspec(3, 1, height_ratios=[0.62, 1.0, 1.0], hspace=0.14,
                              left=0.075, right=0.985, top=0.99, bottom=0.072)

    # ── row 1: node + edge count bars, page split 50/50 ──────────────────────
    bcell = master[0].subgridspec(1, 2, wspace=0.32)
    draw_bars(fig.add_subplot(bcell[0, 0]), "Nodes", NODES, "n nodes", title=False)
    draw_bars(fig.add_subplot(bcell[0, 1]), "Edges", EDGES, "n edges", title=False)

    # ── rows 2-3: one VIEW per row (row2 = top-down, row3 = side), LoF & GoF
    #    side by side; col 0 blank for that row's ChimeraX cartoon ────────────
    blocks = []
    for pi, proj in enumerate(PROJ):
        vgs = master[pi + 1].subgridspec(1, 3, wspace=0.28)
        for di, direction in enumerate(("lof", "gof")):
            blocks.append(draw_view(fig, vgs[0, di + 1], direction, proj, "node",
                                    NODES, row_title=False))
    align_blocks(fig, blocks)

    handles = [Line2D([0], [0], color=COLORS[c], lw=1.8, label=TITLES[c])
               for c in CLASSES]
    fig.legend(handles=handles, frameon=False, ncol=len(CLASSES),
               loc="lower center", fontsize=6, bbox_to_anchor=(0.5, 0.005))
    fig.savefig(OUT / fname, dpi=600)
    fig.savefig(OUT / fname.replace(".pdf", ".png"), dpi=300)
    plt.close()
    print("saved", fname)


composite_panel()
write_chimerax_views()

print("\nAll class node/edge plots ->", OUT)

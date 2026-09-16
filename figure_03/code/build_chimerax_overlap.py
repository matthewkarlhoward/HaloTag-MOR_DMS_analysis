#!/usr/bin/env python3
"""
ChimeraX assets for the EC50/Emax overlap, in the same 8EF6 views as figure 3.

  overlap_classes.cxc      3-class spheres: EC50-only (blue), Emax-only (red),
                           both (purple), at the paper's >1 SD LoF threshold.
  overlap_shared_only.cxc  the shared set alone, to show how much of the
                           receptor moves BOTH parameters.
  potency_selective.cxc    continuous map of the potency-specific residual
                           (blue = EC50 shift beyond the Emax-predicted,
                           red = Emax loss with little EC50 shift).
  render_all.cxc           renders every view to PNG.

Cameras are the captured figure-3 views from
curve_refitting/make_chimerax_param_maps.py, so these panels drop straight into
the existing figure without re-orienting anything.
"""
from pathlib import Path
import sys
import numpy as np
import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parent))
from overlap_core import load, COL, ROOT

OUT  = Path(__file__).resolve().parent / "chimerax"
OUT.mkdir(exist_ok=True)
CIFD = ROOT / "structures/raw/experimental"
CIF  = CIFD / "8ef6.cif"
CH   = "R"                      # receptor chain used by figure 3B/C
LIG  = ":MOI"

# ── captured figure-3 cameras (make_chimerax_param_maps.py) ─────────────────
CAM_SIDE = "-0.66813,0.056198,-0.74192,38.605,0.20656,0.97196,-0.11239,102.95,0.7148,-0.22835,-0.661,35.991"
M6_SIDE = np.array([[-0.85867, -0.19853, 0.47252, 221.61],
                    [0.20054, 0.71828, 0.66622, -76.677],
                    [-0.47166, 0.66683, -0.57695, 181.08]])
CAM_TOP = "-0.011909,0.9134,-0.40688,87.959,-0.25662,0.39051,0.88411,208,0.96644,0.11494,0.22975,156.05"
M_DAMGO_TOP = np.array([[-0.86546, -0.31701, -0.38793, 366.93],
                        [-0.39492, -0.04475, 0.91762, 0.74857],
                        [-0.30826, 0.94737, -0.08646, 59.39]])


def ca(cif, chain):
    m = gemmi.read_structure(str(CIFD / cif))[0]
    return {r.seqid.num: np.array([a.pos.x, a.pos.y, a.pos.z])
            for r in m[chain] for a in r if a.name == "CA"}


def align(cif, chain, ref):
    S = ca(cif, chain); common = sorted(set(S) & set(ref))
    P = np.array([S[i] for i in common]); Q = np.array([ref[i] for i in common])
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


REFCA = ca("8ef6.cif", "R")
_Ra, _ta = align("8efq.cif", "R", REFCA)
_Rd, _td = M_DAMGO_TOP[:, :3], M_DAMGO_TOP[:, 3]
_R6 = _Rd @ _Ra.T
M6_TOP = np.hstack([_R6, (_td - _R6 @ _ta)[:, None]])

def flat(M):
    return ",".join(f"{x:.5f}" for i in range(3) for x in M[i])


# ── data ────────────────────────────────────────────────────────────────────
t, meta = load()
sets = {c: sorted(int(p) for p in t.index[t.cls == c]) for c in ("EC50", "Emax", "both")}


def rng(ids):
    """Collapse a sorted id list into ChimeraX range syntax."""
    out, i = [], 0
    while i < len(ids):
        j = i
        while j + 1 < len(ids) and ids[j + 1] == ids[j] + 1:
            j += 1
        out.append(str(ids[i]) if j == i else f"{ids[i]}-{ids[j]}")
        i = j + 1
    return ",".join(out)


HEAD = f"""close session
open {CIF}
delete ~/{CH}
delete :CLR
delete solvent
hide atoms
show /{CH} cartoons
color /{CH} #f2f2f2
show {LIG} atoms
style {LIG} sphere
color {LIG} #f5e663
color {LIG} byhetero
set bgColor white
lighting soft
graphics silhouettes true
"""

VIEWS = {"side": (CAM_SIDE, flat(M6_SIDE)), "top": (CAM_TOP, flat(M6_TOP))}


def view_block(view):
    cam, mat = VIEWS[view]
    return f"view matrix camera {cam}\nview matrix models #1,{mat}\n"


# ── 1. three-class overlap ──────────────────────────────────────────────────
def classes_body(which=("EC50", "Emax", "both")):
    b = ""
    for c in which:
        ids = rng(sets[c])
        b += (f"# {c}  (n={len(sets[c])})\n"
              f"show /{CH}:{ids} atoms\n"
              f"style /{CH}:{ids} sphere\n"
              f"color /{CH}:{ids} {COL[c]}\n")
    return b


for view in VIEWS:
    (OUT / f"overlap_classes_{view}.cxc").write_text(
        f"# EC50/Emax overlap, 3 classes ({view} view)\n"
        f"# blue EC50-only {len(sets['EC50'])} | red Emax-only {len(sets['Emax'])}"
        f" | purple both {len(sets['both'])}\n"
        + HEAD + classes_body() + view_block(view))
    (OUT / f"overlap_shared_only_{view}.cxc").write_text(
        f"# shared (both-parameter) positions only ({view} view), n={len(sets['both'])}\n"
        + HEAD + classes_body(("both",)) + view_block(view))

# ── 2. potency-specific residual ────────────────────────────────────────────
da = OUT / "potency_residual.defattr"
with open(da, "w") as f:
    f.write("attribute: potsel\nrecipient: residues\n")
    for p, v in t.resid.items():
        f.write(f"\t/{CH}:{int(p)}\t{v:.4f}\n")

LIM = round(float(t.resid.abs().quantile(.98)), 2)
for view in VIEWS:
    (OUT / f"potency_selective_{view}.cxc").write_text(
        f"# potency-specific residual ({view} view)\n"
        f"# blue = EC50 shift beyond the Emax-predicted; red = the reverse; "
        f"white = moves on the shared axis\n"
        + HEAD
        + f"open {da}\n"
        + f"color byattribute r:potsel /{CH} palette "
          f"-{LIM},{COL['Emax']}:0,white:{LIM},{COL['EC50']} "
          f"range -{LIM},{LIM} target c noValueColor #cccccc\n"
        + view_block(view))

# ── 3. master render ────────────────────────────────────────────────────────
names = [f"{s}_{v}" for s in ("overlap_classes", "overlap_shared_only",
                              "potency_selective") for v in VIEWS]
lines = ["# Render every overlap view to PNG",
         "windowsize 1600 1200"]          # match the save aspect, or 2D labels clip
for n in names:
    lines += [f"open {OUT}/{n}.cxc",
              "windowsize 1600 1200",
              f"save {OUT}/{n}.png width 1600 height 1200 supersample 3 transparentBackground true"]
(OUT / "render_all.cxc").write_text("\n".join(lines) + "\n")

print(f"positions: EC50-only {len(sets['EC50'])}  Emax-only {len(sets['Emax'])}  "
      f"both {len(sets['both'])}  (of {len(t)})")
print(f"residual palette limit +/-{LIM} log units")
print("wrote:", ", ".join(sorted(p.name for p in OUT.iterdir())))
for c in ("EC50", "Emax", "both"):
    lab = "  ".join(f"{t.wt[p]}{p}" + (f"({t.motif[p]})" if t.motif[p] else "")
                    for p in sets[c])
    print(f"\n{c} (n={len(sets[c])}):\n  {lab}")

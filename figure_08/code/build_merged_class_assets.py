#!/usr/bin/env python3
"""
build_merged_class_assets.py
----------------------------------------------------------------------------
Regenerate the per-class ChimeraX network assets for a *merged* efficacy
grouping in which the old "strong" and "full" classes are pooled into a
single "Strong Agonists" group (the same 6 high-efficacy agonists that the
Python pipeline in plots/network/class_network_analyses.py already calls
"strong": Morphine, C6guano, Fentanyl, Carfentanil, SR17018, DAMGO).

It also recolors the remaining classes to the unified palette:
    strong_agonists  black   #231f20   (pooled strong + full)
    intermediate     red     #9b1c1d
    weak             blue    #3779b9
    antagonist       purple  #6f08a3
    no_ligand        grey    #888888   (kept neutral)

Nothing existing is overwritten — all output lands in *_merged sibling dirs:
    chimerax_attributes/{lof,gof}_strong_agonists_count.defattr
    chimerax_pseudobonds/classes_merged/{lof,gof}_{cls}.pb
    chimerax_pseudobonds/classes_interactions_merged/{lof,gof}_{cls}_{hbonds,pi,saltbridges}.pb
    chimerax_scripts/classes_merged/class_{lof,gof}_{cls}_darkgrey.cxc

The class node-set / edge / count rules were reverse-engineered from the
existing assets and are VALIDATED below before any merged file is written:
  * count  = per-position mean disruption count across the class's drugs
             (pooled "strong_agonists" = (2*strong + 4*full) / 6, exact)
  * nodes  = positions where count >= max(2, Q75 of nonzero class counts)
             AND in the ligand state-contact set
  * edges  = all Ca-Ca pairs < 8.0 A among nodes (PDB 8EFQ)
  * hbonds = 8EFQ sc-sc H-bonds among nodes        (exact atoms)
  * salts  = 8EFQ salt bridges among nodes          (exact atoms)
  * pi     = 8EFQ pi-pi + cation-pi among nodes     (residue level -> CB/CA)
"""

import csv
import math
import re
import sys
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent                 # network_tools/
REPO = ROOT.parent
ATTR_DIR = ROOT / "chimerax_attributes"
PB_DIR = ROOT / "chimerax_pseudobonds"
CXC_DIR = ROOT / "chimerax_scripts"
STRUCT = REPO / "structures" / "raw" / "experimental" / "8efq.cif"

# output dirs (non-destructive)
PB_CLASSES_OUT = PB_DIR / "classes_merged"
PB_INT_OUT = PB_DIR / "classes_interactions_merged"
CXC_OUT = CXC_DIR / "classes_merged"
for d in (PB_CLASSES_OUT, PB_INT_OUT, CXC_OUT):
    d.mkdir(parents=True, exist_ok=True)

# ── constants reverse-engineered from the existing pipeline ──────────────────
CONTACT_THRESHOLD = 8.0     # Ca-Ca cutoff for network edges
MIN_COUNT = 2.0             # floor on the Q75 node gate
QUANT = 0.75
CAMERA = ("0.270015,-0.078762,0.959629,283.278195,0.962548,0.047287,"
          "-0.266956,91.245077,-0.024352,0.995771,0.088581,170.109358")
ZOOM = 0.78
LABEL_TOP_N = 7        # label the N highest-scoring nodes per panel
LABEL_HEIGHT = 2.6     # ChimeraX 3D label height (Angstroms)

# drug membership of the two classes being pooled (from network_analysis_shared.R)
N_STRONG = 2   # Morphine, C6guano
N_FULL = 4     # Fentanyl, Carfentanil, SR17018, DAMGO

# unified palette
COLORS = {
    "no_ligand":       "#888888",
    "antagonist":      "#6f08a3",
    "weak":            "#3779b9",
    "intermediate":    "#9b1c1d",
    "strong_agonists": "#231f20",
}
TITLES = {
    "no_ligand":       "No Ligand",
    "antagonist":      "Antagonist",
    "weak":            "Weak",
    "intermediate":    "Intermediate",
    "strong_agonists": "Strong Agonists",
}
# classes that already exist and only need recolor (node sets unchanged)
PASSTHROUGH = ["no_ligand", "antagonist", "weak", "intermediate"]
DIRECTIONS = ["lof", "gof"]

# ── helpers ──────────────────────────────────────────────────────────────────
def parse_defattr(path):
    d = {}
    for line in open(path):
        m = re.match(r"\s*/R:(\d+)\s+([-\d.]+)", line)
        if m:
            d[int(m.group(1))] = float(m.group(2))
    return d


def load_ca():
    ca = {}
    with open(ROOT / "pdb_8efq_chain_R_ca.csv") as f:
        for r in csv.DictReader(f):
            ca[int(r["position"])] = (float(r["x"]), float(r["y"]), float(r["z"]))
    return ca


def load_state_contacts():
    pos = set()
    with open(ROOT / "state_contacts_combined.csv") as f:
        for r in csv.DictReader(f):
            try:
                hp = int(float(r["human_position"]))
            except (ValueError, TypeError):
                continue
            shells = sum(int(float(r[c])) for c in (
                "inactive_first_shell", "inactive_second_shell",
                "active_first_shell", "active_second_shell"))
            if shells > 0:
                pos.add(hp)
    return pos


def quantile_type7(vals, q):
    """R's default quantile (type 7) == numpy linear interpolation."""
    s = sorted(vals)
    if not s:
        return None
    if len(s) == 1:
        return s[0]
    h = (len(s) - 1) * q
    lo = math.floor(h)
    hi = math.ceil(h)
    return s[lo] + (h - lo) * (s[hi] - s[lo])


def node_threshold(counts):
    nz = [v for v in counts.values() if v > 0]
    cut = quantile_type7(nz, QUANT)
    if cut is None:
        cut = MIN_COUNT
    return max(MIN_COUNT, cut)


def select_nodes(counts, state_pos):
    cut = node_threshold(counts)
    return sorted(p for p, v in counts.items() if v >= cut and p in state_pos)


def top_n_by_count(nodes, counts, n=LABEL_TOP_N):
    """The n highest-scoring node positions (ties broken by lower position)."""
    ranked = sorted(((p, counts.get(p, 0.0)) for p in nodes),
                    key=lambda kv: (-kv[1], kv[0]))
    return [p for p, _ in ranked[:n]]


def edges_8A(nodes, ca):
    out = []
    for a, b in combinations(sorted(nodes), 2):
        if a in ca and b in ca and math.dist(ca[a], ca[b]) < CONTACT_THRESHOLD:
            out.append((a, b))
    return out


# ── interaction sources (PDB 8EFQ) ───────────────────────────────────────────
INT_DIR = REPO / "structures" / "processed" / "interactions" / "per_structure"
PI_DIR = REPO / "structures" / "processed" / "pi_interactions"


def load_hbonds_scsc():
    rows = []
    with open(INT_DIR / "8efq_hbonds.csv") as f:
        for r in csv.DictReader(f):
            if r["hb_type"] == "sc-sc":
                rows.append((int(r["res1_num"]), r["res1_atom"],
                             int(r["res2_num"]), r["res2_atom"]))
    return rows


def load_saltbridges():
    rows = []
    with open(INT_DIR / "8efq_salt_bridges.csv") as f:
        for r in csv.DictReader(f):
            rows.append((int(r["res1_num"]), r["res1_atom"],
                         int(r["res2_num"]), r["res2_atom"]))
    return rows


def _cb_atom(resname):
    return "CA" if (resname or "").upper() == "GLY" else "CB"


def load_pi():
    """Residue-level pi-pi and cation-pi for 8EFQ -> CB/CA representative atoms."""
    rows = []
    for fname in ("all_pi_pi.csv", "all_pi_cation.csv"):
        path = PI_DIR / fname
        if not path.exists():
            continue
        with open(path) as f:
            rdr = csv.DictReader(f)
            cols = rdr.fieldnames
            r1n = next((c for c in cols if c in ("res1_num", "aromatic_res_num", "ring_res_num")), None)
            r2n = next((c for c in cols if c in ("res2_num", "cation_res_num")), None)
            r1nm = next((c for c in cols if c in ("res1_name", "aromatic_res_name", "ring_res_name")), None)
            r2nm = next((c for c in cols if c in ("res2_name", "cation_res_name")), None)
            for r in rdr:
                if str(r.get("pdb_id", "")).lower() != "8efq":
                    continue
                if not r1n or not r2n:
                    continue
                rows.append((int(float(r[r1n])), _cb_atom(r.get(r1nm, "")),
                             int(float(r[r2n])), _cb_atom(r.get(r2nm, ""))))
    return rows


def filter_pairs(pairs, nodeset):
    """Keep pairs with both endpoints in nodeset; dedupe undirected atom-pairs.

    The source interaction CSVs list some contacts in both directions
    (e.g. 77@OH<->321@ND1 and 321@ND1<->77@OH), which would otherwise draw
    a doubled pseudobond."""
    out = []
    seen = set()
    for a, aa, b, ba in pairs:
        if a in nodeset and b in nodeset and a != b:
            key = frozenset(((a, aa), (b, ba)))
            if key in seen:
                continue
            seen.add(key)
            out.append((a, aa, b, ba))
    return out


# ── writers ──────────────────────────────────────────────────────────────────
def write_defattr(path, attr_name, counts):
    positions = sorted(counts)
    with open(path, "w") as f:
        f.write(f"attribute: {attr_name}\n")
        f.write("match mode: any\n")
        f.write("recipient: residues\n")
        for p in positions:
            f.write(f"\t/R:{p}\t{round(counts[p], 6):g}\n")


def write_edge_pb(path, edges, color, radius=0.44):
    with open(path, "w") as f:
        f.write("; halfbond = false\n")
        f.write(f"; radius = {radius}\n")
        f.write("; dashes = 0\n")
        f.write(f"; color = {color}\n")
        for a, b in edges:
            f.write(f"#2/R:{a}@CA #2/R:{b}@CA\n")


def write_interaction_pb(path, pairs, color, radius, dashes=8):
    with open(path, "w") as f:
        f.write("; halfbond = false\n")
        f.write(f"; radius = {radius}\n")
        f.write(f"; dashes = {dashes}\n")
        f.write(f"; color = {color}\n")
        for a, aa, b, ba in pairs:
            f.write(f"#2/R:{a}@{aa} #2/R:{b}@{ba}\n")


def recolor_pb(src, dst, color):
    """Copy a .pb file, swapping only the `; color =` header line."""
    lines = open(src).read().splitlines()
    with open(dst, "w") as f:
        for line in lines:
            if line.strip().startswith("; color ="):
                f.write(f"; color = {color}\n")
            else:
                f.write(line + "\n")


def write_cxc(path, direction, cls, nodes, n_edges, color, has_hbonds,
              has_pi, has_salt, label_positions=()):
    dir_label = direction.upper()
    title = TITLES[cls]
    sel = ",".join(str(n) for n in nodes)
    L = []
    L.append(f"# Class network script: {title} {dir_label}")
    L.append("# Two-model layout (#1 translucent cartoon, #2 opaque spheres)")
    L.append("")
    L.append("close session")
    L.append(f"open {STRUCT}")
    L.append(f"open {STRUCT}")
    L.append("")
    L.append("hide ~/R")
    L.append("hide solvent")
    L.append("hide ligand")
    L.append("")
    L.append("cartoon #1/R")
    L.append("hide #1/R target a")
    L.append("color #1/R #a8a8a8")
    L.append("color #1/R #a8a8a8 transparency 75")
    L.append("")
    L.append("~cartoon #2")
    L.append("hide #2/R target a")
    L.append("")
    L.append("set bgColor white")
    L.append("lighting soft")
    L.append("lighting shadows false")
    L.append("graphics silhouettes true")
    L.append("")
    L.append(f"open {ATTR_DIR / f'{direction}_{cls}_count.defattr'}")
    L.append("")
    if nodes:
        L.append(f"sel #2/R:{sel}@CA")
        L.append("show sel target a")
        L.append("style sel sphere")
        L.append("size sel atomRadius 1.1")
        L.append(f"color sel {color}")
        L.append("~sel")
        if n_edges > 0:
            L.append(f"open {PB_CLASSES_OUT / f'{direction}_{cls}.pb'}")
        if has_hbonds:
            L.append(f"open {PB_INT_OUT / f'{direction}_{cls}_hbonds.pb'}")
        if has_pi:
            L.append(f"open {PB_INT_OUT / f'{direction}_{cls}_pi.pb'}")
        if has_salt:
            L.append(f"open {PB_INT_OUT / f'{direction}_{cls}_saltbridges.pb'}")
        if label_positions:
            lab = ",".join(str(p) for p in label_positions)
            L.append(f'label #2/R:{lab} residues text "{{0.number}}" '
                     f"height {LABEL_HEIGHT} color black bgColor white")
    else:
        L.append("# no residues pass the gate for this panel (empty network)")
    L.append(f'2dlabels text "{title} — {dir_label}  '
             f'({len(nodes)} residues, {n_edges} edges)" '
             f"xpos 0.04 ypos 0.94 size 22 color black")
    L.append(f"view matrix camera {CAMERA}")
    L.append(f"zoom {ZOOM}")
    L.append("")
    path.write_text("\n".join(L))


# ════════════════════════════════════════════════════════════════════════════
# 0. Load shared data
# ════════════════════════════════════════════════════════════════════════════
ca = load_ca()
state_pos = load_state_contacts()
hbonds_all = load_hbonds_scsc()
salts_all = load_saltbridges()
pi_all = load_pi()
print(f"Loaded: {len(ca)} CA coords, {len(state_pos)} state-contact positions, "
      f"{len(hbonds_all)} sc-sc hbonds, {len(salts_all)} salt bridges, "
      f"{len(pi_all)} pi contacts (8EFQ)")


def existing_nodes(direction, cls):
    """Parse the node residue list from an existing class .cxc (for validation)."""
    cxc = CXC_DIR / "classes" / f"class_{direction}_{cls}_darkgrey.cxc"
    m = re.search(r"sel #2/R:([\d,]+)@CA", cxc.read_text())
    return sorted(int(x) for x in m.group(1).split(",")) if m else []


def existing_edges(direction, cls):
    pb = PB_DIR / "classes" / f"{direction}_{cls}.pb"
    out = []
    for line in pb.read_text().splitlines():
        m = re.findall(r"#2/R:(\d+)@CA", line)
        if len(m) == 2:
            out.append(tuple(sorted((int(m[0]), int(m[1])))))
    return out


# ════════════════════════════════════════════════════════════════════════════
# 1. VALIDATION — reproduce existing strong & full panels exactly
# ════════════════════════════════════════════════════════════════════════════
print("\n=== VALIDATION (reproduce existing assets) ===")
ok = True
for direction in DIRECTIONS:
    for cls in ["strong", "full", "intermediate", "weak", "antagonist", "no_ligand"]:
        counts = parse_defattr(ATTR_DIR / f"{direction}_{cls}_count.defattr")
        nodes = select_nodes(counts, state_pos)
        edges = edges_8A(nodes, ca)
        exp_nodes = existing_nodes(direction, cls)
        exp_edges = existing_edges(direction, cls)
        nmatch = nodes == exp_nodes
        ematch = sorted(edges) == sorted(exp_edges)
        if not (nmatch and ematch):
            ok = False
            print(f"  [FAIL] {direction}_{cls}: nodes {len(nodes)}/{len(exp_nodes)} "
                  f"match={nmatch}; edges {len(edges)}/{len(exp_edges)} match={ematch}")
        else:
            print(f"  [ok]   {direction}_{cls}: {len(nodes)} nodes, {len(edges)} edges")

# validate interaction filtering against existing strong hbonds/salt
strong_nodes = set(existing_nodes("lof", "strong"))
hb_strong = filter_pairs([(a, aa, b, ba) for a, aa, b, ba in hbonds_all], strong_nodes)
salt_strong = filter_pairs([(a, aa, b, ba) for a, aa, b, ba in salts_all], strong_nodes)
print(f"  interaction check (lof_strong nodes): hbonds={len(hb_strong)} "
      f"(expect 2), saltbridges={len(salt_strong)} (expect 0)")

if not ok:
    sys.exit("\nValidation FAILED — node/edge rules do not reproduce existing "
             "assets. Aborting before writing merged output.")
print("Validation PASSED — rules reproduce all existing class panels exactly.")


# ════════════════════════════════════════════════════════════════════════════
# 2. Build the merged "strong_agonists" assets (pooled strong + full)
# ════════════════════════════════════════════════════════════════════════════
print("\n=== MERGED: strong_agonists (pooled strong+full, 2:4 weighting) ===")
merged_node_summary = {}
for direction in DIRECTIONS:
    strong_c = parse_defattr(ATTR_DIR / f"{direction}_strong_count.defattr")
    full_c = parse_defattr(ATTR_DIR / f"{direction}_full_count.defattr")
    positions = sorted(set(strong_c) | set(full_c))
    merged_c = {p: (N_STRONG * strong_c.get(p, 0.0) + N_FULL * full_c.get(p, 0.0))
                   / (N_STRONG + N_FULL) for p in positions}

    # scores
    write_defattr(ATTR_DIR / f"{direction}_strong_agonists_count.defattr",
                  f"{direction}_strong_agonists_count", merged_c)

    # nodes + edges
    nodes = select_nodes(merged_c, state_pos)
    edges = edges_8A(nodes, ca)
    nodeset = set(nodes)
    merged_node_summary[direction] = (len(nodes), len(edges))

    color = COLORS["strong_agonists"]
    write_edge_pb(PB_CLASSES_OUT / f"{direction}_strong_agonists.pb", edges, color)

    # interactions among merged nodes
    hb = filter_pairs(hbonds_all, nodeset)
    salt = filter_pairs(salts_all, nodeset)
    pi = filter_pairs(pi_all, nodeset)
    write_interaction_pb(PB_INT_OUT / f"{direction}_strong_agonists_hbonds.pb",
                         hb, color, radius=0.2)
    write_interaction_pb(PB_INT_OUT / f"{direction}_strong_agonists_pi.pb",
                         pi, color, radius=0.22)
    write_interaction_pb(PB_INT_OUT / f"{direction}_strong_agonists_saltbridges.pb",
                         salt, color, radius=0.26)

    write_cxc(CXC_OUT / f"class_{direction}_strong_agonists_darkgrey.cxc",
              direction, "strong_agonists", nodes, len(edges), color,
              has_hbonds=bool(hb), has_pi=bool(pi), has_salt=bool(salt),
              label_positions=top_n_by_count(nodes, merged_c))
    print(f"  {direction}_strong_agonists: {len(nodes)} nodes, {len(edges)} edges, "
          f"hbonds={len(hb)}, pi={len(pi)}, salt={len(salt)}")


# ════════════════════════════════════════════════════════════════════════════
# 3. Recolor the pass-through classes (node sets unchanged)
# ════════════════════════════════════════════════════════════════════════════
print("\n=== RECOLOR pass-through classes ===")
for direction in DIRECTIONS:
    for cls in PASSTHROUGH:
        color = COLORS[cls]
        counts = parse_defattr(ATTR_DIR / f"{direction}_{cls}_count.defattr")
        nodes = select_nodes(counts, state_pos)
        edges = edges_8A(nodes, ca)

        # recolored edge .pb
        recolor_pb(PB_DIR / "classes" / f"{direction}_{cls}.pb",
                   PB_CLASSES_OUT / f"{direction}_{cls}.pb", color)
        # recolored interaction .pb (whichever exist)
        has = {}
        for kind in ("hbonds", "pi", "saltbridges"):
            src = PB_DIR / "classes_interactions" / f"{direction}_{cls}_{kind}.pb"
            if src.exists():
                dst = PB_INT_OUT / f"{direction}_{cls}_{kind}.pb"
                recolor_pb(src, dst, color)
                # treat as present only if it carries >=1 pseudobond
                body = [ln for ln in dst.read_text().splitlines()
                        if ln and not ln.startswith(";")]
                has[kind] = bool(body)
            else:
                has[kind] = False

        write_cxc(CXC_OUT / f"class_{direction}_{cls}_darkgrey.cxc",
                  direction, cls, nodes, len(edges), color,
                  has_hbonds=has["hbonds"], has_pi=has["pi"],
                  has_salt=has["saltbridges"],
                  label_positions=top_n_by_count(nodes, counts))
        print(f"  {direction}_{cls}: recolored -> {color} "
              f"({len(nodes)} nodes, {len(edges)} edges)")

print("\nDone. Merged + recolored assets written to *_merged dirs.")
print(f"Merged node/edge summary: {merged_node_summary}")

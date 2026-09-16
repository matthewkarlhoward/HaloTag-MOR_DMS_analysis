#!/usr/bin/env python3
"""
Panel D: EC50/Emax overlap as CARTOON coloring on 8EF6, at several LoF cutoffs.
  purple = both, blue = EC50-only, red = Emax-only, grey = neither.
Cutoff = K * (1 SD of the synonymous distribution) for each parameter.
Reuses overlap_core.load() (dec50, eloss, sd_ec, sd_em); figure-3 side view.
Outputs one chimerax/overlap_classes_cartoon_<K>SD.cxc per cutoff.
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from overlap_core import load, COL, ROOT

OUT = Path(__file__).resolve().parent / "chimerax"
CIF = ROOT / "structures/raw/experimental/8ef6.cif"
CAM_SIDE = ("-0.66813,0.056198,-0.74192,38.605,0.20656,0.97196,-0.11239,102.95,"
            "0.7148,-0.22835,-0.661,35.991")
M6_SIDE  = ("-0.85867,-0.19853,0.47252,221.61000,0.20054,0.71828,0.66622,-76.67700,"
            "-0.47166,0.66683,-0.57695,181.08000")
KS = [1.0, 1.5, 2.0]

t, meta = load()

def spec(ids):
    return "/R:" + ",".join(str(i) for i in ids)

for K in KS:
    lof_ec = t.dec50 > K * meta["sd_ec"]
    lof_em = t.eloss > K * meta["sd_em"]
    cls = np.where(lof_ec & lof_em, "both",
          np.where(lof_ec, "EC50", np.where(lof_em, "Emax", "neither")))
    sets = {c: sorted(int(p) for p, cc in zip(t.index, cls) if cc == c)
            for c in ("EC50", "Emax", "both")}
    tag = f"{K:g}".replace(".", "p")
    L = [
        "close session", f"open {CIF}",
        "delete ~/R", "delete :CLR", "delete solvent",
        "hide atoms", "show /R cartoons", "color /R #e8e8e8",
        "show :MOI atoms", "style :MOI sphere", "color :MOI #f5e663", "color :MOI byhetero",
        f"color {spec(sets['both'])} {COL['both']} target c",
        f"color {spec(sets['EC50'])} {COL['EC50']} target c",
        f"color {spec(sets['Emax'])} {COL['Emax']} target c",
        "set bgColor white", "lighting soft", "graphics silhouettes true",
        f"view matrix camera {CAM_SIDE}",
        f"view matrix models #1,{M6_SIDE}",
    ]
    (OUT / f"overlap_classes_cartoon_{tag}SD.cxc").write_text("\n".join(L) + "\n")
    print(f"{K:>4g} SD  ->  both {len(sets['both']):3d}  |  EC50-only {len(sets['EC50']):3d}  "
          f"|  Emax-only {len(sets['Emax']):3d}  ->  overlap_classes_cartoon_{tag}SD.cxc")

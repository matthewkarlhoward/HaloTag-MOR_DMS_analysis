#!/usr/bin/env python3
"""
Position-axis view of morphine vs fentanyl ΔEmax bias.

For each receptor position outside the contact-dominated zone
(Cα > 7 Å from both Morphine and Fentanyl ligand heavy atoms):

  bias = ΔEmax(morphine) − ΔEmax(fentanyl)

Positive → morphine more damaged (= fentanyl-preserving)
Negative → fentanyl more damaged (= morphine-preserving)

Plotted as a lollipop at each position; colored red/blue by bias direction.
In-pocket positions (≤ 7 Å) shown as small grey points along y = 0 for context.

Output: plots/scatter/ec50_emax/morphine_fentanyl_emax_bias_by_position.{pdf,png}
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

ROOT    = Path(__file__).resolve().parents[3]
OUT     = ROOT / "plots" / "scatter" / "ec50_emax"
CIF_DIR = ROOT / "structures" / "raw" / "experimental"

plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 6,
    "axes.labelsize": 6, "axes.titlesize": 6,
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "text.color": "black",
    "axes.edgecolor": "black", "axes.labelcolor": "black",
    "axes.titlecolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "lines.linewidth": 0.5, "patch.linewidth": 0.5,
})

MOR_COL  = "#9b1c1d"
FENT_COL = "#225a9c"
BG_COL   = "#bdbdbd"

MIN_VARIANTS = 5
DIST_THRESH  = 7.0   # Å — Cα to nearest ligand heavy atom

STRUCT = {
    "Morphine": dict(pdb="8ef6", ligand_code="MOI", receptor_chain="R"),
    "Fentanyl": dict(pdb="8ef5", ligand_code="7V7", receptor_chain="R"),
}

def calpha_to_smallmol(pdb_id, ligand_code, receptor_chain):
    parser = MMCIFParser(QUIET=True)
    s = parser.get_structure(pdb_id, CIF_DIR / f"{pdb_id}.cif")
    lig = [a for model in s for chain in model for res in chain
           if res.get_resname() == ligand_code
           for a in res.get_atoms() if a.element != "H"]
    rows = []
    for model in s:
        for chain in model:
            if chain.id != receptor_chain: continue
            for res in chain:
                if not is_aa(res, standard=True) or "CA" not in res: continue
                ca = res["CA"]
                d = min((ca - la) for la in lig)
                rows.append(dict(position=res.id[1], ca_dist=float(d)))
    return pd.DataFrame(rows)

dm = calpha_to_smallmol(STRUCT["Morphine"]["pdb"],
                        STRUCT["Morphine"]["ligand_code"],
                        STRUCT["Morphine"]["receptor_chain"]
                        ).rename(columns={"ca_dist":"ca_mor"})
dfn = calpha_to_smallmol(STRUCT["Fentanyl"]["pdb"],
                         STRUCT["Fentanyl"]["ligand_code"],
                         STRUCT["Fentanyl"]["receptor_chain"]
                         ).rename(columns={"ca_dist":"ca_fen"})
dist = dm.merge(dfn, on="position", how="inner")
dist["ca_min"] = dist[["ca_mor","ca_fen"]].min(axis=1)

# ── Per-position ΔEmax ───────────────────────────────────────────────────────
dms = pd.read_csv(ROOT / "dms_scores/composite_dms_scores.csv")
syn = dms[dms.type == "synonymous"]
mis = dms[dms.type == "missense"]

def per_position_demax(drug):
    emax_col = f"{drug}_emax"
    ct_col   = f"{drug}_curve_type"
    syn_emax = syn.loc[syn[ct_col]=="sigmoid", emax_col].mean()
    mis_sig  = mis[mis[ct_col]=="sigmoid"]
    g = mis_sig.groupby("position").agg(
        emax = (emax_col, "mean"),
        n    = (emax_col, "count"),
        wt   = ("wildtype", "first"),
    ).dropna()
    g = g[g.n >= MIN_VARIANTS]
    g["delta"] = g.emax - syn_emax
    return g

mor = per_position_demax("Morphine").rename(columns={"delta":"dMor"})
fen = per_position_demax("Fentanyl").rename(columns={"delta":"dFen"})
df = mor[["dMor","wt"]].join(fen[["dFen"]], how="inner").reset_index()
df = df.merge(dist[["position","ca_min"]], on="position", how="inner")
df["bias"] = df.dMor - df.dFen
df["outside"] = df.ca_min > DIST_THRESH

# ── Figure ───────────────────────────────────────────────────────────────────
MM = 1/25.4
fig, ax = plt.subplots(figsize=(165*MM, 55*MM))

# Inside-pocket positions (grey, on baseline)
inside = df[~df.outside]
ax.scatter(inside.position, np.zeros(len(inside)),
           s=4, color=BG_COL, edgecolors="none", alpha=0.7, zorder=1)

# Outside-pocket lollipops, colored by bias direction
out = df[df.outside].copy()
out_pos_mor  = out[out.bias > 0]
out_pos_fent = out[out.bias < 0]

# Stems
ax.vlines(out_pos_mor.position,  ymin=0, ymax=out_pos_mor.bias,
          color=MOR_COL,  lw=0.5, alpha=0.85, zorder=2)
ax.vlines(out_pos_fent.position, ymin=0, ymax=out_pos_fent.bias,
          color=FENT_COL, lw=0.5, alpha=0.85, zorder=2)
# Heads
ax.scatter(out_pos_mor.position,  out_pos_mor.bias,
           s=6, color=MOR_COL,  edgecolors="none", alpha=0.95, zorder=3)
ax.scatter(out_pos_fent.position, out_pos_fent.bias,
           s=6, color=FENT_COL, edgecolors="none", alpha=0.95, zorder=3)

# Label the top biasing residues on each side
TOP = 6
top_mor  = out_pos_mor.nlargest(TOP, "bias")
top_fent = out_pos_fent.nsmallest(TOP, "bias")
for _, r in top_mor.iterrows():
    ax.annotate(f"{r.wt}{int(r.position)}",
                (r.position, r.bias),
                xytext=(0, 4), textcoords="offset points",
                fontsize=5, color=MOR_COL, ha="center", va="bottom")
for _, r in top_fent.iterrows():
    ax.annotate(f"{r.wt}{int(r.position)}",
                (r.position, r.bias),
                xytext=(0, -4), textcoords="offset points",
                fontsize=5, color=FENT_COL, ha="center", va="top")

ax.axhline(0, color="black", lw=0.5, zorder=4)
ax.set_xlabel("Residue position")
ax.set_ylabel(r"Bias = $\Delta E_\mathrm{max,Mor} - \Delta E_\mathrm{max,Fent}$")
ax.set_title(rf"Per-position morphine–fentanyl Emax bias  "
             rf"(C$\alpha$ > {DIST_THRESH:.0f} Å, both ligands)",
             color="black")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
for s in ("left","bottom"):
    ax.spines[s].set_linewidth(0.5)

# Legend
handles = [
    plt.Line2D([0],[0], marker="o", color=MOR_COL, markerfacecolor=MOR_COL,
               markeredgecolor="none", markersize=3.5, lw=0.5,
               label="Morphine-biased (Mor LoF stronger)"),
    plt.Line2D([0],[0], marker="o", color=FENT_COL, markerfacecolor=FENT_COL,
               markeredgecolor="none", markersize=3.5, lw=0.5,
               label="Fentanyl-biased (Fent LoF stronger)"),
    plt.Line2D([0],[0], marker="o", color="none", markerfacecolor=BG_COL,
               markeredgecolor="none", markersize=2,
               label=rf"in or near pocket (C$\alpha \leq$ {DIST_THRESH:.0f} Å)"),
]
ax.legend(handles=handles, frameon=False, loc="upper right", fontsize=5)

fig.tight_layout()
out_pdf = OUT / "morphine_fentanyl_emax_bias_by_position.pdf"
out_png = OUT / "morphine_fentanyl_emax_bias_by_position.png"
fig.savefig(out_pdf, dpi=600)
fig.savefig(out_png, dpi=600)
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")
plt.close()

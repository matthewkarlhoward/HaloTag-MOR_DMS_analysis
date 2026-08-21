#!/usr/bin/env python3
"""
Overlay of all three ligand decay curves on a single panel.

For each of Morphine, Fentanyl, DAMGO:
  - per-position mean Δlog[EC50] (sigmoid missense, vs syn mean)
  - Cα distance from receptor residue to nearest ligand heavy atom
  - same filters as distance_vs_delta_ec50.R
        (Fentanyl: mean_ec50_logM ≥ -8.6; Morphine: delta_ec50 ≥ -1)
  - exponential decay fit y = a · exp(-b · x) + c on raw per-position points

Output: plots/scatter/distance_shell/distance_vs_delta_ec50_calpha_overlay.{pdf,png}
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pathlib import Path
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

ROOT = Path(__file__).resolve().parents[3]
OUT  = ROOT / "plots" / "scatter" / "distance_shell"
CIF_DIR = ROOT / "structures" / "raw" / "experimental"

# ── Style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 6,
    "mathtext.default": "regular",
    "axes.labelsize": 6, "axes.titlesize": 6,
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "text.color": "black",
    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "axes.titlecolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "lines.linewidth": 0.5,
    "patch.linewidth": 0.5,
})

# Per-ligand colors — readable when overlaid
LIG_COL = {
    "Morphine": "#9b1c1d",   # deep red
    "Fentanyl": "#225a9c",   # blue
    "DAMGO":    "#2c7c2c",   # green
}

SHIFT  = {"Fentanyl":-13, "Morphine":-12.5, "DAMGO":-10}
STRUCT = {"Morphine":dict(pdb="8ef6", ligand_code="MOI", receptor_chain="R"),
          "Fentanyl":dict(pdb="8ef5", ligand_code="7V7", receptor_chain="R"),
          "DAMGO":   dict(pdb="8efq", ligand_chain="P",  receptor_chain="R")}

MIN_VARIANTS = 5

# ── Cα distance computation ──────────────────────────────────────────────────
def ligand_heavy_atoms_smallmol(structure, ligand_code):
    atoms = []
    for model in structure:
        for chain in model:
            for res in chain:
                hetflag, _, _ = res.id
                if hetflag.strip() == "":
                    continue
                if res.get_resname() == ligand_code:
                    for a in res.get_atoms():
                        if a.element != "H":
                            atoms.append(a)
    return atoms

def ligand_heavy_atoms_peptide(structure, ligand_chain):
    atoms = []
    for model in structure:
        for chain in model:
            if chain.id != ligand_chain:
                continue
            for res in chain:
                for a in res.get_atoms():
                    if a.element != "H":
                        atoms.append(a)
    return atoms

def calpha_distances(drug):
    s = STRUCT[drug]
    pdb_id   = s["pdb"]
    chain_id = s["receptor_chain"]
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id, CIF_DIR / f"{pdb_id}.cif")
    if "ligand_code" in s:
        lig_atoms = ligand_heavy_atoms_smallmol(structure, s["ligand_code"])
    else:
        lig_atoms = ligand_heavy_atoms_peptide(structure, s["ligand_chain"])
    rows = []
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for res in chain:
                if not is_aa(res, standard=True): continue
                if "CA" not in res: continue
                ca = res["CA"]
                d = min((ca - la) for la in lig_atoms)
                rows.append(dict(position=res.id[1], ca_dist=float(d)))
    return pd.DataFrame(rows)

# ── Per-ligand DRC aggregation (matches R-script filters) ────────────────────
dms = pd.read_csv(ROOT/"dms_scores/composite_dms_scores.csv")
syn = dms[dms.type == "synonymous"]
mis = dms[dms.type == "missense"]

def per_drug_df(drug):
    ec50_col = f"{drug}_ec50"
    ct_col   = f"{drug}_curve_type"
    syn_sig = syn[syn[ct_col]=="sigmoid"].copy()
    syn_sig["logM"] = syn_sig[ec50_col] + SHIFT[drug]
    syn_logM = syn_sig.logM.mean()
    mis_sig = mis[mis[ct_col]=="sigmoid"].copy()
    mis_sig["logM"] = mis_sig[ec50_col] + SHIFT[drug]
    pos = mis_sig.groupby("position").agg(
        mean_ec50=("logM","mean"),
        n        =("logM","count"),
    ).dropna()
    pos = pos[pos.n >= MIN_VARIANTS]
    pos["delta_ec50"] = pos.mean_ec50 - syn_logM
    pos = pos.reset_index()
    ca = calpha_distances(drug)
    df = ca.merge(pos, on="position", how="inner")
    if drug == "Fentanyl":
        df = df[df.mean_ec50 >= -8.6]
    if drug == "Morphine":
        df = df[df.delta_ec50 >= -1]
    df = df[df.delta_ec50.notna() & np.isfinite(df.delta_ec50)]
    return df

# ── Exponential decay fit ────────────────────────────────────────────────────
def expdecay(x, a, b, c): return a*np.exp(-b*x) + c

def fit_decay(df):
    x = df.ca_dist.values
    y = df.delta_ec50.values
    popt, _ = curve_fit(expdecay, x, y, p0=[0.5, 0.05, 0.0], maxfev=20000)
    a, b, c = popt
    ypred = expdecay(x, *popt)
    ss_res = np.sum((y-ypred)**2); ss_tot = np.sum((y-y.mean())**2)
    r2 = 1 - ss_res/ss_tot
    return dict(a=a, b=b, c=c, half_decay=np.log(2)/b, r2=r2, n=len(df))

# ── Figure ───────────────────────────────────────────────────────────────────
MM = 1/25.4
fig, ax = plt.subplots(figsize=(60*MM, 40*MM))

fits = {}
for drug in ["Morphine","Fentanyl","DAMGO"]:
    df  = per_drug_df(drug)
    fit = fit_decay(df)
    fits[drug] = fit
    col = LIG_COL[drug]

    # Faint per-position points
    ax.scatter(df.ca_dist, df.delta_ec50,
               s=2, color=col, edgecolors="none",
               alpha=0.25, zorder=2)

    # Fit curve
    xs = np.linspace(df.ca_dist.min(), df.ca_dist.max(), 300)
    ys = expdecay(xs, fit["a"], fit["b"], fit["c"])
    ax.plot(xs, ys, color=col, lw=0.5, zorder=5, label=drug)

ax.axhline(0, color="black", lw=0.5, zorder=1)
ax.set_xlabel(r"C$\alpha$ distance to ligand (Å)")
ax.set_ylabel(r"$\Delta$EC50 (log[M])")
ax.legend(frameon=False, loc="upper right",
          handlelength=1.0, handletextpad=0.4, borderaxespad=0.3)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
for s in ("left","bottom"):
    ax.spines[s].set_linewidth(0.5)

fig.subplots_adjust(left=0.19, bottom=0.195, right=0.985, top=0.98)
out_pdf = OUT / "distance_vs_delta_ec50_calpha_overlay.pdf"
out_png = OUT / "distance_vs_delta_ec50_calpha_overlay.png"
fig.savefig(out_pdf, dpi=600)
fig.savefig(out_png, dpi=600)
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")
print()
print("Fit parameters (Cα distance):")
for drug, fit in fits.items():
    print(f"  {drug:<10}  a={fit['a']:.3f}  b={fit['b']:.3f} Å⁻¹  "
          f"c={fit['c']:.3f}  d½={fit['half_decay']:.2f} Å  "
          f"R²={fit['r2']:.3f}  n={fit['n']}")
crossover = np.mean([np.log(f["a"]/f["c"])/f["b"] for f in fits.values()])
print(f"\nMean crossover (contact = plateau): {crossover:.2f} Å")
plt.close()

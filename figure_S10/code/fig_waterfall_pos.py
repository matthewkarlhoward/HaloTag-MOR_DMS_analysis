"""Position-level waterfall. Compresses the 19 substitutions at each position.

Per-variant reliability on these axes is r_full 0.37; per position it is 0.72,
so the position mean is the defensible unit. Individual substitutions are still
drawn as dots, because a position can be high on the mean for two very different
reasons: most substitutions working a little, or a few working a lot.
"""
import sys
sys.path.insert(0, "/Users/mkh/GitHub/mor_dms_analysis/mor_efficacy/src")
sys.path.insert(0, "src")
import numpy as np, pandas as pd, yaml
from plots import figure, save
from common import EVAL, PROC, REPO, core

MARK = {119: "#c00000", 198: "#1b7837"}      # 198 = the literature TM4 serine
lig = yaml.safe_load(open(REPO / "mor_efficacy" / "configs" / "ligands.yaml"))
sc = core(pd.read_parquet(PROC / "scores_long.parquet"))
S = sc.pivot(index="variant_id", columns="ligand", values="score")
AGO = [l for l in S.columns if not lig[l]["is_holdout"] and lig[l]["emax_trupath"]]
LOW = [l for l in AGO if lig[l]["emax_trupath"] <= 70]

d = pd.DataFrame({
    "variant_id": S.index,
    "switch": S[["Naloxone", "Naltrexone"]].mean(1).values - S["FSK"].values,
    "amp": S[LOW].mean(1).values})
d["position"] = [int(v[1:-1]) for v in d.variant_id]

fig = figure(170, 70)
for i, (col, ttl, lab) in enumerate([
        ("switch", "A  antagonist-to-agonist switch",
         "antagonist response minus own basal"),
        ("amp", "B  efficacy amplifier", "gain at weak agonists")]):
    ax = fig.add_axes([0.075 + i * 0.505, 0.20, 0.40, 0.64])
    P = d.groupby("position")[col].agg(["mean", "size"]).reset_index()
    P = P.sort_values("mean", ascending=False).reset_index(drop=True)
    P["rank"] = np.arange(len(P))
    rank_of = dict(zip(P.position, P["rank"]))
    ax.fill_between(P["rank"], 0, P["mean"], color="#d0d0d0", lw=0, zorder=1)
    # individual substitutions behind the position mean
    dd = d.assign(r=d.position.map(rank_of))
    ax.scatter(dd.r, dd[col], s=0.7, color="#9e9e9e", alpha=0.35, lw=0,
               rasterized=True, zorder=2)
    ax.plot(P["rank"], P["mean"], color="#4a4a4a", lw=0.7, zorder=3)
    ax.axhline(0, color="black", lw=0.5, zorder=4)
    for pos, c in MARK.items():
        if pos not in rank_of:
            continue
        r = rank_of[pos]
        sub = dd[dd.position == pos]
        ax.scatter(sub.r, sub[col], s=5, color=c, lw=0, alpha=0.85, zorder=5)
        m = float(P.loc[P.position == pos, "mean"].iloc[0])
        ax.scatter([r], [m], s=26, marker="v", color=c, edgecolor="black",
                   linewidth=0.3, zorder=6)
        ax.annotate(f"{pos}  (#{r + 1} of {len(P)})", (r, m), fontsize=5,
                    color=c, xytext=(7, 9 if pos == 119 else -13),
                    textcoords="offset points")
    ax.set_xlabel(f"positions ranked by mean  (n = {len(P)})")
    ax.set_ylabel(lab, fontsize=5.5)
    ax.set_title(ttl, fontsize=6.5, loc="left")
    ax.set_xlim(-6, len(P) + 6)
fig.text(0.5, 0.03, "grey area/line = position mean; small grey dots = individual "
         "substitutions; triangle = highlighted position mean", ha="center", fontsize=5)
save(fig, EVAL / "fig_waterfall_positions.pdf")
save(fig, EVAL / "fig_waterfall_positions.png")

for col in ["switch", "amp"]:
    P = d.groupby("position")[col].mean().sort_values(ascending=False).reset_index()
    for pos in MARK:
        r = int(P[P.position == pos].index[0]) + 1
        print(f"{col:7s} position {pos}: mean {float(P[P.position==pos][col].iloc[0]):+.3f}  rank {r} of {len(P)}")
    print(f"  top 8: {P.head(8).position.tolist()}")

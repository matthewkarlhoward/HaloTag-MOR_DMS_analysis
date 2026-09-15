"""Identification of A119L: the position waterfall plus every substitution at 119.

Panel A (top-left cell of the grid) is the position-level waterfall on the
antagonist-to-agonist switch axis, which is what nominates position 119 in the
first place. Panels B are the 19 substitutions at 119, sorted by fitted efficacy
slope (most negative = strongest "everything becomes a fuller agonist"
phenotype). Naloxone and naltrexone are drawn separately rather than stacked at
the same x, so that agreement between the two antagonists is visible as a
consistency check.

Basal (no ligand, FSK) sits on its own segment of the x axis, broken away from
the ligand series, rather than being carried across the panel as a level line.

Canvas is 165 x 175 mm, i.e. the final supplemental figure, not a panel to be
assembled elsewhere.
"""
import sys
sys.path.insert(0, "/Users/mkh/GitHub/mor_dms_analysis/mor_efficacy/src")
sys.path.insert(0, "src")
import numpy as np, pandas as pd, yaml
import matplotlib as mpl
from plots import figure, save
from common import PROC, REPO, core

POS = 119
MIN_N = 5            # substitutions a position needs to enter the waterfall
AGREE = 0.05         # naloxone/naltrexone must agree to within this
BASAL_C = "#7b3fa0"
MARK_C = "#c00000"

# --- canvas geometry, millimetres ------------------------------------------
W, H = 165.0, 175.0
NCOL, NROW = 5, 4
LEFT, RIGHT = 11.5, 3.0
BLOCK = (W - LEFT - RIGHT - 3.5 * (NCOL - 1)) / NCOL     # panel block width
PITCH_X = BLOCK + 3.5
W_BAS, W_GAP = 2.4, 1.6                                  # basal stub, axis break
W_MAIN = BLOCK - W_BAS - W_GAP
TOP0, PITCH_Y, PH = 163.0, 37.0, 25.0                    # row-0 top, pitch, height
W_WF = BLOCK - 4.0   # waterfall is narrowed to clear the next panel's y labels


def ax_mm(x, y, w, h):
    """Axes placed by lower-left corner in mm on the 165 x 175 canvas."""
    return fig.add_axes([x / W, y / H, w / W, h / H])


lig = yaml.safe_load(open(REPO / "mor_efficacy" / "configs" / "ligands.yaml"))
sc = core(pd.read_parquet(PROC / "scores_long.parquet"))
S = sc.pivot(index="variant_id", columns="ligand", values="score")
E = sc.pivot(index="variant_id", columns="ligand", values="score_se")
expr = pd.read_parquet(PROC / "expression.parquet").set_index("variant_id")

AGO = sorted([l for l in S.columns if not lig[l]["is_holdout"]
              and lig[l]["emax_trupath"] is not None],
             key=lambda l: lig[l]["emax_trupath"])
X = np.array([lig[l]["emax_trupath"] for l in AGO], float)
Z = (X - X.mean()) / X.std()

vs = [v for v in S.index if v[1:-1] == str(POS)]
slope = {v: np.polyfit(Z, S.loc[v, AGO].values, 1)[0] for v in vs}
vs = sorted(vs, key=lambda v: slope[v])          # most negative slope first

# --- waterfall table (panel A) ---------------------------------------------
wf = pd.DataFrame({
    "switch": S[["Naloxone", "Naltrexone"]].mean(1).values - S["FSK"].values,
    "ok": (S["Naloxone"] - S["Naltrexone"]).abs().values <= AGREE},
    index=S.index)
wf["position"] = [int(v[1:-1]) for v in wf.index]
P = (wf[wf.ok].groupby("position").switch.agg(["mean", "size"]).reset_index()
       .query(f"size >= {MIN_N}")
       .sort_values("mean", ascending=False).reset_index(drop=True))
P["rank"] = np.arange(len(P))

fig = figure(W, H)
cmap = mpl.cm.get_cmap("RdBu_r")

# ============================ panel A: waterfall ===========================
ax = ax_mm(LEFT, TOP0 - PH, W_WF, PH)
ax.fill_between(P["rank"], 0, P["mean"], color="#d4d4d4", lw=0, zorder=1)
ax.plot(P["rank"], P["mean"], color="#3a3a3a", lw=0.7, zorder=2)
ax.axhline(0, color="black", lw=0.5, zorder=3)
r119 = int(P.loc[P.position == POS, "rank"].iloc[0])
m119 = float(P.loc[P.position == POS, "mean"].iloc[0])
ax.scatter([r119], [m119], s=16, marker="v", color=MARK_C, edgecolor="black",
           linewidth=0.3, zorder=5)
ax.annotate(f"position {POS}\n#{r119 + 1} of {len(P)}", (r119, m119), fontsize=5,
            color=MARK_C, xytext=(7, 2), textcoords="offset points")
ax.set_xlim(-4, len(P) + 4)
ax.set_xlabel(f"positions, ranked  (n = {len(P)})", fontsize=5.5, labelpad=1.5)
ax.set_ylabel("antagonist response\nminus own basal", fontsize=5.5, labelpad=1.5)
ax.tick_params(labelsize=5, pad=1.5)
fig.text((LEFT - 9.0) / W, (TOP0 + 1.0) / H, "A", fontsize=8, fontweight="bold",
         va="bottom")

# ====================== panel B: every substitution ========================


def break_mark(x_mm, y_mm):
    """A pair of diagonal slashes sitting in the gap between the two x segments."""
    for dx in (-0.45, 0.45):
        fig.add_artist(mpl.lines.Line2D(
            [(x_mm + dx - 0.35) / W, (x_mm + dx + 0.35) / W],
            [(y_mm - 0.7) / H, (y_mm + 0.7) / H],
            transform=fig.transFigure, color="black", lw=0.5,
            solid_capstyle="butt", zorder=10))


for i, v in enumerate(vs):
    r, c = divmod(i + 1, NCOL)            # cell 0 is the waterfall
    x0 = LEFT + c * PITCH_X
    y0 = TOP0 - r * PITCH_Y - PH
    leftmost = (c == 0) or (i == 0)
    bottom = (r == NROW - 1)

    bax = ax_mm(x0, y0, W_BAS, PH)                        # basal stub
    ax = ax_mm(x0 + W_BAS + W_GAP, y0, W_MAIN, PH)        # ligand series

    col = cmap(0.5 + 0.5 * np.sign(-slope[v]) *
               min(1.0, 0.35 + abs(slope[v]) / 0.13 * 0.65))

    # --- basal segment: no ligand, its own piece of x axis
    fsk, fsk_e = S.loc[v, "FSK"], E.loc[v, "FSK"]
    bax.axhline(0, color="black", lw=0.4, ls=(0, (3, 2)), zorder=1)
    bax.errorbar([0], [fsk], yerr=[fsk_e], fmt="none", ecolor=BASAL_C,
                 elinewidth=0.4, zorder=3)
    bax.plot([0], [fsk], marker="X", ms=3.4, color=BASAL_C, mew=0, zorder=4)
    bax.set_xlim(-1, 1)
    bax.set_xticks([0])
    bax.set_xticklabels(["basal"] if bottom else [], fontsize=4.2)
    bax.spines["right"].set_visible(False)
    bax.tick_params(axis="x", length=1.5, pad=1.2)
    bax.tick_params(axis="y", labelsize=5, pad=1.5)
    if not leftmost:
        bax.set_yticklabels([])
    else:
        bax.set_ylabel("DMS score", fontsize=5.5, labelpad=1.5)

    # --- ligand segment
    ax.axhline(0, color="black", lw=0.4, ls=(0, (3, 2)), zorder=1)
    ax.axvspan(-3, 30, color="#f2f2f2", lw=0, zorder=0)
    ax.errorbar([8], [S.loc[v, "Naloxone"]], yerr=[E.loc[v, "Naloxone"]],
                fmt="o", ms=2.4, color="#444444", mfc="white", mew=0.7,
                elinewidth=0.4, zorder=4)
    ax.errorbar([20], [S.loc[v, "Naltrexone"]], yerr=[E.loc[v, "Naltrexone"]],
                fmt="s", ms=2.4, color="#444444", mfc="white", mew=0.7,
                elinewidth=0.4, zorder=4)
    ax.errorbar(X, S.loc[v, AGO].values, yerr=E.loc[v, AGO].values, fmt="o",
                ms=1.8, color=col, lw=0, elinewidth=0.4, zorder=3)
    y = S.loc[v, AGO].values
    b1, b0 = np.polyfit(Z, y, 1)
    xs = np.linspace(X.min(), X.max(), 20)
    ax.plot(xs, b0 + b1 * (xs - X.mean()) / X.std(), color=col, lw=1.2, zorder=2)
    ax.set_xlim(-3, 106)
    ax.set_xticks([0, 50, 100])
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, labelleft=False)
    ax.tick_params(axis="x", labelsize=5, pad=1.2)
    if not bottom:
        ax.set_xticklabels([])

    for a in (bax, ax):
        a.set_ylim(-0.6, 0.75)
        a.set_yticks([-0.5, 0.0, 0.5])
    break_mark(x0 + W_BAS + W_GAP / 2, y0)

    ex = expr.expression_score.get(v, np.nan)
    ax.set_title(f"{v}   slope {b1:+.3f}", fontsize=5.5, pad=1.5, x=0.44)
    ax.text(0.03, 0.04, f"expr {ex:+.2f}", transform=ax.transAxes, fontsize=4.2,
            color="#666666")
    if i == 0:
        fig.text((x0 - 8.5) / W, (y0 + PH + 1.0) / H, "B", fontsize=8,
                 fontweight="bold", va="bottom")

fig.text(0.5, (H - 6.0) / H, "Identification of A119L: position 119 tops the "
         "antagonist-to-agonist switch, and leucine is the strongest "
         "substitution at the position", ha="center", va="top", fontsize=6.5)

CAP = [
    "A  positions ranked by mean antagonist response minus that variant's own basal. Naloxone and naltrexone are required to agree to within 0.05 and a position needs",
    "\u2265 5 surviving substitutions to enter the ranking; red triangle = position 119.",
    "B  one panel per substitution at position 119, ordered by efficacy slope (most negative first). Coloured points = the agonist series placed by WT ligand efficacy,",
    "solid line = fit to the agonists only. Open circle = naloxone, open square = naltrexone, drawn apart so the two antagonists can be compared.",
    "Purple X on the separate left segment of the x axis = basal, no ligand. Antagonist points above it are ligand-driven; level with it, constitutive.",
]
for j, line in enumerate(CAP):
    fig.text(0.5, (17.0 - j * 3.4) / H, line, ha="center", fontsize=4.5)

XL = "WT ligand efficacy (TRUPATH Emax, % DAMGO)"
fig.text((LEFT + 2 * PITCH_X + W_BAS + W_GAP + W_MAIN / 2) / W,
         (TOP0 - 3 * PITCH_Y - PH - 5.2) / H, XL, ha="center", fontsize=5.5)

save(fig, "eval/figF_position119_all_substitutions.pdf")
save(fig, "eval/figF_position119_all_substitutions.png")

print(f"waterfall: position {POS} rank {r119 + 1} of {len(P)}, mean {m119:+.3f}")
print(f"{len(vs)} substitutions at {POS}, by slope:")
for v in vs:
    print(f"  {v}  slope {slope[v]:+.3f}   Nlx {S.loc[v,'Naloxone']:+.2f}  "
          f"Ntx {S.loc[v,'Naltrexone']:+.2f}  FSK {S.loc[v,'FSK']:+.2f}  "
          f"expr {expr.expression_score.get(v, np.nan):+.2f}")

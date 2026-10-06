#!/usr/bin/env python3
"""
Extract the raw BRET replicates behind the K100 / V175 / R278 TruPath panels.

Those panels (figure "G protein interface", panels d-f) were drawn in Prism from
raw 515/410 BRET, never normalised. The replicates live in two Prism projects on
the Desktop (not in the repo); this script pulls them into a long CSV so the
downstream normalisation no longer needs them:

  A  20260125_mor_ICL_mutants.prism     "<Ligand> read 2 comp"   6 reps
       WT / K100D / K100N / V175E / V175N
       ligands DAMGO, Morphine, Butorphanol, Fentanyl, PZM21, Nalbuphine
       -> panel e (V175E/V175N), and an independent K100 replicate set

  B  20260125_mor_ICL_mutants.prism     "20260427 data: DAMGO R278"   6 reps
       a verbatim copy of the RAW 20260427 doubles plate (the three
       "20260427 data: ..." sheets are byte-identical); 32 sets =
       {WT, 119L, 100D, 100N, 278D, 119L+X} x {DAMGO, PZM21, Nalbuphine,
       Naloxone}. This is what panels d and f actually plot, and figures/
       doubles1_points.csv is its first FOUR replicates after normalisation.

  C  20260121_mkh_mkh_mutant_trupath.prism  "<Ligand> read 2 comp"   6 reps
       WT / R278D / I280D / I280K / I280P, ligands DAMGO, Buprenorphine,
       Morphine, MP, Fentanyl, Nalbuphine -- an earlier, independent R278D run
       (no PZM21), kept for the ligands the 20260427 plate does not cover.

Output: figures/singles_raw_points.csv
  source, run, ligand, variant, logM, rep, bret
Dose 0 (the no-ligand well) is kept and flagged by logM == 0; the -14 point is
kept as measured. Downstream normalisation uses -13 .. -4 only.
"""
import csv
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "figures" / "singles_raw_points.csv"

PRISM_ICL = "/Users/mkh/Desktop/mor_trupath/20260125_mor_ICL_mutants.prism"
PRISM_278 = "/Users/mkh/Desktop/mor_trupath/20260121_mkh_mkh_mutant_trupath.prism"

ICL_LIGANDS = ["DAMGO", "Morphine", "Butorphanol", "Fentanyl", "PZM21",
               "Nalbuphine"]
P278_LIGANDS = ["DAMGO", "Buprenorphine", "Morphine", "MP", "Fentanyl",
                "Nalbuphine"]


def _txt(t):
    """Prism titles are either plain strings or rich-text dicts."""
    if isinstance(t, dict):
        return t.get("string", "")
    return t if isinstance(t, str) else str(t)


def read_prism(path):
    """{sheet title: (set titles, replicatesCount, doses, matrix)}."""
    zf = zipfile.ZipFile(path)
    doc = json.loads(zf.read("document.json"))
    out = {}
    for uid in doc["sheets"]["data"]:
        sh = json.loads(zf.read(f"data/sheets/{uid}/sheet.json"))
        tbl = sh["table"]
        sets = [_txt(json.loads(zf.read(f"data/sets/{ds}.json")).get("title", ""))
                for ds in tbl["dataSets"]]
        raw = zf.read(f"data/tables/{tbl['uid']}/data.csv").decode()
        doses, mat = [], []
        for row in csv.reader(raw.splitlines()):
            if not row or row[0] == "":
                continue
            doses.append(float(row[0]))
            mat.append([float(v) if v not in ("", "nan") else np.nan
                        for v in row[1:]])
        out[_txt(sh.get("title", ""))] = (sets, tbl.get("replicatesCount"),
                                          np.array(doses), np.array(mat))
    return out


def long_rows(sets, reps, doses, mat, split, source, run, keep_ligand=None):
    """Melt one comp sheet. `split` maps a set title -> (ligand, variant)."""
    rows, seen = [], set()
    for i, title in enumerate(sets):
        if not title:
            continue
        ligand, variant = split(title)
        if keep_ligand is not None and ligand != keep_ligand:
            continue                       # reference column from another drug
        if (ligand, variant) in seen:
            continue                       # duplicated reference column
        seen.add((ligand, variant))
        block = mat[:, i * reps:(i + 1) * reps]
        for d, dose in enumerate(doses):
            for r in range(reps):
                y = block[d, r]
                if np.isnan(y):
                    continue
                rows.append(dict(source=source, run=run, ligand=ligand,
                                 variant=variant, logM=dose, rep=r,
                                 bret=float(y)))
    return rows


def main():
    rows = []

    icl = read_prism(PRISM_ICL)
    for lig in ICL_LIGANDS:                                   # source A
        sets, reps, doses, mat = icl[f"{lig} read 2 comp"]
        rows += long_rows(sets, reps, doses, mat,
                          split=lambda t: tuple(t.split(" ", 1)),
                          source="ICL", run="20260125", keep_ligand=lig)

    sets, reps, doses, mat = icl["20260427 data: DAMGO R278"]  # source B
    rows += long_rows(sets, reps, doses, mat,
                      split=lambda t: tuple(t.split("_", 1)),
                      source="doubles_raw", run="20260427")

    p278 = read_prism(PRISM_278)                              # source C
    for lig in P278_LIGANDS:
        sets, reps, doses, mat = p278[f"{lig} read 2 comp"]
        rows += long_rows(sets, reps, doses, mat,
                          split=lambda t: tuple(t.split("-", 1)),
                          source="R278_run1", run="20260121", keep_ligand=lig)

    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    for src, g in df.groupby("source"):
        print(f"{src:12s} run {g.run.iloc[0]}  "
              f"{g.ligand.nunique()} ligands x {g.variant.nunique()} variants, "
              f"{g.rep.nunique()} reps, doses "
              f"{sorted(g.logM.unique())[0]:.0f}..{sorted(g.logM.unique())[-1]:.0f}"
              f"  ({len(g)} rows)")
    print(f"-> {OUT.relative_to(REPO)}  ({len(df)} rows)")


if __name__ == "__main__":
    main()

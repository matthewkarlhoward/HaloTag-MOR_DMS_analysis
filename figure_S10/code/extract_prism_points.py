#!/usr/bin/env python3
"""
Extract the raw normalised dose-response replicates behind doubles1/doubles2.

The Prism fit tables in figures/doubles{1,2}.xlsx were fit to these sheets
(matched by Prism's recorded #Y-values: 40 = 4x10, 54 = 6x9):

  doubles1.xlsx  <-  20260427_doubles.prism : "plates 1-4 norm"  (4 reps, doses -4..-13)
  doubles2.xlsx  <-  20260421_doubles.prism : "comp norm"        (6 reps, doses -5..-13)

Both live on the external drive; edit SOURCES if the mount path differs. Output
is written into the repo so downstream plotting no longer needs the drive:

  figures/doubles1_points.csv , figures/doubles2_points.csv
      long format: set, ligand, variant, partner, is_double, logM, rep, response
"""
import csv
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "figures"

# (output stem, .prism path, data-sheet title)
SOURCES = {
    "doubles1": ("/Volumes/MKH_DRIVE/drive2/20260427/20260427_doubles.prism",
                 "plates 1-4 norm"),
    "doubles2": ("/Volumes/MKH_DRIVE/drive/20260421/20260421_doubles.prism",
                 "comp norm"),
}


def read_prism_sheet(prism_path, sheet_title):
    """Return (variant_names, reps, doses[np], matrix[rows x names*reps])."""
    zf = zipfile.ZipFile(prism_path)
    names = zf.namelist()

    def load(path):
        return json.loads(zf.read(path))

    doc = load("document.json")
    target = None
    for uid in doc["sheets"]["data"]:
        sh = load(f"data/sheets/{uid}/sheet.json")
        if sh["title"] == sheet_title:
            target = sh
            break
    if target is None:
        raise SystemExit(f"sheet {sheet_title!r} not found in {prism_path}")

    tbl = target["table"]
    reps = tbl["replicatesCount"]
    variants = [load(f"data/sets/{ds}.json")["title"] for ds in tbl["dataSets"]]

    csv_bytes = zf.read(f"data/tables/{tbl['uid']}/data.csv").decode()
    doses, mat = [], []
    for row in csv.reader(csv_bytes.splitlines()):
        if not row or row[0] == "":
            continue
        doses.append(float(row[0]))
        mat.append([float(v) if v not in ("", "nan") else np.nan
                    for v in row[1:]])
    return variants, reps, np.array(doses), np.array(mat)


def to_long(variants, reps, doses, mat, set_name):
    rows = []
    for i, variant in enumerate(variants):
        block = mat[:, i * reps:(i + 1) * reps]         # rows(dose) x reps
        ligand, _, tail = variant.partition("_")        # DAMGO_119L_100D
        is_double = "119L_" in variant and variant != f"{ligand}_119L"
        partner = tail.replace("119L_", "") if is_double else tail
        for d, dose in enumerate(doses):
            for r in range(reps):
                y = block[d, r]
                if np.isnan(y):
                    continue
                rows.append(dict(set=set_name, ligand=ligand, variant=variant,
                                 partner=partner, is_double=is_double,
                                 logM=dose, rep=r, response=y))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    for stem, (prism, sheet) in SOURCES.items():
        variants, reps, doses, mat = read_prism_sheet(prism, sheet)
        df = to_long(variants, reps, doses, mat, stem)
        dest = OUT / f"{stem}_points.csv"
        df.to_csv(dest, index=False)
        print(f"{stem}: {len(variants)} variants x {reps} reps x {len(doses)} "
              f"doses -> {dest.relative_to(REPO)}  ({len(df)} rows)")

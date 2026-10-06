#!/usr/bin/env python3
"""Extract the 20260502 MOR WT TruPath run (DAMGO / buprenorphine / SR-17018) into tidy CSVs.

Reads only from raw/ in this folder, so it runs with MKH_DRIVE unplugged.

Outputs (written next to this script):
  sr17018_trupath_wells.csv        per-well donor/acceptor/ratio, both reads
  sr17018_trupath_dose_summary.csv per ligand x dose means, plus donor QC
  sr17018_trupath_prism_values.csv the BRET ratios exactly as stored in the .prism
  sr17018_trupath_prism_fits.csv   Prism's own 3-parameter fit table
  sr17018_trupath_qc.csv           the donor-collapse / donor-ratio-coupling check

Plate map (reverse engineered from the raw reads, then validated against the .prism
value-for-value -- see validate() and README.md):
  rows    DAMGO A/F/K, buprenorphine B/C/G/H/L/M, SR-17018 D/E/I/J/N/O, row P blanks
  columns dose pairs, 1-2 = 1e-4 M down to 21-22 = 1e-14 M, 23-24 = vehicle
  the .prism analysis uses READ2 only; READ1 is the earlier read of the same plate
"""

import csv
import io
import json
import re
import statistics as st
import zipfile
from pathlib import Path

import openpyxl

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
PRISM = RAW / "20260502_mor_WT_DAMGO_SR_BUP.prism"
READS = {"READ1": RAW / "20260502_READ1.XLSX", "READ2": RAW / "20260502_READ2.xlsx"}

LIGAND_ROWS = {"DAMGO": "AFK", "Buprenorphine": "BCGHLM", "SR17018": "DEIJNO"}
ROW_LIGAND = {r: lig for lig, rows in LIGAND_ROWS.items() for r in rows}
# log10[ligand] -> the column pair holding it
DOSE_COLS = {
    -4: (1, 2), -5: (3, 4), -6: (5, 6), -7: (7, 8), -8: (9, 10), -9: (11, 12),
    -10: (13, 14), -11: (15, 16), -12: (17, 18), -13: (19, 20), -14: (21, 22),
    0: (23, 24),  # vehicle; encoded as 0 to match the .prism x column
}
COL_DOSE = {c: d for d, cols in DOSE_COLS.items() for c in cols}
PRISM_LIGAND_ORDER = ["DAMGO", "Buprenorphine", "SR17018"]
PRISM_READ = "READ2"
BLANK_ROW = "P"


def read_plate(path):
    """Return {'donor': {row: [24 vals]}, 'acceptor': {...}} from a CLARIOstar export.

    The two exports differ in layout: one starts row labels in column A, the other in
    column B, so the label column is detected per row rather than assumed.
    """
    ws = openpyxl.load_workbook(path, read_only=True, data_only=True).worksheets[0]
    out, channel = {}, None
    for row in ws.iter_rows(values_only=True):
        vals = list(row)
        joined = " ".join(str(v) for v in vals if v is not None)
        if "Raw Data (410" in joined:
            channel = "donor"
            out[channel] = {}
            continue
        if "Raw Data (515" in joined:
            channel = "acceptor"
            out[channel] = {}
            continue
        if channel is None:
            continue
        for j, v in enumerate(vals[:3]):
            if isinstance(v, str) and len(v.strip()) == 1 and v.strip().isalpha():
                out[channel][v.strip()] = [
                    x for x in vals[j + 1:] if isinstance(x, (int, float))
                ]
                break
    return out


def prism_data_table():
    """BRET ratios as stored in the .prism 'Data 1' sheet -> {(ligand, x): [values]}."""
    z = zipfile.ZipFile(PRISM)
    uid = None
    for name in z.namelist():
        if name.startswith("data/sheets/") and name.endswith("sheet.json"):
            sheet = json.loads(z.read(name).decode())
            if sheet.get("title") == "Data 1":
                uid = sheet["table"]["uid"]
                break
    if uid is None:
        raise RuntimeError("no 'Data 1' sheet in the .prism file")
    rows = list(csv.reader(io.StringIO(z.read(f"data/tables/{uid}/data.csv").decode())))
    table = {}
    for row in rows:
        if not row or not row[0].strip():
            continue
        x = float(row[0])
        cells = row[1:]
        for i, ligand in enumerate(PRISM_LIGAND_ORDER):
            vals = [float(c) for c in cells[i * 12:(i + 1) * 12] if c.strip()]
            if vals:
                table.setdefault((ligand, x), []).extend(vals)
    return table


def prism_fit_table():
    """Prism's 'Table of results' rows for the 3-parameter fit."""
    z = zipfile.ZipFile(PRISM)
    for name in z.namelist():
        if name.startswith("data/sheets/") and name.endswith("sheet.json"):
            sheet = json.loads(z.read(name).decode())
            if "results" in str(sheet.get("title", "")).lower():
                uid = sheet["table"]["uid"]
                text = z.read(f"data/tables/{uid}/data.csv").decode()
                return [r for r in csv.reader(io.StringIO(text)) if any(c.strip() for c in r)]
    return []


def wells(plates):
    """Flatten both reads into per-well records."""
    out = []
    for read, plate in plates.items():
        for row in sorted(plate["donor"]):
            donor, acceptor = plate["donor"][row], plate["acceptor"][row]
            for i, (dv, av) in enumerate(zip(donor, acceptor), start=1):
                out.append({
                    "read": read,
                    "plate_row": row,
                    "plate_col": i,
                    "ligand": "blank" if row == BLANK_ROW else ROW_LIGAND.get(row, ""),
                    "log_dose_M": "" if row == BLANK_ROW else COL_DOSE.get(i, ""),
                    "is_vehicle": int(row != BLANK_ROW and COL_DOSE.get(i) == 0),
                    "donor_410": dv,
                    "acceptor_515": av,
                    "bret_ratio": round(av / dv, 6) if dv else "",
                    "used_in_prism_fit": int(read == PRISM_READ and row != BLANK_ROW),
                })
    return out


def validate(plates, prism):
    """Check the plate map reproduces the .prism values exactly. Returns report rows."""
    plate = plates[PRISM_READ]
    report = []
    for (ligand, x), pv in sorted(prism.items()):
        expected = sorted(round(
            plate["acceptor"][r][c - 1] / plate["donor"][r][c - 1], 6)
            for r in LIGAND_ROWS[ligand] for c in DOSE_COLS[x])
        got = sorted(round(v, 6) for v in pv)
        exact = len(got) == len(expected) and all(
            abs(a - b) < 1e-4 for a, b in zip(got, expected))
        actual = x if exact else None
        if not exact:
            # find which raw dose (or pair of doses) the prism row really holds
            for cand in DOSE_COLS:
                ref = sorted(round(
                    plate["acceptor"][r][c - 1] / plate["donor"][r][c - 1], 6)
                    for r in LIGAND_ROWS[ligand] for c in DOSE_COLS[cand])
                if len(ref) == len(got) and all(abs(a - b) < 1e-5 for a, b in zip(got, ref)):
                    actual = cand
                    break
            else:
                for c1 in DOSE_COLS:
                    for c2 in DOSE_COLS:
                        if c1 >= c2:
                            continue
                        ref = sorted(round(
                            plate["acceptor"][r][c - 1] / plate["donor"][r][c - 1], 6)
                            for cand in (c1, c2)
                            for r in LIGAND_ROWS[ligand] for c in DOSE_COLS[cand])
                        if len(ref) == len(got) and all(
                                abs(a - b) < 1e-5 for a, b in zip(got, ref)):
                            actual = f"{c1}+{c2}"
                            break
                    if actual is not None:
                        break
        report.append({
            "ligand": ligand,
            "prism_log_dose_label": x,
            "n": len(pv),
            "matches_label": int(exact),
            "actual_raw_log_dose": actual if actual is not None else "unresolved",
        })
    return report


def main():
    plates = {k: read_plate(v) for k, v in READS.items()}
    prism = prism_data_table()

    well_rows = wells(plates)
    write_csv("sr17018_trupath_wells.csv", well_rows)

    # dose summary from the validated raw mapping (correct dose labels for all three)
    summary = []
    plate = plates[PRISM_READ]
    for ligand, rows in LIGAND_ROWS.items():
        for dose in sorted(DOSE_COLS, reverse=True):
            ratios, donors = [], []
            for r in rows:
                for c in DOSE_COLS[dose]:
                    dv, av = plate["donor"][r][c - 1], plate["acceptor"][r][c - 1]
                    ratios.append(av / dv)
                    donors.append(dv)
            summary.append({
                "ligand": ligand,
                "log_dose_M": dose,
                "n": len(ratios),
                "bret_mean": round(st.mean(ratios), 6),
                "bret_sd": round(st.stdev(ratios), 6),
                "bret_sem": round(st.stdev(ratios) / len(ratios) ** 0.5, 6),
                "donor_410_median": round(st.median(donors), 1),
            })
    write_csv("sr17018_trupath_dose_summary.csv", summary)

    prism_rows = [
        {"ligand": lig, "prism_log_dose_label": x, "replicate": i + 1,
         "bret_ratio": round(v, 6)}
        for (lig, x), vals in sorted(prism.items()) for i, v in enumerate(vals)
    ]
    write_csv("sr17018_trupath_prism_values.csv", prism_rows)

    with open(HERE / "sr17018_trupath_prism_fits.csv", "w", newline="") as fh:
        csv.writer(fh).writerows(prism_fit_table())

    write_csv("sr17018_trupath_qc.csv", qc(plate))

    val = validate(plates, prism)
    write_csv("sr17018_trupath_prism_label_check.csv", val)
    bad = [r for r in val if not r["matches_label"]]
    print(f"plate map validated: {len(val) - len(bad)}/{len(val)} prism rows match their label")
    for r in bad:
        print(f"  MISLABELLED  {r['ligand']} labelled 1e{r['prism_log_dose_label']:.0f}"
              f" is really raw dose 1e{r['actual_raw_log_dose']}")


def qc(plate):
    """Donor collapse and donor/ratio coupling, per ligand and dose band."""
    bands = {"top_dose": (1, 2), "mid_dose": (11, 12), "low_dose": (19, 20),
             "vehicle": (23, 24)}
    out = []
    for ligand, rows in LIGAND_ROWS.items():
        for label, cols in bands.items():
            donors, ratios = [], []
            for r in rows:
                for c in cols:
                    dv, av = plate["donor"][r][c - 1], plate["acceptor"][r][c - 1]
                    donors.append(dv)
                    ratios.append(av / dv)
            out.append({
                "ligand": ligand,
                "band": label,
                "plate_cols": f"{cols[0]}-{cols[1]}",
                "n": len(donors),
                "bret_mean": round(st.mean(ratios), 6),
                "donor_410_median": round(st.median(donors), 1),
                "r_donor_vs_ratio": round(pearson(donors, ratios), 3),
            })
    # plate-wide reference: is ratio coupled to donor at all on this plate?
    donors, ratios = [], []
    for r in plate["donor"]:
        if r == BLANK_ROW:
            continue
        for dv, av in zip(plate["donor"][r], plate["acceptor"][r]):
            if dv > 200:
                donors.append(dv)
                ratios.append(av / dv)
    out.append({
        "ligand": "ALL", "band": "plate_wide", "plate_cols": "1-24",
        "n": len(donors), "bret_mean": round(st.mean(ratios), 6),
        "donor_410_median": round(st.median(donors), 1),
        "r_donor_vs_ratio": round(pearson(donors, ratios), 3),
    })
    return out


def pearson(a, b):
    ma, mb = st.mean(a), st.mean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** 0.5
    return num / den if den else float("nan")


def write_csv(name, rows):
    if not rows:
        return
    with open(HERE / name, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {name}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()

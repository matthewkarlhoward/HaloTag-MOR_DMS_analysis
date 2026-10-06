# Archiving this repository to Zenodo

The plan: a **standalone Zenodo record**, not linked to GitHub. The DOI is reserved
now so it can go in the preprint, the record stays private until the preprint posts,
then it is published and opened. GitHub gets its own separate link in the preprint.

> Zenodo's interface has changed since the InvenioRDM migration, so button labels may
> not match exactly. The concepts (reserve DOI, access level, versioning) are stable.

---

## Why not the GitHub integration

Zenodo's GitHub integration only lists and archives **public** repositories, and it is
driven by a release webhook it cannot receive for a private repo. Since this repo stays
private until the preprint is up, the integration is not usable for the preprint DOI.

Note also that a manual record and a GitHub-integration record are **separate records
with separate concept DOIs**. A manual upload cannot later be folded into a GitHub
integration version chain, so the route chosen here is the one to stay on.

---

## Step 1 — upload artifact

Already built from the `v1.0.0` tag:

```
HaloTag-MOR_DMS-v1.1.0.zip     ~205 MB, ~1,410 files
```

Rebuild it any time with:

```bash
git archive --format=zip --prefix=HaloTag-MOR_DMS-v1.1.0/ -o HaloTag-MOR_DMS-v1.1.0.zip v1.1.0
```

Well under Zenodo's 50 GB per-record limit.

## Step 2 — create the record and reserve the DOI

1. Sign in at <https://zenodo.org> (link ORCID while you are there; it propagates to the
   record and to your ORCID profile).
2. **New upload**, drag in the zip. Wait for the upload to finish before anything else.
3. Set **Upload type** to *Dataset*. This repo is code *and* data; *Dataset* is the
   better fit given `data/` and the supplementary workbook dominate. *Software* is
   defensible if you would rather emphasise the code.
4. Find the DOI field and choose **"Reserve DOI"** (sometimes shown as *Get a DOI now!*).
   This is the step that matters — it hands you `10.5281/zenodo.22051140` immediately,
   without publishing anything.
5. Set **Access** to **Restricted**, and leave the record as a **saved draft**. Do not
   press Publish yet.
6. **Save.**

You now have a DOI string for the preprint, and nothing is publicly visible.

## Step 3 — metadata to paste

`.zenodo.json` in this repo holds the intended metadata, but Zenodo **only reads that
file through the GitHub integration** — for a manual upload you must enter it by hand.
The values are there to copy from. Summary:

| Field | Value |
|---|---|
| Title | Code and data for: Distinct activation mechanisms underlie ligand efficacy at a GPCR |
| Upload type | Dataset |
| License | MIT (code) — note the CC BY 4.0 data split in the description |
| Authors | the 15 authors, in paper order, from `CITATION.cff` |
| Keywords | deep mutational scanning; GPCR; mu-opioid receptor; OPRM1; allostery; pharmacological efficacy; opioids |
| Related identifier | `https://github.com/matthewkarlhoward/HaloTag-MOR_DMS_analysis` — relation *is supplement to* |

Once the preprint has its own DOI, add it as a second related identifier with relation
*is supplement to*.

## Step 4 — in the preprint

Two separate links, since Zenodo is not tracking GitHub:

```
Code and data are archived at Zenodo (DOI: 10.5281/zenodo.22051140) and are also
available at https://github.com/matthewkarlhoward/HaloTag-MOR_DMS_analysis.
```

**The reserved DOI does not resolve until you publish the record.** Anyone clicking it
before Step 5 gets an error page. Publish the record at or before the moment the
preprint goes live.

## Step 5 — when the preprint posts

1. On Zenodo, open the draft and **Publish**. The DOI goes live.
2. Change **Access** from Restricted to **Open**.
3. Make the GitHub repo public:
   ```bash
   gh repo edit matthewkarlhoward/HaloTag-MOR_DMS_analysis --visibility public
   ```
4. Add the DOI badge to `README.md` (placeholder already in place — swap the number in
   both spots):
   ```markdown
   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22051140.svg)](https://doi.org/10.5281/zenodo.22051140)
   ```

Direction matters: opening a restricted record is something you can do yourself;
**closing an open one requires contacting Zenodo support.** Do not publish as Open
early.

## Updating later

Published files are immutable. To ship a revision, use **New version** on the record —
it mints a new version DOI under the same **concept DOI**. Cite the concept DOI in the
paper so the citation always resolves to the latest version.

If you tag a new release here first, rebuild the zip from that tag so the archive and
the repository stay in step.

## Still outstanding

Sheet `07_validation_pharmacology` of the supplementary workbook covers the double
mutants only. The single-mutant BRET behind figures 7d–f and the flow cytometry behind
supplemental figure 11b are still Prism-only. Worth resolving **before** publishing the
Zenodo record, since files cannot be edited afterwards without a new version.

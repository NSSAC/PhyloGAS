# Mirroring the synthetic populations to Zenodo

> **Provenance record.** The scripts described here are operational one-offs
> and are not tracked in this repository (see `.gitignore`); they live on the
> cluster where they run. This document is kept so the reasoning behind the
> generated metadata survives.
>
> **Status: Virginia is published** at
> [10.5281/zenodo.23088534](https://doi.org/10.5281/zenodo.23088534) with all
> 13 files. Verified post-publication: 12 creators, CC-BY-4.0, CDC grant
> CK22-2204, and all four related identifiers landed correctly.

The UVA Dataverse has been unreliable (website 502, metadata API 500, file API
roughly 50% 503 as of 2026-09-30). Zenodo has been solid for the EpiHiper
replicates, so the synthetic populations are being mirrored there too.

Two steps: **stage the files on the cluster**, then **run the publisher**.

---

## Licensing

The Dataverse deposit is **CC BY 4.0**. Three independent places in the
captured landing page agree:

| source | value |
|---|---|
| schema.org JSON-LD `license` | `http://creativecommons.org/licenses/by/4.0` |
| rendered page text | "Creative Commons Attribution 4.0 International License." |
| badge / link | `CC BY 4.0`, `licenses/by/4.0` |

No custom Terms of Use, Terms of Access, restrictions, or special permissions
are attached — the only matches for those phrases are Dataverse's own site
boilerplate and its unused "request access" widget.

**So the mirror is clearly permitted**, and CC-BY-4.0 on the Zenodo records is
the correct, matching choice.

CC BY's one obligation is attribution. Two things satisfy it:

1. **The 12 original authors are carried over verbatim** as the Zenodo
   creators — Mortveit, Adiga, Baek, Bhattacharya, Eubank, Machi, Marathe,
   Porebski, Swarup, Venkatramanan, Wilson, Xie.
2. **The description states that the record is a mirror** and asks users to
   cite the original Dataverse deposit as well.

One gap worth closing: pass `--dataverse-doi` so the source deposit is recorded
as a machine-readable `isDerivedFrom` relation rather than only as prose. The
Virginia DOI is `10.18130/V3/5LSDCY`; the other five are on their respective
Dataverse pages.

```bash
python scripts/zenodo/publish_populations.py --states ga \
    --base-dir /project/... --dataverse-doi 10.18130/V3/<GA-ID>
```

Worth noting the authorship distinction, since it affects citation: the
**populations** are authored by the Biocomplexity synthetic-population group
and derived from Dataverse, while the **EpiHiper replicates**
(10.5281/zenodo.23067670) are authored by you. Different provenance, which is
why they are separate records rather than one combined deposit.

---

## Part 1 — Staging (for whoever has cluster access)

### Layout

One directory per state code, under a common root:

```
/project/bii_nssac/biosurveillance/dataverse/
├── va/
├── ca/
├── ga/
├── ma/
├── mn/
└── wa/
```

Virginia is already staged at `/project/bii_nssac/biosurveillance/dataverse/va`.
Follow that pattern for the rest.

### Files per state

Thirteen files, matching the Dataverse deposit. `{st}` is the lowercase code:

| file | needed by PhyloGAS? |
|---|---|
| `{st}_persontrait_epihiper.txt.xz` | **yes** — demographics join |
| `{st}_person.csv.xz` | **yes** — race / hispanic |
| `{st}_household.csv.xz` | **yes** — links person to residence |
| `{st}_residence_locations.csv.xz` | **yes** — the only source of lat/long |
| `{st}_contact_network_epihiper.txt.xz` | no — EpiHiper input only |
| `{st}_persontrait_epihiper_db.tar.gz` | no — EpiHiper input only |
| `{st}_population_network.csv.xz` | no |
| `{st}_activity_assignment_day.csv.xz` | no |
| `{st}_location_assignment_day.csv.xz` | no |
| `{st}_activity_locations.csv.xz` | no |
| `{st}_location_info.csv.xz` | no |
| `activity_types.csv` | no — shared lookup |
| `DP_US-Populations_v_2_4_0_...pdf` | no — data dictionary |

The first four are what PhyloGAS actually reads. The rest are mirrored for
completeness so the Zenodo record stands on its own.

Get the exact list for a state:

```bash
python scripts/zenodo/state_metadata.py --state ga --list-files
```

### Notes

* **Keep the `.xz`.** Nothing needs decompressing; pandas reads `.xz` natively
  and the compressed files are ~6x smaller to transfer.
* **Partial staging is fine.** The publisher uploads what it finds and lists
  what is missing, so you can stage incrementally.
* `activity_types.csv` and the PDF are identical across states — copy the same
  ones into each directory.
* **One caveat on `activity_types.csv`:** the Dataverse description says it
  documents v2.4.0 populations, but some of these are v2.5.0. Worth fixing the
  wording before publishing; the text lives in
  `scripts/zenodo/state_metadata.py` under `FILE_DESCRIPTIONS`.

---

## Part 2 — Publishing

### Setup

```bash
# https://zenodo.org/account/settings/applications/tokens/new/
# scopes: deposit:write, deposit:actions
export ZENODO_TOKEN=...
```

### Rehearse first

```bash
# 1. what would be uploaded, and from where
python scripts/zenodo/publish_populations.py \
    --states va --base-dir /project/bii_nssac/biosurveillance/dataverse --dry-run

# 2. a real run against the sandbox (separate account and token)
export ZENODO_SANDBOX_TOKEN=...
python scripts/zenodo/publish_populations.py \
    --states va --base-dir /project/... --sandbox
```

The sandbox is a full copy of Zenodo that mints fake DOIs. Worth one pass
before touching the live site.

### Live

```bash
python scripts/zenodo/publish_populations.py \
    --states va ga ma mn wa ca \
    --base-dir /project/bii_nssac/biosurveillance/dataverse
```

This creates **drafts**. Nothing is published, no DOI is minted. Review each in
the web UI, then publish by hand.

### If it is interrupted

Re-run the same command. Progress is journalled to
`.zenodo_publish_state.json`: records already created are reused, files already
uploaded and checksum-verified are skipped. Safe on a flaky connection, which
matters for the 400 MB-4 GB contact networks.

### Flags

| flag | effect |
|---|---|
| `--dry-run` | report only; create nothing |
| `--sandbox` | target sandbox.zenodo.org (uses `ZENODO_SANDBOX_TOKEN`) |
| `--keep-going` | continue past a failed file or state |
| `--dataverse-doi` | record the original deposit as `isDerivedFrom` |
| `--publication-date` | override (defaults to today) |

Recording provenance is worth doing:

```bash
python scripts/zenodo/publish_populations.py --states va \
    --base-dir /project/... --dataverse-doi 10.18130/V3/5LSDCY
```

---

## What the records will contain

Generated by `scripts/zenodo/state_metadata.py`, with provenance:

| field | source |
|---|---|
| title | `Synthetic Population for {State}, US (ver. 2.4.0)` |
| description | the Dataverse JSON-LD, plus a note explaining the mirror |
| creators | the 12 Dataverse authors, verbatim |
| license | CC-BY-4.0 |
| funding | CDC PGCoE CK22-2204, from the existing PhyloGAS ABM record |
| related | cites the two PhyloGAS papers; `isSupplementTo` the ABM deposit and the GitHub repo |
| per-file descriptions | `dataverse_files.txt` |

Field-for-field parity with the published ABM record (Zenodo 23067670) was
checked, and every vocabulary ID verified live:

```
resource_type  dataset        -> 200
license        cc-by-4.0      -> 200
funder         042twtr12      -> 200  (CDC)
relations      cites, issupplementto, isderivedfrom -> all present
```

Inspect before running:

```bash
python scripts/zenodo/state_metadata.py --state ga | jq .metadata.title
python scripts/zenodo/state_metadata.py --all --outdir /tmp/zmeta
```

---

## After publishing

Each record gets its own DOI. Add them to
`src/phylogas/sources.py` so `phylogas fetch-data` prefers Zenodo over
Dataverse. The file-ID table there is currently Dataverse-only; mirroring will
let it fall back, or switch over entirely.

Tell me the DOIs and I will wire them in.

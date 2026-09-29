# Data acquisition

> **Status: implemented.** The evaluation below led to
> `phylogas fetch-data` and `phylogas build-demographics`. The original
> `download_pgcoe_dataverse-v2.sh` remains as reference for its file-ID table.

Assessed 2026-09-28 against the live UVA Dataverse API. Every file ID in the
script was probed; a partial download was performed and verified.

---

## Verdict

**The script is accurate, works, and should become the basis of
`phylogas fetch-data`.** Its file-ID table is the single most valuable thing in
it — that mapping is not derivable from the DOI alone and would be tedious to
reconstruct.

Three things must change before it is wired in: the **per-state DOI discovery**,
the **`unxz -k` disk blowup**, and the fact that it stops one step short of the
file the pipeline actually reads.

---

## What was verified

| Check | Result |
|---|---|
| API endpoint reachable | Yes — `HTTP/2 303` redirect to S3 |
| File IDs resolve | Yes, all probed IDs |
| ID -> filename mapping correct | **5/5 spot-checks exact** (`va_household.csv.xz`, `va_person.csv.xz`, `wa_persontrait_epihiper.txt.xz`, `mn_person.csv.xz`, `ga_contact_network_epihiper.txt.xz`) |
| Real download works | Yes — 2 MB byte-range fetch returned valid `XZ compressed data, checksum CRC64` |
| `"${@:-${DEFAULT_STATES[@]}}"` idiom | Correct: 0 args -> all 6 states; N args -> those N |
| `set -euo pipefail`, `curl -f --retry 3`, skip-if-exists | All sound |

### Important discovery: the data spans **six DOIs**, not one

`config.template.yaml` and the README both name a single
`doi:10.18130/V3/5LSDCY`. That is only Virginia. Resolved from the redirects:

| State | DOI |
|---|---|
| VA | `10.18130/V3/5LSDCY` |
| CA | `10.18130/V3/A7DQWM` |
| GA | `10.18130/V3/4DWNRH` |
| MA | `10.18130/V3/ZB0SGL` |
| MN | `10.18130/V3/SB2PWT` |
| WA | `10.18130/V3/PNGMRJ` |

**This is why the hardcoded ID table matters** — and also why the config's
single-DOI assumption needs fixing.

---

## Problems to fix before adopting

### 1. `unxz -f -k` roughly doubles an already large footprint

`-k` keeps the compressed original. Actual sizes from the Dataverse API:

| VA file | Compressed |
|---|---|
| `va_contact_network_epihiper.txt.xz` | 907.2 MB |
| `va_persontrait_epihiper.txt.xz` | 25.5 MB |
| `va_person.csv.xz` | 25.4 MB |
| `va_household.csv.xz` | 20.2 MB |
| **total** | **978 MB** |

Decompressed at ~6x, keeping both: **~6.8 GB per state, ~40 GB for all six.**

Worse, the default is *all six states*, so a bare `./download_pgcoe_dataverse-v2.sh`
pulls ~40 GB. Recommend: drop `-k` (or make it a flag), and require an explicit
state list rather than defaulting to everything.

### 2. The 907 MB contact network is downloaded but never used

Nothing in PhyloGAS reads `*_contact_network_epihiper.txt`. It is an EpiHiper
*input*, and PhyloGAS consumes EpiHiper *output* (`output.csv.gz`). It is 93% of
the per-state download.

Make it opt-in (`--with-contact-network`) for users who intend to run EpiHiper
themselves. That alone cuts the default from ~40 GB to ~2.8 GB.

### 3. It stops one step short of what the pipeline reads

The script produces `va_persontrait_epihiper.txt` and `va_person.csv`.
The config wants:

```yaml
population:
  persontrait_file: "{data_dir}/va_2_4_0_demographics.csv"
```

Those are not the same file. The gap is bridged by
`src/phylogas/popprep/merge_persontrait.py` (rescued from
`synthetic_biosurveillance`):

```
va_persontrait_epihiper.txt  ─┐
va_person.csv                 ├─> merge_persontrait.py ─> va_2_4_0_demographics.csv
county FIPS lookup           ─┘
```

The script fetches exactly the two Dataverse inputs this needs — good — but the
merge is not run, so the pipeline's declared input never materialises.
`phylogas fetch-data` should invoke the merge as a post-step.

### 4. Missing FIPS lookup

`merge_persontrait.py` needs a `fips_lookup_file`, and
`popprep/county_centroids.py` reads a hardcoded `state_county_centroids.csv`
that **is not in this repo**. `TwinSampler/Data/county_fips.csv` (57 KB) looks
like the right artifact.

Small file — commit it to `data/` rather than downloading it.

### 5. Minor

- `${state_lc^^}` is a bashism; fine given the `#!/usr/bin/env bash` shebang,
  but it will break under `sh`.
- No checksum verification. Dataverse exposes MD5 per file in the dataset JSON;
  worth checking since a truncated 900 MB download fails late and confusingly.
- Unquoted `for file_name in $files` relies on word splitting. Works here (no
  spaces in names) but an array would be safer.

---

## Remaining data acquisition gaps

With the script adopted and the merge wired in, here is the full picture.

### Resolved by the script (+ merge step)

| Input | Source |
|---|---|
| `population.persontrait_file` | Dataverse + `merge_persontrait.py` |
| `population.household_file` | Dataverse (`*_household.csv.xz`) |

### Still unresolved

| Config key | What it is | Where it comes from | Status |
|---|---|---|---|
| `population.rucc_file` | `Ruralurbancontinuumcodes2023.csv` | **USDA ERS**, not Dataverse | Needs its own fetch; small, could be committed |
| `genetic_painter.align_fasta` | `clean_va_delta_sequences.fasta` — the training MSA | GISAID-derived | **Cannot be redistributed.** Needs a documented recipe via `phylogas prep-seeds`, or ship only the derived entropy outputs |
| `genetic_painter.entropy_thresholds` | `run.03.threshold.file` | Derived from the MSA above | Ship pre-computed — it is small (318 KB) and avoids the GISAID problem entirely |
| `genetic_painter.probability_matrix` | `run.03.base.threshold.df.npy` | Same | Ship pre-computed (3.6 MB) |
| `genetic_painter.seed_fasta` | `Virginia_B_1_617_2_seed_sequences.fasta.gz` | GISAID-derived (13 MB) | Redistribution unclear; `prep-seeds` can rebuild from Cov-Spectrum |
| `genetic_painter.reference_fasta` | `reference.fasta` (Wuhan-Hu-1) | Public | **Commit it** — 30 KB |
| `epihiper.output_csv` | Simulation replicates (~550 MB each) | Your pending deposit | Pending |
| `ascertainment.parameters` | `ascertainment_parameters.yaml` | TwinSampler repo | Ships with that package |
| `ascertainment.schedule_input` | `Virginia_importation_schedule.csv` | TwinSampler `Data/` | Ships with that package |

### The GISAID constraint is the real blocker

Three inputs derive from GISAID sequences that generally cannot be
redistributed. The practical resolution:

1. **Ship the derived artifacts, not the sequences.** `run.03.threshold.file`
   (318 KB) and `run.03.base.threshold.df.npy` (3.6 MB) are aggregate entropy
   statistics, not sequences. Committing these lets users run
   `phylogas paint` without ever touching GISAID, and `phylogas train` becomes
   optional for anyone reproducing published results.
2. **Document the recipe** for users who want to retrain on their own MSA.
3. **`prep-seeds` already pulls from Cov-Spectrum**, which has different terms
   than GISAID — worth confirming what that permits.

This is the single highest-value decision for external adoption: with the two
derived files committed, a new user needs only the EpiHiper replicate plus the
demographics to run the pipeline end to end.

---

## Recommended shape for `phylogas fetch-data`

```
phylogas fetch-data --config config.yaml --states va
    [--with-contact-network]   # opt in to the 907 MB EpiHiper input
    [--keep-compressed]        # opt in to unxz -k behaviour
    [--no-merge]               # skip building the demographics CSV
```

Steps:
1. Resolve file IDs (keep the table; optionally verify against the per-state
   DOI via the datasets API).
2. Download with `curl -f --retry`, skip if present, verify MD5 from the
   dataset JSON.
3. Decompress **without** `-k` by default.
4. Run `merge_persontrait.py` to produce `<state>_2_4_0_demographics.csv`.
5. Print a `validate-config`-style summary of what is now present and what is
   still missing.

Implementation note: keep the hardcoded ID table. Discovering IDs by querying
each DOI is more elegant but adds six API calls and a failure mode, and the IDs
are stable. Record the per-state DOI next to each ID group as documentation.

### Also update

- `config.template.yaml` and `README.md` currently advertise one DOI. Replace
  with the six-DOI table, or point at the PGCOE collection.
- `cli.py`'s `fetch-data` default `doi:10.18130/V3/5LSDCY` is Virginia-only.


---

## Resolution (2026-09-28)

### Contact network: confirmed unused, now opt-in

Searched all three repos (`PhyloGAS`, `TwinSampler`, `BeyondBaseline`) across
`.py`, `.sh`, `.yaml`, `.json`, `.md` and the EpiHiper configs in `cfg/exp1/`:
**zero references to `contact_network`.** It is an EpiHiper *input*; PhyloGAS
consumes EpiHiper *output* (`output.csv.gz`).

It is also 88% of the download:

| | all six states |
|---|---|
| default (`fetch-data`) | **0.91 GB** |
| with `--with-epihiper-inputs` | 8.84 GB |

### The demographics file is a derived join

Confirmed by inspecting the live files. The v2.4.0 persontrait:

* opens with an EpiHiper JSON schema line,
* has `county_fips`, not `county`,
* has `hispanic`, not `latino`,
* **has no `home_latitude`/`home_longitude` at all.**

v1.9.0 carried the coordinates directly, which is why the original
`merge_persontrait_with_person.py` never joined for them. Under v2.4.0 they are
two hops away:

```
persontrait --hid--> household --rlid--> residence_locations(latitude, longitude)
```

This added `<state>_residence_locations.csv.xz` to the file table (missing from
the original script). Verified end to end on Minnesota: 5,218,000 rows,
**100% coordinate coverage**, values inside the state's real bounding box.

### New file IDs

| state | residence_locations | size |
|---|---|---|
| va | 120568 | 32.9 MB |
| ca | 120621 | 116.7 MB |
| ga | 120596 | 39.8 MB |
| ma | 120590 | 22.4 MB |
| mn | 120580 | 27.2 MB |
| wa | 120606 | 30.2 MB |

### Decompression

`.xz` is now kept by default. pandas reads it natively, and these files expand
~6x. `--decompress` is available for tools that need plain text.

### Painter no longer fails silently

Given a raw persontrait file the painter previously caught the `KeyError`,
disabled augmentation, and produced a complete-looking run with empty
demographic columns. It now detects the JSON header, validates the requested
`--add_metadata` columns up front, and exits with the exact command needed to
build the right file.

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

The table below is the running list; rows marked **Resolved** are kept
for the record rather than deleted.

| Config key | What it is | Where it comes from | Status |
|---|---|---|---|
| `population.rucc_file` | `Ruralurbancontinuumcodes2023.csv` | **USDA ERS**, not Dataverse | **Resolved.** A pinned copy is committed at `data/`; `fetch-data --with-rucc` refreshes it from the agency. `fetch-data` only reports it as outstanding when no copy is found |
| `genetic_painter.align_fasta` | the training MSA, under `data/training_sequences/` | Cov-Spectrum, which serves **open** sequences | **Resolved.** `fetch-data --with-training-sequences` runs `seq_prep.py` bulk mode over `training.date_from`/`date_to` (default: the painter's own window) and writes this path. No redistribution constraint — we store no sequences ourselves and the query returns public data. Left out of git only because it is large and exactly reproducible from the configured window |
| `genetic_painter.entropy_thresholds` | stage 00 under `results_dir` | Written by `phylogas train` | Output, not input. A pre-trained map could also be shipped in `data/example_data/` and pointed at — worth doing as a convenience, not as a licensing workaround |
| `genetic_painter.probability_matrix` | stage 00 under `results_dir` | Same | Same; 3.6 MB for the Virginia Delta map |
| `genetic_painter.seed_fasta` | `Virginia_B_1_617_2_seed_sequences.fasta` in `seeds.output_folder` | Cov-Spectrum (open sequences) | **Resolved.** `prep-seeds --seed-mode` builds it; the CLI derives the filename from state + lineage and finds it without a config edit |
| `genetic_painter.reference_fasta` | `reference.fasta` (Wuhan/Hu-1/2019) | Public | **Committed.** The id matches ncov's `refine.root`, so an emitted copy is self-rooting — 30 KB |
| `epihiper.output_csv` | Simulation replicates (~550 MB each) | Your pending deposit | Pending |
| `ascertainment.parameters` | `ascertainment_parameters.yaml` | TwinSampler repo | Ships with that package |
| `ascertainment.schedule_input` | `Virginia_importation_schedule.csv` | TwinSampler `Data/` | Ships with that package |

### Sequence redistribution is not the blocker (resolved 2026-10-02)

Earlier revisions of this document treated GISAID terms as the main obstacle
to external adoption. That framing was wrong on both counts:

* **PhyloGAS stores no sequences.** Both sequence inputs are fetched on
  demand, never committed.
* **Cov-Spectrum serves open sequences by default**, and that is what
  `seq_prep.py` queries in both of its modes — `--seed_mode` for importation
  seeds and bulk mode for the training alignment. So neither
  `genetic_painter.seed_fasta` nor `genetic_painter.align_fasta` carries a
  redistribution constraint; they are simply large and reproducible, which is
  why `.gitignore` excludes them.

What remains worth doing, as a convenience rather than a workaround:

1. **Ship a pre-trained mutation model.** The entropy map is small — 318 KB of
   thresholds plus a 3.6 MB matrix — and training it is the slowest part of a
   first run. A copy in `data/example_data/` that the config can point at
   would let a new user go straight to `phylogas paint`. Not yet wired in.
2. **Document the retraining recipe** for users bringing their own MSA. Mostly
   covered by `fetch-data --with-training-sequences` plus the `training`
   config block.

So a new user needs the EpiHiper replicate and the demographics table; both
sequence inputs come down from Cov-Spectrum on demand.

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


---

## EpiHiper simulation replicates (Zenodo, added 2026-09-30)

The transmission networks `phylogas paint` consumes are now published:

**"Phylogeographic Analysis Similars - Agent Based Simulations"**
`doi:10.5281/zenodo.23067670` (concept `10.5281/zenodo.23067669`), open access.

| state | file | size | md5 (first 12) |
|---|---|---|---|
| va | `va_replicate_0_output.csv.gz` | 530.3 MB | `2971f450fde9` |
| ga | `ga_replicate_0_output.csv.gz` | 583.7 MB | `508eab5b6a04` |
| ma | `ma_replicate_0_output.csv.gz` | 431.8 MB | `bfdcc576f006` |
| mn | `mn_replicate_0_output.csv.gz` | 428.1 MB | `ac6f078664fe` |
| wa | `wa_replicate_0_output.csv.gz` | 425.2 MB | `9be12c1f8726` |

No California replicate is published; `fetch-data` reports this rather than
failing or silently skipping.

```bash
phylogas fetch-data --states va --with-simulations
```

Opt-in, because the simulations roughly quintuple the per-state download
(0.10 GB -> 0.62 GB for VA).

### Checksum verification

All Zenodo downloads are MD5-verified against the record metadata. Behaviour:

* file absent -> download, verify, keep
* file present and matching -> skipped without re-downloading
* file present and **not** matching -> reported, deleted, re-downloaded

That last case matters: a truncated earlier attempt used to be silently
reused. Verified by truncating a 425 MB file to 400 MB and re-running:

```
checksum MISMATCH on existing file (got e83fc0853999..., expected 9be12c1f8726...); re-downloading
checksum OK
```

`--no-verify` skips it if you need the speed.

### Retry / backoff

Both downloaders retry up to 4 times with exponential backoff (5s, 10s, 20s)
on 5xx and network errors. Added after a live Dataverse `503` aborted a
multi-state fetch mid-run.

### Verified end to end

Fetched `wa` (4 population files + the 425 MB replicate), MD5 verified, then:

```
phylogas paint --input-graph-csv wa/wa_replicate_0_output.csv.gz \
    --painted-prefix E2 --start-tick 124 --num-ticks 100 ...
  Wrote 390,885 records in 0.9 min (6,865 rec/s)
  Dropped 399 importations that had no seed sequence (3,322 supplied for 3,721)

phylogas benchmark sequence --painted out/wa_test.sequences.fasta.xz ...
  mean_substitutions    0.204985
  max_substitutions     4
  pairs_over_8          0
```

Two notes from that run:

* **The replicates are two-variant.** WA carries both `E1` (3.37M rows, ticks
  40-359) and `E2` (8.82M rows, ticks 124-699). `--painted-prefix` selects
  which wave to paint; the VA example data behaves the same way.
* **Seed shortfall is real for non-VA states.** WA has 3,721 `E2` importations
  but the bundled seed FASTA is Virginia Delta with 3,322 sequences. The
  shortfall warning fires correctly and names the affected tick range. Per-state
  seed sets (via `phylogas prep-seeds`) are needed for full coverage.


---

## Switched to Zenodo (2026-10-01)

All six synthetic populations are now mirrored on Zenodo and `phylogas
fetch-data` downloads from there. Dataverse is no longer contacted.

| state | record | core files |
|---|---|---|
| va | [23088534](https://doi.org/10.5281/zenodo.23088534) | 104.0 MB |
| ca | [23088872](https://doi.org/10.5281/zenodo.23088872) | 455.4 MB |
| ga | [23088629](https://doi.org/10.5281/zenodo.23088629) | 123.3 MB |
| ma | [23088723](https://doi.org/10.5281/zenodo.23088723) | 79.1 MB |
| mn | [23088773](https://doi.org/10.5281/zenodo.23088773) | 74.9 MB |
| wa | [23088830](https://doi.org/10.5281/zenodo.23088830) | 91.0 MB |

Each record holds all 13 Dataverse files; `fetch-data` pulls the four needed
to build the demographics table, and `--with-epihiper-inputs` adds the contact
network and persontrait database.

**Why:** measured 2026-09-30, Dataverse returned 502 on its website, 500 on
its metadata API, and roughly 50% 503 on its file API. Zenodo has been
reliable for the ABM replicates throughout.

**Checksums are now real.** The filenames, sizes and MD5s in
`src/phylogas/sources.py` were harvested from the Zenodo API rather than
transcribed, so a corrupt or truncated download is always caught.

Verified end to end on Minnesota: four files downloaded and MD5-verified, then
joined into 5,218,000 demographics rows with 100% coordinate coverage and
coordinates inside the state bounding box. A re-run verifies the existing
files and downloads nothing.

The Dataverse DOIs remain in `sources.py` as `DATAVERSE_DOIS` for citation and
manual fallback.


---

## Benchmark inputs resolve automatically (2026-10-01)

`benchmark.allevents` and `benchmark.truth_mugration` default to `"auto"`.
TwinSampler names its outputs deterministically, so the derivation is exact
rather than a guess:

```
ascertainment.output = results/02_simulated_linelists/linelist.csv
  -> results/02_simulated_linelists/linelist_allevents.csv.xz
  -> results/02_simulated_linelists/linelist_mugration.json
```

Each key accepts `"auto"`, `"none"`, or an explicit path. When a file is
absent the dependent benchmark is skipped with a stated reason and the
Snakemake rule is not created, so nothing fails and nothing is silently
scored.

`phylogas validate-config` reports what resolved, including size and mtime so
a stale file from an earlier run is visible:

```
Benchmarking inputs:
  [okay] all-events (ABM truth)    results/.../linelist_allevents.csv.xz
         1.9 MB, modified 0h ago  [auto from ascertainment.output]
  [ -- ] mugration truth           skipped -- auto: not found at results/.../linelist_mugration.json
```

Note that TwinSampler writes the mugration file only when the events carry a
`county` column, so `auto` resolving to nothing there can be legitimate.

`benchmark.infections` was renamed `benchmark.allevents`. It is the *input* to
`phylogas assign-variants`, not the variant-annotated output; the old name
invited pointing it at that rule's own output. The previous key is still read
as a fallback.

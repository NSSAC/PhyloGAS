# PhyloGAS
## Sequence Acquisition / Preparation

### `seq_prep.py`

This script is a flexible tool for acquiring FASTA sequences from the Cov-Spectrum API. It operates in two main modes, controlled by the `--seed_mode` flag.

#### Mode 1: Seed Finding Mode (`--seed_mode`)

Reads a cluster table (from the UCSC SARS-CoV-2 Genome Browser) and writes two
things that must agree with each other:

* **the importation schedule** -- how many importations land on each day, one
  per cluster, dated by that cluster's first sample
* **the seed FASTA** -- the real sequence that founds each of those
  importations, written in schedule order

They are derived from the same rows, in `cluster_seeds.py`, because they used
to be produced separately and disagreed. The schedules EpiHiper ran came from
`src/importation_analysis.py`; seeds came from here and took
`samples.split(",")[0]`, but that list is ordered **newest first**. For 79% of
Virginia's Delta clusters the chosen sequence was the cluster's *last* sample,
up to 209 days after the importation it founded. The painter assigns seed
records positionally as importations come up, so early importations received
late genomes (rank correlation between the two dates: -0.54), inflating
root-to-tip divergence early in the study.

Now the seed of a cluster is its **earliest** sample, and every seed's
collection date equals its importation date exactly.

**Retrieval and substitution.** Sequences are requested in batches (one POST
per `--batch_size` strains), then the response is checked to see which strains
came back. An importation whose sequence is missing falls back to its cluster's
next oldest sample and is retried, up to `--max_fetch_rounds` times. A cluster
that exhausts every sample is dropped and reported -- and because the schedule
is written from the surviving clusters, its importation leaves the schedule
too. That matters: a seed FASTA with one fewer record than the schedule expects
would shift every later pairing by one.

**Surrogates.** A cluster with no sample in Cov-Spectrum's open data at all --
397 of Washington's 3,739 Delta clusters, against about none for Virginia --
borrows the retrieved sequence of the nearest cluster of the same sublineage,
and the manifest's `surrogate_of` column names the donor. Without that the seed
set is short of the importations the ABM made, so the painter has no genome for
those chains and skips every descendant. One donor may serve several
importations. `--no_surrogates` drops them instead.

**Ticks.** The schedule's own `tick` column counts from the first importation.
With `--abm_config` (or `--tick_zero`) the schedule is also written on the
ABM's absolute ticks -- `<State>_<lineage>_ticks.csv`, the `tick,count` form
EpiHiper's seeding template consumes. Reading tick 0 from the ABM's config
(EpiHiper's `config.json`, key `tickZero`) keeps one source of truth; the
bundled experiment uses 2020-11-30 for every state, so each state's ticks sit
on a shared calendar.

**Outlier filtering.** `--outlier_method` now also removes the cluster's
importation from the schedule, not just its seed. The published EpiHiper
schedules were built with **no** outlier filtering, so use `none` to reproduce
them.

**Usage:**
```bash
python seq_prep.py \
  --state "Virginia" \
  --pango "B.1.1.7,B.1.617.2" \
  --output_folder ../../../data/importations/sequences \
  --seed_mode \
  --outlier_method none \
  --abm_config ../../../cfg/exp1/config.json
```

Writes, per lineage, `Virginia_<lineage>_seed_strains.txt`,
`Virginia_<lineage>_seed_sequences.fasta` (in importation order),
`Virginia_<lineage>_seed_manifest.csv` and `Virginia_<lineage>_ticks.csv`,
plus the consolidated `Virginia_schedule.csv`.

**The manifest** lists, for each FASTA record in order, the strain, its
cluster, its importation date and (with a tick 0) its ABM tick. `phylogas
paint` passes it to the painter automatically when it sits beside the seed
FASTA, and the painter then checks two things before painting: that record
*i* is the strain listed for importation *i*, and that importations per day in
the transmission log match seeds per day. Either mismatch is a loud warning
(`--strict_seed_pairing` makes it fatal). Run against the old Virginia seed
FASTA, the order check fires on the first record.

Verified against what the simulations actually ran: with `--outlier_method
none`, the generated schedules reproduce all five
`data/importations/schedules/<State>_schedule.csv` files (every row and
importation count) and all ten `data/seed_schedule/<st>-<n>.csv` files
byte-for-byte.

One correction comes with that: the committed `<State>_schedule.csv` files have
a misaligned `seq_count` column. It was assigned from a `groupby("tick")` into
rows grouped by `(tick, variant)`, so from the first day carrying both variants
onward the values belong to the wrong row (and the last 44 rows are empty). The
`importations` column -- the one EpiHiper consumes -- is correct, and
regenerating fixes `seq_count`.

#### Mode 2: Bulk Download Mode (Default)

This is the default mode. It bypasses the cluster analysis and directly queries the Cov-Spectrum API for all available sequences that match the specified metadata (state, Pango lineage, and an optional date range).

**Usage:**
```bash
python seq_prep.py \
  --state "California" \
  --pango "JN.1" \
  --output_folder ./bulk_sequences/ \
  --date_from "2023-11-01" \
  --date_to "2024-02-29"
```
This command will download all JN.1 sequences from California between the specified dates and save them to `./bulk_sequences/California_JN_1_2023-11-01_2024-02-29.fasta`. If the date parameters are omitted, it will download all sequences for all time.


### `run_all_states.sh`

This script automates the bulk download of sequences using `seq_prep.py`. It iterates through a predefined list of all US states and territories and runs a download job for each one.

**Usage:**
1.  **Configure:** Open the script and set the `PANGO_LINEAGE` variable. You can also optionally set a `DATE_FROM` and `DATE_TO`.
2.  **Make Executable:** `chmod +x run_all_states.sh`
3.  **Run:** `./run_all_states.sh`

The script will create a directory (e.g., `jn1_sequences_by_state`) and populate it with FASTA files, one for each state.



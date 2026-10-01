# Command map: old invocations to the `phylogas` CLI

Every old command still works. The package restructure moved files but kept the
underlying scripts runnable, and the mutational-model imports fall back to flat
imports when a script is executed directly. This table is the migration guide,
not a deprecation notice.

Last verified: 2026-09-28.

---

## Quick reference

| Stage | Old | New |
|---|---|---|
| 0. Seeds | `python seq_prep.py --seed_mode ...` | `phylogas prep-seeds --config config.yaml` |
| 1. Train | `python genetic_painter.py --analysis_type entropy_analysis ...` | `phylogas train --config config.yaml` |
| 2. Paint | `python genetic_painter.py --analysis_type generate_sequence ...` | `phylogas paint --config config.yaml` |
| 3. Linelist | `python simulate_linelist.py ...` | `simulate_linelist ...` (TwinSampler) |
| 4. Sampling | `python3 run_all_scenarios.py ...` | `scenarios-runner ...` (BeyondBaseline) |
| 5. Subset | `python subset_fasta_streaming.py -m ... -f ... -o ...` | `phylogas subset-fasta -m ... -f ... -o ...` |
| 6. Nextstrain | `snakemake --snakefile .../ncov/Snakefile ...` | rule `nextstrain_build` (set `nextstrain.enabled: true`) |
| all | *(run each by hand)* | `phylogas run --config config.yaml --cores all` |

---

## File relocations

| Old path | New path |
|---|---|
| `src/genetic_painter/genetic_painter.py` | `src/phylogas/painter/genetic_painter.py` |
| `src/genetic_painter/mutational_models/` | `src/phylogas/painter/mutational_models/` |
| `src/genetic_painter/cosine_entro.py` | `src/phylogas/painter/cosine_entro.py` |
| `src/genetic_painter/sampling.py` | `src/phylogas/painter/sampling.py` |
| `src/genetic_painter/test/` | `src/phylogas/painter/test/` |
| `src/seq_prep/` | `src/phylogas/seqprep/` |
| `src/filter_variant_time.py` | `src/phylogas/filter_variant_time.py` |

Moved with `git mv`, so `git log --follow` still works.

Earlier relocation, for context: the painter came from
`synthetic_biosurveillance/src/genetic_painter.py`. That copy is now stale
(984 lines vs 1227) and should be deleted from that repo — see
`docs/salvage_audit.md`.

---

## Stage 0: seed acquisition

**Old**
```bash
python seq_prep.py \
    --state Virginia \
    --pango B.1.617.2 \
    --output_folder ../data/importations/sequences \
    --seed_mode \
    --outlier_method chaining
```

**New**
```bash
phylogas prep-seeds --config config.yaml --seed-mode
```

Config keys: `population.state_name`, `variant.pango`, `seeds.output_folder`,
`seeds.outlier_method`.

---

## Stage 1: entropy training

**Old** (`run.03.vadelta.a`)
```bash
python ../genetic_painter.py \
    --analysis_type entropy_analysis \
    --random_number_seed 43 \
    --threshold_file  run.03.vadelta.output/run.03.threshold.file \
    --align_fasta     /project/.../clean_va_delta_sequences.fasta \
    --base_threshold_df run.03.vadelta.output/run.03.base.threshold.df
```

**New**
```bash
phylogas train --config config.yaml
```

Config keys: `genetic_painter.align_fasta`, `.entropy_thresholds`,
`.probability_matrix`, `random_seed`.

Note: `np.save` appends `.npy`, so `base_threshold_df` written as
`...df` becomes `...df.npy` on disk. The config uses the `.npy` form for both
read and write; the loader tolerates either.

---

## Stage 2: painting

**Old** (the full production invocation)
```bash
python ../src/genetic_painter/genetic_painter.py \
    --analysis_type          both \
    --random_number_seed     43 \
    --threshold_file         "run.03.vadelta.output/run.03.threshold.file" \
    --base_threshold_df      "run.03.vadelta.output/run.03.base.threshold.df.npy" \
    --start_date             2021-04-07 \
    --input_graph_csv        ".../replicate_0/output.csv.gz" \
    --seed_fasta             ".../Virginia_B_1_617_2_seed_sequences.fasta" \
    --align_fasta            ".../clean_va_delta_sequences.fasta" \
    --output_prefix          ".../run_03_vadelta_2026_09_10_128to428" \
    --compression            xz \
    --persontrait_file       ".../va_2_4_0_demographics.csv" \
    --add_metadata           "gender,county,home_latitude,..." \
    --input_graph_painted_prefix "E2" \
    --proportional \
    --reference              ".../reference.fasta" \
    --location               '{"country":"USA",...}' \
    --reference_location     '{"country":"China",...}' \
    --rate_limit \
    --initial_viral_load     10 \
    --start_tick             128 \
    --num_ticks              300
```

**New**
```bash
phylogas paint --config config.yaml
```

Config mapping:

| Old flag | Config key |
|---|---|
| `--threshold_file` | `genetic_painter.entropy_thresholds` |
| `--base_threshold_df` | `genetic_painter.probability_matrix` |
| `--input_graph_csv` | `epihiper.output_csv` |
| `--seed_fasta` | `genetic_painter.seed_fasta` |
| `--align_fasta` | `genetic_painter.align_fasta` |
| `--output_prefix` | `genetic_painter.output_prefix` |
| `--reference` | `genetic_painter.reference_fasta` |
| `--persontrait_file` | `population.persontrait_file` |
| `--add_metadata` | `genetic_painter.add_metadata` |
| `--input_graph_painted_prefix` | `genetic_painter.painted_prefix` |
| `--rate_limit` | `genetic_painter.mutation_model: rate_limit` |
| `--poor` | `genetic_painter.mutation_model: poor` |
| *(neither)* | `genetic_painter.mutation_model: simple` |
| `--proportional` / `--neutral` | `genetic_painter.proportional: true/false` |
| `--initial_viral_load` | `genetic_painter.initial_viral_load` |
| `--start_date` / `--start_tick` / `--num_ticks` | `genetic_painter.start_date` / `.start_tick` / `.num_ticks` |
| `--compression` | `genetic_painter.compression` |
| `--random_number_seed` | `random_seed` |
| `--location` | `genetic_painter.location` |

Any flag may still be passed explicitly to override the config:
```bash
phylogas paint --config config.yaml --num-ticks 30 --compression-level 0
```

Use `--dry-run` to print the underlying `genetic_painter.py` command without
running it — useful when checking that a config maps to the flags you expect.

### New flags with no old equivalent

| Flag | Config key | Why |
|---|---|---|
| `--compression-level` | `genetic_painter.compression_level` | lzma preset. Default is now **1**, not 6. Level 6 runs at ~243 rec/s vs ~9,000 at level 1 — roughly 6 extra hours on a 5.3M-record run for <5% size gain. |
| `--compression-threads` | `genetic_painter.compression_threads` | Routes output through external `xz -T<n>` instead of single-threaded in-process lzma. |

---

## Stage 3: ascertainment (TwinSampler)

Unchanged apart from packaging. Once `pip install -e .` has been run in the
TwinSampler repo, `simulate_linelist` is on `PATH`:

```bash
simulate_linelist \
  --epihiper .../output.csv.gz \
  --people   .../va_persontrait_epihiper.txt \
  --households .../va_household.csv \
  --rucc     .../Ruralurbancontinuumcodes2023.csv \
  --ascertain .../ascertainment_parameters.yaml \
  --start_date 2021-04-07 --start_tick 128 --stop_tick 428 \
  --out .../linelist.csv --seed 42 --output_all_events \
  --schedule_input .../Virginia_importation_schedule.csv --variant_mode 2
```

Driven by the `simulate_linelist` rule in the Snakefile.

---

## Stage 4: adaptive sampling (BeyondBaseline)

**Old**
```bash
python3 run_all_scenarios.py --linelist ... --population ... --algorithms surs
```

**New** (after `pip install -e .` in BeyondBaseline)
```bash
scenarios-runner --linelist ... --population ... --algorithms surs
```

Full entry-point list:

| Command | Purpose |
|---|---|
| `scenarios-runner` | Full scenario sweep (`run_all_scenarios.py`) |
| `scenarios-weekly` | Operational weekly loop — the "lite" mode (`run_weekly_sampling.py`) |
| `scenarios-replicates` | Drive multiple replicates |
| `scenarios-lasso-greedy` | LASSO greedy KL-minimisation sweep |
| `scenarios-lasso-stratified` | LASSO stratified sweep |
| `scenarios-aggregate` | Aggregate metrics across replicates |

### The operational "lite" mode

`vision.txt` specifies a health-department-facing command:

```bash
beyond_baseline.py -l linelist.txt -g already_sequenced.txt -p population.txt \
    -b <budget> --time_budget 4 --current_date <t> --algorithm SURS > recommendation.txt
```

The closest existing implementation is `scenarios-weekly`:

```bash
scenarios-weekly \
  --linelist linelist.csv \
  --population population.txt \
  --already-sequenced weekly_results/history_combined_all.csv \
  --target "LL,P" \
  --current-date 2021-07-12 \
  --batch-size 200 \
  --no-replacement \
  --algorithms surs
```

Argument names differ from the vision (`--linelist` vs `-l`,
`--batch-size` vs `-b`, `--already-sequenced` vs `-g`). Adding short aliases and
a `beyond-baseline` console alias would make the vision's line literally
runnable; see `UNCOMMITTED_CHANGES.md` in that repo.

---

## Stage 5: FASTA subsetting

**Old**
```bash
python ./seq_prep/subset_fasta_streaming.py \
  -m .../linelist.csv__seed42_scenario1_SURS_samples.csv.xz \
  -f .../run_03_vadelta.sequences.fasta.xz \
  -o ../outputs/run_03_vadelta.SURS.fasta.xz
```

**New**
```bash
phylogas subset-fasta -m <samples.csv.xz> -f <full.fasta.xz> -o <subset.fasta.xz>
```

The streaming implementation is used (single pass, no index), so it is safe on
multi-GB FASTA. The `pyfastx`-indexed variant remains at
`src/phylogas/seqprep/subset_fasta.py` for random-access use.

---

## Stage 6: Nextstrain

**Old**
```bash
snakemake --snakefile .../ncov/Snakefile \
          --configfile .../run_03_vadelta.SURS.yaml \
          --cores all --rerun-incomplete
```

**New** — set `nextstrain.enabled: true` and point `nextstrain.snakefile` /
`.configfile` at your ncov checkout, then it runs as part of `phylogas run`.
Kept opt-in because ncov is a large external dependency.

---

## SLURM submission scripts

`src/phylogas/painter/test/` still holds the hand-written scripts:

| Script | Role |
|---|---|
| `asw_run_latest.sh` | sbatch wrapper |
| `run.03.vadelta.a` | entropy training only |
| `run.03.vadelta.b` | painting only |
| `run.03.vadelta.both` | both stages |

They now anchor to their own location (`BASH_SOURCE`), so submission directory
no longer matters. They remain useful for one-off runs; `phylogas run --profile
slurm` is the path for sweeps across many replicates.

---

## Not yet implemented

| Promised | Status |
|---|---|
| `phylogas fetch-data` | **Implemented.** Downloads the four core population files per state from the six per-state Dataverse deposits, then builds the demographics table. Contact networks are opt-in (`--with-epihiper-inputs`). |
| `phylogas build-demographics` | **Implemented.** Joins persontrait + person + household + residence_locations + FIPS into the table the painter reads. |
| `phylogas benchmark` | **Implemented.** `mugration` / `sequence` / `compare`. Moved from BeyondBaseline; verified bit-identical. |
| Docker/Apptainer image | `Dockerfile` not yet written. |
| `nextstrain_build` rule | Written but untested — needs a local ncov checkout. |


---

## Data acquisition

**Old** (manual, or `download_pgcoe_dataverse-v2.sh`)
```bash
./download_pgcoe_dataverse-v2.sh va ga
# then, separately and by hand, build the demographics table
```

**New**
```bash
phylogas fetch-data --states va ga            # download + build demographics
phylogas fetch-data --states va --dry-run     # show the plan and sizes first
```

Key differences from the shell script:

| | script | `fetch-data` |
|---|---|---|
| contact network (7.9 GB / 6 states) | always | opt-in `--with-epihiper-inputs` |
| `residence_locations` | **missing** | included (required for coordinates) |
| decompression | always, `-k` keeps both (~40 GB) | off by default; pandas reads `.xz` |
| demographics build | not performed | automatic |
| default state set | all six | from config, or `--states` |

To build the demographics table on its own:
```bash
phylogas build-demographics \
    --persontrait va/va_persontrait_epihiper.txt.xz \
    --person      va/va_person.csv.xz \
    --household   va/va_household.csv.xz \
    --residence   va/va_residence_locations.csv.xz \
    --fips        data/county_fips.csv \
    --out         data/va/va_2_4_0_demographics.csv
```

### Config key rename

`population.persontrait_file` -> `population.demographics_file`, because the
file is a derived join and never the raw EpiHiper persontrait. The old key is
still honoured as a fallback.


---

## Benchmarking (moved from BeyondBaseline, 2026-09-30)

**Old** — one command, but required ABM truth even for demographic metrics:
```bash
python3 run_all_scenarios.py --linelist ... --infections ... \
    --abm_mugration linelist_mugration.json --save-samples
```

**New** — selection and scoring separated:
```bash
beyond-baseline-sweep --linelist ... --population ... --save-samples --outdir runs/
phylogas benchmark mugration --truth <abm.json> \
    --samples 'runs/*_samples.csv.xz' --infections <allevents.csv.xz>
```

Or both at once:
```bash
phylogas compare-strategies --truth <abm.json> --infections <allevents.csv.xz> \
    --linelist ... --population ... --algorithms surs LASSO-Greedy
```

| Old | New |
|---|---|
| `--abm_mugration` on the sweep | `phylogas benchmark mugration` |
| `Mugration_Metrics.csv` (BB) | same filename, written by PhyloGAS |
| *(nothing)* | `phylogas benchmark sequence` — parent->child divergence |
| *(nothing)* | `phylogas benchmark compare` — parsimony vs augur |

Why: benchmarking needs ABM ground truth; sampling does not. BeyondBaseline
must stay runnable by a health department on a real linelist. Full migration
guide in `BeyondBaseline/cste_instructions.txt`.


---

## Sublineage scope in seed acquisition

`seq_prep.py` previously hardcoded a trailing `*` on the lineage, so
`--pango B.1.617.2` silently meant "Delta and every descendant". That is
usually what you want, but it was not stated and could not be turned off.

Now controlled by `--include-sublineages` / `--no-include-sublineages`,
**default ON** (preserving the previous behaviour).

```bash
phylogas prep-seeds --config config.yaml --seed-mode                          # B.1.617.2*
phylogas prep-seeds --config config.yaml --seed-mode --no-include-sublineages # B.1.617.2 exactly
```

Config: `variant.include_sublineages: true`

The flag affects both code paths, which implement the scope differently:

| mode | with sublineages | exact only |
|---|---|---|
| bulk (Cov-Spectrum query) | `pangoLineage=B.1.617.2*` | `pangoLineage=B.1.617.2` |
| seed (cluster TSV) | descendants rolled up via `pango_aliasor` | no roll-up |

Measured on Virginia Delta:

| | sublineages ON | exact only |
|---|---|---|
| Cov-Spectrum samples | 21,034 | 836 |
| bulk download, 2021-06-01..15 | 22 seqs | 10 seqs |
| seed strains (cluster TSV) | 3,320 | 171 |

Bulk-mode output filenames gain an `_exact` suffix when sublineages are
excluded, so the two variants cannot overwrite each other.

### Seed acquisition during data fetch

```bash
phylogas fetch-data --states wa --with-simulations --with-seeds
```

Needed for states other than Virginia: the bundled seed FASTA is Virginia
Delta (3,322 sequences), and WA alone has 3,721 E2 importations, so the
painter reports a shortfall.

**Known upstream issue:** `--seed-mode` downloads a cluster TSV from
`clustertracker.gi.ucsc.edu`, whose TLS certificate expired 2025-07-02
(`notAfter=Jul  2 23:59:59 2025 GMT`). Three workarounds, in order of
preference:

```bash
# 1. local copy -- safest; a snapshot ships in data/importations/
phylogas prep-seeds --config config.yaml --seed-mode \
    --input_file data/importations/sarscov2_clusters_2024_11_12_filtered.tsv.gz

# 2. bypass verification -- fetches current data from UCSC
phylogas prep-seeds --config config.yaml --seed-mode --insecure-download

# 3. bulk mode -- queries Cov-Spectrum directly, unaffected by the cert
phylogas prep-seeds --config config.yaml
```

`--insecure-download` disables certificate verification **for that one
request**. The transfer stays encrypted, but the server's identity is not
authenticated, so it prints a warning each time. Secure by default; also
settable via `seeds.insecure_download: true` for unattended runs.

Verified 2026-09-30: with the flag, UCSC serves a current 254 MB /
467,045-row table, yielding 3,323 Virginia Delta seed strains (versus 3,320
from the bundled November-2024 snapshot).


---

## Painting a subset of infections

The painter must walk the whole transmission tree -- a child's genome derives
from its parent's -- but it does not have to *write* all of it.

```bash
# default: every painted infection, as before
phylogas paint --config config.yaml

# also emit one FASTA per sampling strategy, in the same pass
phylogas paint --config config.yaml \
    --linelist-filter results/02_simulated_linelists/linelist.csv \
    --linelist-filter surs=runs/surs_samples.csv.xz \
    --linelist-filter lasso50=runs/LASSO-Greedy_samples.csv.xz
```

```
Output sets:
  all                     29,009 records  out.sequences.fasta.xz
  linelist                 8,412 records  out.linelist.sequences.fasta.xz
  surs                       500 records  out.surs.sequences.fasta.xz
  lasso50                    300 records  out.lasso50.sequences.fasta.xz
```

`--linelist-filter` is repeatable and accepts `LABEL=FILE` or bare `FILE`
(label taken from the filename). Identifiers are read from whichever of
`alias_pid`, `infection_id`, `strain`, `sim_pid` or `pid` the file provides,
so TwinSampler linelists, BeyondBaseline samples files and painter metadata
all work unchanged.

Verified: filtered sets are strict subsets of the unfiltered run with
byte-identical sequences, and an unfiltered run is byte-identical to the
previous behaviour.

### Post-hoc alternative

```bash
phylogas subset-fasta -m runs/surs_samples.csv.xz \
    -f all_infections.fasta.xz -o surs.fasta.xz
```

Produces exactly the same output as the in-pass filter (verified). Use
`--linelist-filter` to avoid writing the large FASTA at all; use
`subset-fasta` to re-cut an existing one without repainting.

### `alias_pid` in metadata

Painter metadata now carries `alias_pid` (`{pid}.{tick}`) alongside `strain`,
so joins against TwinSampler and BeyondBaseline outputs are a direct key match
instead of parsing the strain ID.

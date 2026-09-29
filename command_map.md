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
| `phylogas benchmark` | Missing. Topological F1 / mugration cosine exist in `BeyondBaseline/scripts/scenarios_simulation/mugration_station.py` and should be promoted rather than rewritten. |
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

# Command map: old invocations to the `phylogas` CLI

Every old command still works. The package restructure moved files but kept the
underlying scripts runnable, and the mutational-model imports fall back to flat
imports when a script is executed directly. This table is the migration guide,
not a deprecation notice.

Last verified: 2026-10-02.

---

## Quick reference

| Stage | Old | New |
|---|---|---|
| 0. Seeds | `python seq_prep.py --seed_mode ...` | `phylogas prep-seeds --config config.yaml --seed-mode` |
| 0b. Training seqs | `python seq_prep.py ...` (bulk mode, by hand) | `phylogas fetch-data --with-training-sequences` |
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

---

## Panel letters retired (2026-10-02)

The scenario sweeps labelled every metric with a single letter -- "panel A",
"panel F" -- and the letter went into the `eval_type` column and into output
filenames, so a results directory held
`lasso_all_scenarios_K_coverage_size_100.csv` with nothing to say what K was.
The letters are gone. The new name is the old one with the letter stripped,
except `A_targets`, which became `kl_targets`.

| Old `eval_type` | New `eval_type` | Computed by |
|---|---|---|
| `A_targets` | `kl_targets` | BeyondBaseline |
| `B_cumulative_infections` | `cumulative_infections` | PhyloGAS (needs ABM truth) |
| `C_stride_window_infections` | `stride_window_infections` | PhyloGAS (needs ABM truth) |
| `E_stride_variant_prevalence_error` | `stride_variant_prevalence_error` | PhyloGAS (needs ABM truth) |
| `F_stride_component_coverage` | `stride_component_coverage` | PhyloGAS (needs ABM truth) |
| `I_coverage_size_0` | `coverage_size_0` | PhyloGAS (needs ABM truth) |
| `J_coverage_size_10` | `coverage_size_10` | PhyloGAS (needs ABM truth) |
| `K_coverage_size_100` | `coverage_size_100` | PhyloGAS (needs ABM truth) |
| `L_coverage_size_1000` | `coverage_size_1000` | PhyloGAS (needs ABM truth) |
| `M_8_week_rolling_tree_coverage` | `8_week_rolling_tree_coverage` | PhyloGAS (needs ABM truth) |
| `N_equity_<stratifier>` | `equity_<stratifier>` | BeyondBaseline |

`D_`, `G_` and `H_` appeared in the notebooks' prefix tuples but were never
emitted by any code in either repository. They are dropped, not renamed.

Output filenames followed the rename:

| Old | New |
|---|---|
| `lasso_all_scenarios_B_cumulative_infections.csv` | `lasso_all_scenarios_cumulative_infections.csv` |
| `lasso_all_scenarios_K_coverage_size_100.csv` | `lasso_all_scenarios_coverage_size_100.csv` |
| `A_table3_targets_1xN.png` | `kl_targets_table_1xN.png` |
| `N_equity_heatmap_<age>.png` | `equity_heatmap_<age>.png` |
| `Ranking_Selected_Panels_ACFHKM.csv` | `Ranking_Selected_Metrics.csv` |

### Where the definitions live

`BeyondBaseline/scripts/scenarios_simulation/eval_names.py` is the registry of
record. It declares, per metric, whether a higher value is better, which
metric family it belongs to, and whether it needs ABM ground truth. Two things
that used to be inferred from the letter prefix -- direction of improvement
and metric family -- were spelled differently in each notebook
(`('G_','H_','I_','J_','K_','L_','M_')` in one, `('F_','I_',...)` in another);
they now come from one place.

`canonicalize()` maps any old letter-coded name onto the new one, and the
aggregator and the notebooks apply it when reading result CSVs, so existing
`AUC_rankings.csv` / `Aggregated_Median_AUC_Rankings.csv` files still load.

### Why the split is by ground truth, not by "sweep vs evaluation"

BeyondBaseline keeps `kl_targets` and `equity_*` because they need only the
line list it was handed -- a health department can compute them on real data.
Everything that needs the ABM's true infection counts or its transmission
graph is a PhyloGAS metric. The sweep itself needs no ground truth at all:
`--infections` is only ever read by the truth metrics.

The LASSO subgroup sweeps had broken on this, because they imported the moved
functions at module scope. They now import them optionally:

```python
try:
    from phylogas.benchmark.truth_metrics import (...)
    _HAVE_TRUTH = True
except ImportError:
    _HAVE_TRUTH = False
```

so `beyond-baseline-lasso-greedy` runs against a line list with no PhyloGAS
installed, printing which metrics it skipped.

---

## Status and input resolution (2026-10-02)

`phylogas status` changed in three ways after a cluster run sent the user in
circles:

1. **Training sequences are a tracked stage.** `genetic_painter.align_fasta`
   was never checked, so nothing told you to run bulk mode before `train`.
2. **Every missing input is reported**, not just the first. The old report
   named one next step -- `train` -- while the seed FASTA was also absent, so
   the rest had to be found by backtracking.
3. **Seeds are resolved, not just looked up.** `prep-seeds` writes into
   `seeds.output_folder` and `fetch-data --with-seeds` writes into
   `<data_dir>/<state>/seeds`, while the painter reads
   `genetic_painter.seed_fasta`. All three readers (`status`,
   `validate-config`, `paint`) now check the configured key and then those two
   locations, reporting where the file was found:

```
  [okay] Seed sequences    0.0 MB  data/importations/sequences/Virginia_B_1_617_2_seed_sequences.fasta
         found via seeds.output_folder; config points elsewhere
```

   so a fresh `prep-seeds` run no longer needs the YAML hand-edited.

`fetch-data` also stopped announcing `Ruralurbancontinuumcodes2023.csv` as
"still required from elsewhere" unconditionally. It now looks at
`population.rucc_file`, the download destination, `./data/` and the copy
committed in the repository, and says nothing when any of them has it.

---

## Training sequences for the mutational model

`seq_prep.py` has two independent modes. `--seed_mode` picks importation seeds
from a UCSC cluster TSV; without it, **bulk mode** queries Cov-Spectrum over a
date range and writes one FASTA. Bulk mode is what produces the training
alignment the entropy model is fitted to, and it was missing from the CLI: the
date-range flags existed in the script but `prep-seeds` never passed them, so
a bulk run would have queried the lineage's whole history.

```bash
phylogas fetch-data --states va --with-training-sequences
phylogas fetch-data --states va --with-training-sequences --dry-run   # show the window first
```

`--training-sequences` is accepted as an alias. The window is:

| config | window used |
|---|---|
| `training.date_from` + `training.date_to` set | exactly those |
| unset | `genetic_painter.start_date` .. `start_date + num_ticks` |

The derived default is the simulated window, matching what the existing runs
did. The command always prints the window and where it came from, so this is
never silent:

```
Training sequences (Cov-Spectrum bulk mode):
  target: data/clean_va_delta_sequences.fasta
  window: 2021-04-07 .. 2022-02-01   (genetic_painter.start_date + num_ticks; set training.date_from/date_to to widen)
```

A mutational model is usually better fitted to a window that *starts earlier*
than the wave being simulated, so set `training.date_from` explicitly when
that matters -- the default is a convenience, not a recommendation.

Bulk mode names its own output (`<State>_<lineage>_<from>_<to>.fasta`), so the
command moves it onto the configured `genetic_painter.align_fasta` path. Leave
that key unset and it derives
`<data_dir>/training_sequences/<State>_<lineage>_<from>_<to>.fasta` and tells
you what to set. An existing file is never re-downloaded.

---

## Output layout: results scoped by project (2026-10-02)

Two rules now decide where a file goes:

* `data_dir` holds **inputs** — anything fetched or handed to the pipeline.
* `results_dir` holds everything **derived**, under numbered stages.

`results_dir` carries the project name, so parallel projects cannot write over
each other and one run is a single directory to archive or delete:

```yaml
results_dir: "results/{project_name}"
```

Every stage inherits that without further change:

```
results/va_delta_wave/00_mutation_model/     <- new: the entropy map
results/va_delta_wave/01_synthetic_genomes/
results/va_delta_wave/02_simulated_linelists/
results/va_delta_wave/03_sampled_datasets/
results/va_delta_wave/04_nextstrain_builds/
```

### What moved

| Config key | Old default | New default |
|---|---|---|
| `genetic_painter.entropy_thresholds` | `{data_dir}/example_data/run.03.threshold.file` | `{results_dir}/00_mutation_model/{project_name}.thresholds.txt` |
| `genetic_painter.probability_matrix` | `{data_dir}/example_data/run.03.base.threshold.df.npy` | `{results_dir}/00_mutation_model/{project_name}.base_threshold_df.npy` |
| `genetic_painter.align_fasta` | `{data_dir}/clean_va_delta_sequences.fasta` | `{data_dir}/training_sequences/va_delta.fasta` |
| `genetic_painter.seed_fasta` | `{data_dir}/example_data/Virginia_..._seed_sequences.fasta.gz` | `{seeds.output_folder}/Virginia_..._seed_sequences.fasta` |

The entropy map was the clear error: it is the *output* of `train`, but the
default wrote it into `data/example_data/`, which is gitignored scratch and
named as if it were sample input. The `run.03.*` names were a label from the
cluster test runs, inherited by every new project.

Existing configs keep working — these are template defaults, and a config
pointing at a staged or pre-trained map is resolved as before.

### Directories are created now

The painter does not create its own output directories, and nothing upstream
did either, so `train` failed outright when the configured folder did not
exist yet. `train` and `paint` now create what they are about to write into
(skipped under `--dry-run`).

### Also fixed

`status` printed `fetch-data` next-steps without `--config`, so following its
advice produced

```
ERROR: --with-training-sequences needs --config to know the date range and the target path.
```

And `--with-training-sequences` now says plainly that it renamed the file:
`seq_prep.py` prints its own output name, then the command moves it onto
`genetic_painter.align_fasta`, which previously looked like two conflicting
paths scrolling past.

`.gitignore` gained `data/training_sequences/` and `data/clean_*_sequences.fasta`.
The training alignment had been landing in a tracked path, so a ~100 MB FASTA
was one `git add -A` away from being committed.

---

## The strain id as a cross-repo contract (2026-10-03)

Nextstrain joins sequences to metadata on the `strain` string. The painter
writes it as a FASTA header; TwinSampler writes it into the line list's
`strain` column, from a separate implementation. If they disagree by one
character ncov drops the sequence, sometimes without an error.

### What was wrong

`simulate_linelist.format_final_linelist` hardcoded the geography as Python
default arguments, while the painter read the same four values from
`genetic_painter.location`:

```python
def format_final_linelist(raw_linelist_df, country="USA",
                          region="North America",
                          division="Virginia", divisionAbbr="VA"):
```

Both call sites used the bare defaults. So a Washington run would have
produced line list strains reading `USA/VA-EHip-...` against painter headers
reading `USA/WA-EHip-...` -- no sequence matching its metadata. It had not
fired only because every run so far has been Virginia.

Worse, `phylogas paint` never passed `--location` at all, so
`genetic_painter.location` was dead config and the painter always used its own
argparse defaults. `command_map.md` claimed otherwise.

`alias_pid` was also built in three places with three type treatments -- no
cast, `astype(str)`, and `int()` -- which agree only while pid and tick are
integral. A float tick gave `123.45.0` from one and `123.45` from another.
(Two of those live in TwinSampler: `simulate_linelist.py` builds it inside the
strain id, `label_components.py` builds the column itself.)

### What changed

| | before | after |
|---|---|---|
| `phylogas paint --location` | never passed | from `genetic_painter.location` (and `.reference_location`) |
| `simulate_linelist` geography | hardcoded defaults | `--country/--region/--division/--division_abbr` |
| Snakefile | start_date/start_tick shared only | geography shared too, via `LOCATION` |
| `alias_pid` / `strain` | 3 implementations | `phylogas/ids.py` canonical; TwinSampler keeps a pinned copy |

`src/phylogas/ids.py` is the canonical implementation. The painter routes
through it; output is byte-identical for integral input, and non-integral
input now raises instead of producing a second spelling.

TwinSampler deliberately keeps its own copy -- it must run without PhyloGAS
installed -- so both are pinned to one fixture table:

```
PhyloGAS/tests/data/strain_ids.json          <- the table
PhyloGAS/tests/test_strain_id_contract.py
twin_sampler/.../test_data/strain_ids.json   <- same table, copied
twin_sampler/.../test_strain_id_contract.py
```

Change either formula and at least one test fails. Verified: all five fixture
cases agree between the two implementations, across VA/WA/GA geographies, and
a WA run now yields WA strains on both sides.

### Checking a build before running it

```bash
phylogas check-join -f <painted_or_subset.fasta.xz> -m <linelist.csv.xz>
phylogas check-join -f out.surs.sequences.fasta.xz -m runs/surs_samples.csv.xz --column strain
```

Exits non-zero if any FASTA header lacks a metadata row, and reports the
unmatched ids. When nothing matches at all it prints an example from each side
and names the two causes that shift every id at once -- geography and the
date anchor. Spare metadata rows are reported but fine; ncov filters them.

This is a precondition for `nextstrain build`, not a diagnostic, and
`nextstrain-config` will call it before rendering.

### Which metadata feeds Nextstrain

TwinSampler's line list, because it carries `component_id`, `variant_label`,
`rucc_code`, `ses_category` and the rest of the analysis columns. The painter's
metadata has none of those.

The painter's metadata is kept all the same: it is the only record of
`real_strain`, which links a simulated lineage back to the real-world seed
sequence it descended from, and it is the authoritative `strain` <-> `alias_pid`
mapping written by the same code that wrote the FASTA headers. Nothing else
reads it -- not TwinSampler (`--epihiper/--people/--households/--rucc/
--ascertain/--schedule_input`), not `assign-variants` (which reads the
allevents file). In the Snakefile `PAINTED_META` appears only as a paint output
and in the `all` / `paint_only` target lists.

Nothing in TwinSampler consumes its own `strain` column either. The only
consumer anywhere was `filter_fasta_by_metadata` in `utils/sampling.py`, whose
sole caller is its own `__main__`; TwinSampler deleted that file in the
BeyondBaseline split, and `phylogas subset-fasta` supersedes it. Two stale
copies remain, in BeyondBaseline and PhyloGAS -- a salvage-audit item.

So the `strain` column exists for exactly one reason: to be the `strain`
column ncov joins on.

---

## Nextstrain integration (2026-10-03)

**Old** -- a hand-written ncov config per run, with the run name typed into
four places and absolute `/home/anwarren/...` paths in three:

```bash
snakemake --snakefile .../ncov/Snakefile \
          --configfile .../run_03_vadelta_2026_03_22_128to428.SURS.yaml \
          --cores all --rerun-incomplete
```

**New**

```bash
phylogas nextstrain-config --config config.yaml --build-type strategy --algo surs
nextstrain build /path/to/ncov --configfile data/phylogas/surs/config.yaml
```

or as part of `phylogas run` with `nextstrain.enabled: true`.

### Two arms, differing in inputs as well as subsampling

| | `strategy` (default) | `all_infections` |
|---|---|---|
| sequences | the BeyondBaseline subset FASTA | the full painted FASTA |
| metadata | that arm's `_samples.csv.xz` | the all-events line list |
| scheme | `strategy_focal_context` | `country_17k` |
| focal set | taken **whole** | capped at 17,000 |

The default arm never re-subsamples the focal set: BeyondBaseline already chose
it, and capping it again would measure ncov's subsampler instead of the
strategy. The 17k cap is bound to the all-infections input, where the cap *is*
the experiment, so it cannot be applied by accident to a sampled arm.

### base.yaml carries deviations only

ncov loads `defaults/parameters.yaml` and deep merges the user's `--configfile`
over it, so `cfg/nextstrain/ncov/base.yaml` states only what differs:
`coalescent` skyline (ncov: opt) and `clock_filter_iqd` 4 (ncov: 8), plus
`traits`, `files` and the two subsampling schemes.

Four things turned up while merging the a/b configs:

| finding | consequence |
|---|---|
| `filter.group_by` is never read by ncov -- `config["filter"]` is only indexed by input name, `min_length` and `skip_diagnostics` | the long group_by lists in the old configs were inert; grouping belongs in the subsampling scheme |
| `traits.columns` had no `county` in either new config | `phylogas benchmark mugration` reads `models['county']`, so those builds could not have been scored. Now pinned, and `nextstrain-config` refuses a base without it |
| omitting `clock_rate` does **not** make TreeTime infer it | ncov's default 0.0008 applies instead. The refine rule passes `--clock-rate` unconditionally and reads the key with `[]`, so inference is impossible from inside ncov -- hence `phylogas benchmark clock` as a separate step on ncov's own intermediates |
| `reference_id.txt` read `Wuhan-Hu-1/2019` | the reference files and `refine.root` use `Wuhan/Hu-1/2019`; `--include` matched nothing, so the root could be subsampled away. Fixed, and `21L` added |

Also: `USA/WA1/2020` in ncov's defaults is `reference_node_name`, which ncov's
own docs list as **Unused**. The root is still `Wuhan/Hu-1/2019` and the
alignment reference is still `MN908947`, so there was nothing to migrate.

### Why inputs are staged rather than referenced

`nextstrain build <dir>` mounts only that directory. Under the docker and
singularity runtimes an absolute path outside it is invisible, so
`nextstrain-config` hardlinks (or copies) the arm's sequences and metadata into
`<ncov>/data/phylogas/<arm>/` and the shared support files into
`<ncov>/data/phylogas/_shared/`, then writes ncov-relative paths. Referencing
absolute paths would work under `ambient`/`conda` and fail under the runtime
most likely on a cluster.

This is also why the old configs had everything under `.../ncov/data/` -- that
was a container requirement, not an accident.

### Validation is a precondition, not a diagnostic

`nextstrain-config` refuses to write a config when the sequences do not all
join to their metadata, and names the two causes that shift every id at once.
Nothing is staged on failure. `phylogas check-join` runs the same check
standalone.

### Installing it

The CLI and the workflow are separate; installing the CLI does not bring ncov,
which has no `nextstrain-pathogen.yaml` and so is a clone-and-build workflow.
See the README's "Optional: Nextstrain" section. `phylogas status` reports the
checkout and the CLI, distinguishing "not installed" from "installed but no
runtime set up".

---

## Clock rate: three quantities, two modes (2026-10-03)

Full decision record in `docs/clock_modes.md`. The mechanics:

```bash
phylogas benchmark clock --config config.yaml --arm all_infections --mode all
phylogas benchmark clock --config config.yaml --mode truth      # no ncov run needed
phylogas benchmark clock --config config.yaml --arm surs --mode inferred
```

| quantity | from | comparable to real data |
|---|---|---|
| `mu_truth` | root-to-tip vs date, Hamming from the reference, tree-free | yes |
| `mu_sim` / `mu_real` | `augur refine` with no `--clock-rate` | yes |
| `mu_operational` | the rate ncov assumed, read back from `branch_lengths.json` | n/a -- assumed, not measured |

Why root-to-tip and not the substitutions the painter placed: root-to-tip is
what TempEst and TreeTime regress and what the ~8e-4 subs/site/year literature
figure derives from, so it is the only form comparable to real data. Counting
placed substitutions counts *events*, including back-mutations and repeat hits
no real root-to-tip plot can see. Their ratio is a saturation diagnostic.

Validated against synthetic sequences evolved at a known 8.0e-4: recovered
7.94e-4, a 0.7% error, R2 0.76 over a 300-day span.

### Why mu_truth exists at all

`rate_limited` has no clock parameter. It is per replication cycle
(`mutation_rate_per_cycle = 3.40e-6`, Poisson draws per transmission, burst
sizes 10-1000), so the per-site-per-year rate is emergent from the painter's
parameters *and* EpiHiper's transmission timing -- the generation interval does
half the work and lives outside PhyloGAS. The rate has to be measured; there is
no number to look up, and nothing to pin `clock_rate` to.

With mu_truth in hand the comparison decomposes. `mu_sim` vs `mu_truth` asks
whether TreeTime recovers a known clock; `mu_truth` vs `mu_real` asks whether
the simulation's evolution resembles reality. Comparing `mu_sim` to `mu_real`
directly -- the obvious two-way comparison -- conflates them.

### Two things the mode changes besides the prior

`clock_filter_iqd` defaults to **0** in the inferred step. It measures
deviation from the fitted line, so with the rate being fitted it prunes valid
branches, and by differing amounts per arm -- manufacturing a difference
between sampling strategies that is an artefact of filtering.

**The arm matters as much as the mode.** Root-to-tip regression depends on
which tips are present, so a rate fitted on a `strategy` arm is partly a
property of the sampler. Fit on `all_infections` to describe the simulation;
fit per `strategy` arm and the spread is itself a benchmark result -- how much
each sampling design biases the rate estimate, which an empirical study cannot
ask because real surveillance has one realised scheme.

### Also fixed

`cmd_benchmark` imported `benchmark.runner` up front, which pulls in networkx
and the mugration machinery. `benchmark clock` needs neither, so a clock
estimate was impossible in an environment that had augur but not networkx. The
imports are now per branch.

---

## Wiring the declared config into `phylogas run` (2026-10-03)

Four config keys existed but nothing read them, so the commands worked by hand
while `phylogas run` only did the strategy arms with the operational clock.

| key | was | now |
|---|---|---|
| `nextstrain.builds` | declarative; the Snakefile hardcoded `strategy`/`{algo}` | read to derive the arm list |
| `nextstrain.clock.mode` | nothing read it | maps to `benchmark clock` invocations |
| the `all_infections` arm | `nextstrain-config` accepted it; no rule existed | its own rule, with its own inputs |
| `benchmark clock` | not a rule | `rule clock_estimates`, a target of `rule all` |

`nextstrain.builds` expands a `strategy` entry to one arm per algorithm
(inheriting `sampling.algorithms` when the entry names none) and an
`all_infections` entry to a single arm. An empty or absent list falls back to
every sampling arm, so enabling nextstrain alone still does something.

The two config rules are separate rather than one wildcard rule because the
arms differ in their *inputs*, not only their subsampling: strategy arms take
the sampled FASTA and that arm's samples CSV, `all_infections` takes the
painted FASTA and the all-events line list. The strategy rule carries a
`wildcard_constraints` on the configured algorithms so the two cannot be
ambiguous.

`rule clock_estimates` is one rule rather than one per arm, because the CLI
accumulates into a single CSV and replaces any existing row for the same
(arm, quantity) -- so there is no ragged per-arm concatenation to do. The
builds are inputs only when the mode needs them: `truth` is tree-free and
reads the painted sequences directly.

**A bug this turned up:** an unrecognised `clock.mode` left `rule all` asking
for `clock_estimates.csv` while creating no rule to produce it -- a DAG error
with a confusing message. The mode-to-invocations mapping is now module-level
and shared by both, and an unknown mode raises at parse time naming the valid
ones. Verified for all three modes that the requested targets and the rules
that exist agree.

---

## The painter no longer emits the reference record

`--emit_reference` (off by default) restores it.

It was written only to the unfiltered output set, never to the
`--linelist-filter` sets, so the painted FASTA carried a record with no row in
the TwinSampler line list -- which `phylogas check-join` correctly reports as
one unmatched sequence on the `all_infections` arm. Removing it makes that
join clean.

Nothing needed it:

| worry | why it holds |
|---|---|
| does Nextstrain still get a root? | yes -- its own `reference_data` input carries `Wuhan/Hu-1/2019`, which is what `refine.root` names. The emitted record was named after the header of `data/reference/reference.fasta` (`Wuhan-Hu-1/2019`, hyphenated), so it never satisfied the root anyway |
| can mu_truth still find the reference? | yes -- it reads `genetic_painter.reference_fasta`, which is committed to this repository, not the FASTA's first record |
| can we get it back? | `--emit_reference`, or `genetic_painter.emit_reference: true` |

Correcting something stated earlier in this file's history: `--reference` is
**not** the ancestral genome the painter mutates from. `align_ref` is used only
in the emission block; the ancestral sequences come from `seed_fasta`. With
emission off, `--reference` matters only to the clock metrics.

---

## The reference id now matches ncov's root (2026-10-03)

`data/reference/reference.fasta` was headed `Wuhan-Hu-1/2019`. ncov's
`refine.root` names `Wuhan/Hu-1/2019`. One character, and it meant the
reference PhyloGAS ships could never serve as a tree root.

Renamed to the slash form. The sequence is unchanged (md5 `105c82802b67`,
29,903 bp) -- only the header.

What this buys: with `genetic_painter.emit_reference: true` the painted FASTA
is **self-rooting**. augur or TreeTime can root it directly, with no
Nextstrain installation and no separate reference input, which is what makes
the mutation-rate work possible standalone.

### Three consequences, all handled

**1. The two reference copies are the same sequence, masked differently.**
`cfg/nextstrain/references_sequences.fasta` masks the first 100 and last 50
bases to `N` -- standard ncov practice for unreliable genome termini. Ours does
not. Identical at all 29,750 unmasked positions. So they are interchangeable in
substance, but not byte-identical, and now they share a name.

**2. Emitting it would have failed the join check.** The TwinSampler line list
has no row for the reference -- it is analysis context, not a simulated
infection -- so `check-join` counted it as an unmatched sequence, which
`nextstrain-config` treats as a blocking error. Both now exempt the ids in
`cfg/nextstrain/reference_id.txt`, which is the same file ncov's subsampling
`include` reads, so the two cannot drift. `--ignore-ids` overrides.

```
Matched  : 2 / 2  (100.00%)
Context  : 1 exempt (reference/clade anchor, no line-list row by design): Wuhan/Hu-1/2019
```

A genuine mismatch still fails, and the denominator excludes the exempt
records rather than flattering the percentage.

**3. Emitting it into an ncov build duplicates the root.** Both inputs would
carry `Wuhan/Hu-1/2019`; ncov dedupes by strain and keeps one, so which
masking survives is not obvious. `nextstrain-config` now warns when it sees
the root in an arm's FASTA and points at `emit_reference: false`.

So: **off when running through ncov** (reference_data supplies a
terminus-masked root), **on for standalone** augur/TreeTime work. Either way
the clock metrics read `reference_fasta` directly, so `mu_truth` is unaffected.

Also renamed in the gitignored `data/example_data/reference.fasta`, for local
runs that point there. Historical files -- `experiments/*.snakemake.log`,
notebook outputs, BeyondBaseline's `development_dialog.json` -- keep the old
spelling as records of past runs.

---

## The reference is now byte-identical to ncov's (2026-10-03)

`data/reference/reference.fasta` carries ncov's exact
`references_sequences.fasta` record: the same 29,903 bases with the first 100
and last 50 masked to `N` (md5 `bdb4ec6a5b30`). Previously ours was unmasked
(`105c82802b67`), identical at every unmasked position but not byte-equal.

So emitting the reference into an ncov build is now harmless rather than
merely warned about -- whichever copy ncov keeps after deduping by strain is
the same sequence. `nextstrain-config` says so instead of cautioning.

### This changed the rate arithmetic

`hamming_to_reference` skips any position that is not an unambiguous base on
both sides, so divergence is measured over the **29,750 unmasked sites**, not
29,903. The per-site rate therefore divides by the comparable-site count:
dividing by the raw genome length would deflate it by 0.51% -- a systematic
error, not noise.

`comparable_sites` is now a function in `benchmark/clock.py`, and both numbers
appear in the result row. Verified: the unit conversion round-trips exactly,
and a synthetic run at a known 8.0e-4 recovers 7.85e-4 over a 300-day span
(within seed-to-seed variation; an earlier run with different noise gave
7.94e-4).

### Doc correction

`docs/clock_modes.md` claimed `mu_placed` -- the painter's own count of
substitutions placed, as a saturation diagnostic -- was computed alongside
`mu_truth`. It is not. Implementing it needs per-record accounting the painter
does not do, so it is now marked **not implemented** in both the quantity
table and comparison row 5.

## Install guidance: two layers, both required (2026-10-03)

The README's "Optional: Nextstrain" section was correct but understated, and
the failure mode it did not guard against is the common one: install the CLI,
watch `nextstrain check-setup` pass, and then be baffled that there is no
workflow to run. `check-setup` selects and validates a **runtime**; it never
obtains a **pathogen workflow**. ncov ships no `nextstrain-pathogen.yaml`, so
`nextstrain setup` will not fetch it either. The clone is mandatory.

It is also not a formality. Checked against a local ncov clone:

| Inherited from ncov | Extent |
| --- | --- |
| `files` entries | `base.yaml` overrides 2 of 13 -- we take `include`, `exclude` (360 KB), `reference_seq.gb`, `reference_seq.fasta`, `annotation.gff`, `color_ordering.tsv` (669 KB), `color_schemes.tsv` (3.9 MB), `description.md`, `clades.tsv`, `clade_display_names.yml`, `sites_ignored_for_tree_topology.txt` |
| top-level config keys | 15 of 20 untouched (`ancestral`, `genes`, `mask`, `tree`, `priorities`, `frequencies`, `nextclade_dataset`, `sanitize_metadata`, `strip_strain_prefixes`, ...) |
| `defaults/` | 5.6 MB |
| scripts and rules | 32 scripts plus the workflow itself |

So `cfg/nextstrain/ncov/base.yaml` is a deviations file in the strict sense,
and it only works because ncov loads the user `--configfile` first and merges
`defaults/parameters.yaml` underneath it.

### Three things now stated that were not

- **`ambient` vs `conda`.** `ambient` is "whatever environment Nextstrain CLI
  is itself running in" -- user-managed. Someone who ran
  `nextstrain setup conda` wants `conda`; `ambient` is for someone who
  conda-installed `augur`/`auspice` themselves. Relevant for cluster testing,
  where conda is the realistic runtime.
- **`nextstrain.runner: "snakemake"`.** Already implemented (Snakefile:338,
  `cli.py:588`) but undocumented. Bypasses the CLI entirely and calls
  `snakemake` in the ncov directory, for an already-activated environment or a
  site Snakemake profile.
- **Build-time network.** `workflow/snakemake_rules/main_workflow.smk:473` runs
  `nextclade dataset get --name sars-cov-2`. A compute node without egress will
  fail there, so the dataset has to be fetched on a login node first. This is
  not a PhyloGAS requirement and cannot be configured away from our side.

Recommended sequence, now in the README: installer -> `nextstrain setup
<runtime>` -> `check-setup --set-default` -> `git clone ncov` -> set
`nextstrain.dir`.

## The three conda environments, and the snakemake runner fix (2026-10-03)

Writing the install guidance surfaced a question we had not answered: the
instructions now mention two conda environments (PhyloGAS's and the one
`nextstrain setup conda` builds), so should they be unified? No -- and
checking ncov turned up a third, plus a bug.

| # | Environment | Owner | Holds |
| --- | --- | --- | --- |
| 1 | PhyloGAS | us | `phylogas`, pandas, scipy, sklearn, networkx, snakemake |
| 2 | `$NEXTSTRAIN_HOME/runtimes/conda/env` | Nextstrain CLI | augur, auspice, nextclade |
| 3 | `<ncov>/workflow/envs/nextstrain.yaml` | ncov | `augur=22.4.0`, `nextclade=3.9.0`, `iqtree=2.2.0.3`, `epiweeks=2.1.2` |

Environment 3 was the one we had missed. Every rule in ncov's
`main_workflow.smk` carries `conda: config["conda_environment"]`, which
`ncov/Snakefile:132` absolutizes against the Snakefile's own directory. Those
directives are **inert unless snakemake is run with `--use-conda`**, and
`nextstrain build` does not pass it -- so under the CLI route environment 2
serves augur and environment 3 is never materialized.

Not merging is the right call: environment 2 is a managed artifact that
`nextstrain update` rebuilds, and PhyloGAS never imports augur -- it shells out
to `nextstrain build`. A subprocess boundary gains nothing from a shared
environment and would force augur's numpy/pandas/biopython pins into the same
solve as scipy and sklearn. The `ambient` runtime is the supported way to have
one environment, and it is the trade we decline.

### Two bugs in the snakemake runner

The `runner: "snakemake"` branch had never worked:

| Bug | Why |
| --- | --- |
| missing `--directory` | ncov writes `auspice/{prefix}_{build}.json` and `results/...` relative to **CWD**, not to the Snakefile. Run from the PhyloGAS project directory, ncov's outputs scatter into our tree and the rule's `cp {ns_dir}/auspice/*_{arm}.json` finds nothing. `nextstrain build <dir>` chdirs for us, which is why only this branch was affected |
| missing `--use-conda` | without it the `conda:` directives are ignored and `augur`/`nextclade`/`iqtree` must be on the `PATH` of the PhyloGAS environment -- i.e. the branch silently *required* the env merge we just argued against. With it, ncov's own pins are used and environment 1 needs only snakemake plus conda |

The `--configfile` is also now passed as `{ns_dir}/{rel_config}`, since
`--configfile` resolves against the invocation CWD rather than `--directory`.

### New key: `nextstrain.conda_prefix`

`--use-conda` re-solves environment 3 per working directory, which is a
multi-GB solve we do not want repeated per arm. The key maps to
`--conda-prefix` (the same thing `SNAKEMAKE_CONDA_PREFIX` sets); empty means
let snakemake decide (`.snakemake/conda` under the ncov checkout). Only read by
the snakemake runner.

`NEXTSTRAIN_HOME` is documented in the README for the same class of reason --
environment 2 is multi-GB and `~/.nextstrain` is the wrong place for it on a
cluster with a quota'd or non-shared `$HOME`.

### README trimmed to the user's view (2026-10-03)

The previous README pass explained the three environments, the merge
argument and the `--use-conda` mechanics -- the reasoning above. That is
maintainer material; a user needs to know which commands to run. The README
section now says only: you stay in the PhyloGAS environment; clone ncov; pick
the CLI (Option A) or snakemake-only (Option B); two cluster notes. 137 lines
to 57. The reasoning stays here.

Two corrections made on the way:

- The nextclade dataset is an output file
  (`data/sars-cov-2-nextclade-defaults.zip`), so it downloads on the first
  build and is reused -- not on every build, as the previous text said.
- `phylogas status` ignored `nextstrain.runner` and reported the CLI as
  missing under `runner: "snakemake"`. It now checks for `snakemake` and
  `conda`/`mamba` in that case, and its setup hint is the one-step
  `nextstrain setup --set-default conda`.

## Verification run: DAG builds, one flaky test fixed (2026-10-03)

`snakemake -n` on the cluster builds the DAG cleanly -- 6 jobs
(`simulate_linelist`, `sample_scenarios` x2, `subset_fasta` x2, `all`), the
`{algo}` wildcard expanding over `surs` and `lasso_greedy` as configured. The
~150 Snakefile lines added over this session had never had a DAG built before;
no syntax or wildcard errors. The Nextstrain rules are absent from the DAG
because the test config has `nextstrain.enabled: false`, so `nextstrain_build`
and `clock_estimates` are still unexercised.

The "missing provenance/metadata" note for `paint_network` and `train_entropy`
is expected: those outputs were produced by earlier manual runs, so Snakemake
has no recorded provenance to compare against. Not an error.

### `test_rate_limited_model_produces_valid_sequences` was flaky, not broken

`pytest tests` gave 66 passed, 1 failed, with
`implausible mutation frequency: 0/300`. The model is fine; the test was.

The mutation count is `Poisson(mutation_rate_per_cycle * n * cycles)`. With the
test's `n = 2000`:

| quantity | value |
| --- | --- |
| `_lambda_per_cycle` | `3.40e-6 * 2000` = 0.0068 |
| `cycles` (burst 10..1000) | 3 for burst<=31, 2 for 32..999, 1 for 1000; mean 2.02 |
| mean mutations per call | 0.0137 |
| expected mutating calls | ~4 of 300 |
| **P(zero across 300)** | **1.6%, i.e. ~1 run in 62** |

Confirmed both analytically and by running the test shape across 40 seeds: 1 of
40 came back zero. The failure was the 1-in-62.

Two causes, both fixed:

- **Nondeterminism.** The test seeded a local `default_rng`, which only built
  the sequence. The model reads two *globals* -- `random.randint` for burst
  size and `np.random.randint` to seed its own Generator -- so results varied
  run to run. Both are now seeded.
- **No statistical power.** The thresholds do not affect the mutation *count*
  (they only pick which site is hit), so `n` is the only lever. The test now
  uses the real genome length, 29903, where the mean is ~0.21 per call and ~54
  of 300 calls mutate. Across 40 seeds: min 39, max 70, no zeros.

The `~82% of transmissions draw zero mutations` comment in `rate_limited.py`
refers to the full genome and checks out: `exp(-3.40e-6 * 29903 * 2.02)` =
0.815. The toy `n` in the test was what pushed that to 98.6%.

## First real clock measurement, and the span it reported (2026-10-04)

`phylogas benchmark clock --mode truth` on the painted VA Delta genomes:

```
mu_truth  8.8734e-04 subs/site/year   R2=0.124  n=20,000  span=537d
          slope=0.0723 subs/genome/day   mean_divergence=45.13
```

The rate is good -- ~11% above the canonical ~8e-4 -- and it is internally
consistent with the painter's own parameters, which is the part worth keeping.
The painter draws `Poisson(3.40e-6 * 29903 * 2.02)` = 0.206 mutations per
transmission, so 26.4 subs/genome/year needs ~128 transmissions/year along a
chain, implying a **2.84-day generation interval**. That is right for Delta.
The mutation engine and the epidemic dynamics agree without being fitted to
each other.

`mean_divergence = 45.1` is the Delta-from-Wuhan offset. It inflates the
intercept, not the slope, so it does not bias the rate.

### But span=537 described a 70-day study

The simulation runs 2021-04-07 to ~2021-06-15. The 537 days came from one
record dated 2019-12-26 -- Wuhan-Hu-1's real collection date. The run predated
the `emit_reference` gate, so the reference was in the painted output.

This matters because `date_span_days` is the field the docs tell you to check
to judge a fit, and it was describing a period the fit barely covered. The
rate itself was barely affected: at n=20,000 the anchor carries only 3% of the
x-variance, so it inflated the slope ~0.7%.

Fixed in `root_to_tip`, which now takes `date_min`/`date_max` and
`exclude_ids`:

| Default | Source | Why |
| --- | --- | --- |
| window `>= start_date` | `genetic_painter.start_date` | records outside the simulated period are not observations of its clock |
| drop reference ids | `cfg/nextstrain/reference_id.txt` (via `_context_ids`) | divergence from itself is structurally 0, at a date 500 days early |

Neither applies to `--kind real`. Overrides: `--date-min`, `--date-max`,
`--exclude-id`, `--keep-context`. New columns: `skipped_out_of_window`,
`skipped_excluded_id`, `window`, `first_date`, `last_date`.

Validated against synthetic data at a known 8.0e-4:

| window | anchor left in | fix applied |
| --- | --- | --- |
| 400 days | 7.97e-4 (-0.4%), R2 0.862, span 868d | 7.94e-4 (-0.8%), R2 0.861, span 400d |
| 70 days | 8.46e-4 (+5.8%), R2 0.212, span 538d | 7.87e-4 (-1.7%), R2 0.163, span 70d |

The 70-day row reproduces the real run closely (R2 0.16 vs 0.124, span 538 vs
537), which confirms the diagnosis -- and shows the low R2 is a property of
the short window, not a fault. Over 70 days the clock contributes only ~5
substitutions against a comparable seed-to-seed spread.

### Two smaller fixes in the same pass

- `benchmark clock --config` defaulted to `None` while `status` defaults to
  `config.yaml`. Run without it, every lookup fell through to a hardcoded
  default and the CSV landed in `results/05_benchmarks/` instead of the
  project's `results/<name>/05_benchmarks/` -- where the Snakefile's
  `clock_estimates` rule expects it. Now defaults to `config.yaml`; an
  explicitly named missing file is still an error, the default degrades to no
  config so the command still works outside a project.
- `--max-records` help now says it truncates in file order rather than
  sampling. It breaks out of the read loop, so on a tick-ordered FASTA it
  biases toward early infections.

## `21L` removed from the reference context (2026-10-04)

`cfg/nextstrain/references_*` and `reference_id.txt` carried a second record,
`21L` (Omicron BA.2, dated 2021-11-01), added as a "clade anchor" when these
files came over from the Omicron-era ncov configs (see the `reference_id.txt`
row above).

In a Delta-wave build it does not belong. Neither reference record has a
`division`, so `strategy_focal_context`'s context query
`division != '{division}'` selects both, and the tree gets a lone BA.2 tip five
months after the simulation ends -- distorting the topology and pulling on the
`inferred` clock rate.

Removed from all three files. The root needs nothing extra: ncov's
`combine_samples` always appends `files.include` (`defaults/include.txt`,
which lists `Wuhan/Hu-1/2019`), and the Wuhan record is unchanged and still
byte-identical to `data/reference/reference.fasta` (md5 `bdb4ec6a5b30...`).

If a future build wants a clade anchor, it should belong to the variant being
simulated; the comment in `reference_id.txt` now says so.

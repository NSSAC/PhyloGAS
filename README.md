# PhyloGAS
## PhyloGeographic Analysis Similars


[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.github.io)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11-blue.svg)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**PhyloGAS** is a digital-twin-based framework for benchmarking public health genomic surveillance strategies. 

Genomic surveillance faces constrained sequencing capacity and ascertainment biases, distorting the continuous global transmission process into fragmented clusters. Because the "ground truth" of an outbreak is rarely known in the real world, quantifying the performance of deployed phylodynamic reconstruction methods remains a significant challenge. 

PhyloGAS bridges this gap by integrating **agent-based epidemiological simulations** with a **calibrated mutation engine**, allowing public health officials and researchers to test their phylodynamic pipelines against a known, synthetic ground truth.

---

## 🧬 Pipeline Architecture

PhyloGAS is orchestrated via **Snakemake**. The pipeline consists of five primary stages:

1. **Data Acquisition:** Fetches synthetic population demographics from the UVA Dataverse and EpiHiper simulation replicates from Zenodo, MD5-verified.
2. **The Genetic Painter (`src/phylogas/painter/`):** Overlays biologically plausible SARS-CoV-2 genomes onto the transmission tree, starting from real imported sequences (the seed FASTA) rather than the reference. A "Speedometer" (calibrated Poisson rate-limiting) dictates *how many* mutations occur per transmission; a "Map" (MSA-derived entropy and substitution matrices) dictates *where* they fixate.
3. **Ascertainment Simulation ([TwinSampler](https://github.com/NSSAC/TwinSampler)):** Applies severity- and demographics-dependent detection to every infection, producing a biased line list plus the full infection record. It models *who* gets detected, not reporting delay: dates are symptom/state onset.
4. **Adaptive Sampling ([BeyondBaseline](https://github.com/NSSAC/BeyondBaseline)):** Selects which cases to sequence under a weekly budget — SURS, stratified, greedy and LASSO-based strategies, across stride / pool-window / target scenarios. Each (scenario × algorithm) pair is a **recipe**, e.g. `4S__surs` or `4S-4_LL-P__lasso_greedy`.
5. **Reconstruction & Benchmarking:** Builds a Nextstrain (augur/TreeTime) tree per chosen recipe, then scores against the ground truth: geographic flow (cosine similarity, topological F1), variant prevalence, painted-sequence sanity, and the molecular clock (the rate the painter actually produced vs. the rate Nextstrain assumed or inferred — see [`docs/clock_modes.md`](docs/clock_modes.md)).

---

## 🧪 Beta testers

Start with **[`BETA_TESTING.md`](BETA_TESTING.md)** — install, first run,
known gaps, and what to report.

---

## 🛠️ Installation

PhyloGAS relies on a master environment file that automatically installs its dependencies, including satellite repositories, via Git URLs. 

### Option A: Conda/Mamba (Recommended for Developers)
```bash
git clone https://github.com/NSSAC/PhyloGAS.git
cd PhyloGAS

# Create the environment and install all dependencies
mamba env create -f environment.yml
conda activate phylogas_env

# Install PhyloGAS itself (editable)
pip install -e .
```
*(`environment.yml` also pulls `TwinSampler`, `BeyondBaseline` and
`pango_aliasor` directly from GitHub. If that pip step fails, report it rather
than working around it.)*

### Option B: pip only
```bash
pip install -e ".[all]"     # or ".[fast]" for just pyarrow
```

Verify the install:
```bash
phylogas --version
phylogas --help
```

### Option C: Docker / Apptainer
> **Status: not yet available.** The `Dockerfile` has not been written and no
> image is published. Use Option A or B for now.

### Optional: Nextstrain, for the phylodynamic stage

Only needed for `nextstrain.enabled: true`.

You keep working in your PhyloGAS environment the whole time. Nextstrain's own
tools (augur, nextclade, …) live in a separate environment that is built for
you and switched into automatically whenever a build runs — you never activate
it, and you should not install augur into the PhyloGAS environment.

**1. The ncov workflow clones itself.** There is no path to configure and
nothing to clone by hand. On first use `phylogas run` (or `phylogas
ncov-checkout`) clones ncov into the project's own results directory at a
pinned commit:

```yaml
nextstrain:
  enabled: true
  # The defaults, shown for reference -- you do not need to set any of these.
  dir: "{results_dir}/ncov"
  repo: "https://github.com/nextstrain/ncov.git"
  ref: "3432c85760b6c8cd0f815a85167ebb46c8131abb"
```

**One checkout per project, not one shared.** ncov writes
`results/combined_metadata.tsv.xz`, `results/combined_sequences_for_subsampling.fasta.xz`
and `results/combined_sequence_index.tsv.xz` at fixed paths that every build in
a checkout shares. Two projects building in one checkout at the same time
therefore overwrite each other's inputs, and the symptoms are nasty: a lock
error if you are lucky, another state's tips in your tree if you are not. A
checkout is ~25 MB, so a copy per project is the cheap way out. It stays until
you delete the project's results, because the clock benchmark reads the build's
`results/`.

`ref` is pinned so every project, and every rerun, builds with the same ncov.
Bump it deliberately; existing checkouts are left alone, so delete them (or the
project's results) to pick up a new ref.

*No internet on the machine?* Clone ncov once by hand and point `source_dir` at
it. It is cloned from there instead, and its already-downloaded Nextclade
dataset is copied across so the first build needs no network either:

```yaml
nextstrain:
  source_dir: "/shared/ncov"
```

**2. Choose how builds run.** Pick one.

*Option A — Nextstrain CLI (default, recommended).* Install it once; it builds
and maintains the tools environment.

```bash
curl -fsSL --proto '=https' https://nextstrain.org/cli/installer/linux | bash   # or .../mac
nextstrain setup --set-default conda      # use docker instead if you have it
```

*Option B — Snakemake only, no CLI.* Snakemake builds the tools environment
from ncov's own recipe on the first run. Needs `conda` or `mamba` on your
`PATH`.

```yaml
nextstrain:
  runner: "snakemake"
```

**3. Check.** `phylogas status` reports whether the workflow and the runner are
in place and what is missing.

#### On a cluster

- **Put the tools environment somewhere roomy.** It is several GB. Option A:
  `export NEXTSTRAIN_HOME=/scratch/$USER/nextstrain` before running setup.
  Option B: set `nextstrain.conda_prefix: "/scratch/$USER/snakemake-conda"`,
  which also means it is built once rather than per run.
- **The clone and the first build need internet.** The checkout is cloned from
  GitHub, and ncov then downloads a Nextclade dataset into it on first use and
  reuses it afterwards. If compute nodes have no outbound access, either run
  `phylogas ncov-checkout --config <config>` and the first build from a login
  node, or set `nextstrain.source_dir` to a checkout you cloned earlier.
- **Budget disk for a checkout per project.** ~25 MB of code, plus the
  Nextclade dataset and whatever ncov's own `results/` grows to for each
  build.

`phylogas nextstrain-config` (run for you by `phylogas run`) writes one ncov
config per **recipe** listed in `sampling.nextstrain_recipes` and copies its
inputs into the project's ncov checkout under
`data/phylogas/<project_name>/<recipe>/`, alongside a `_shared/` directory for
the files every recipe uses (lat/longs, the Auspice config, the reference). It
stops with an error if any sequence lacks a metadata row. To run it by hand:

```bash
phylogas nextstrain-config --config config.yaml --build-type strategy --recipe 4S__surs
```

Builds use a **fixed** clock rate (`clock_rate: 0.0008` in
`cfg/nextstrain/ncov/base.yaml`, ncov's SARS-CoV-2 default) and
`coalescent: opt`. The painter's measured rate is lower (about 2.8–2.9e-4 on
the Virginia Delta run), so dated trees are scaled to the assumed rate, not the
simulated one. `nextstrain.clock` controls the comparison; see
`docs/clock_modes.md`.

---

## 🚀 Quickstart & Usage

PhyloGAS is driven entirely by a configuration file: `config.yaml`.

**1. Configure your run:**
```bash
cp config.template.yaml config.yaml
$EDITOR config.yaml
```

**2. Ask what to do next — at any point:**
```bash
phylogas status
```
```text
PhyloGAS pipeline status
====================================================
  [okay] Configuration            config.yaml
  [MISS] Synthetic population     data/va/va_2_4_0_demographics.csv
  [okay] Entropy map                   3.4 MB  data/run.03.base.threshold.df.npy
  [okay] Seed sequences               12.5 MB  data/Virginia_B_1_617_2_seed_sequences.fasta.gz
  [MISS] EpiHiper output          data/example_data/.../output.csv.gz
  [MISS] Painted genomes          results/01_synthetic_genomes/va_delta_wave.sequences.fasta.xz

Next step:
  phylogas fetch-data --states va
```
`status` reports every stage and prints the single next command. It is the
intended entry point — run it whenever you are unsure where you are.

**3. Acquire the data:**
```bash
phylogas fetch-data --states va --dry-run   # preview: what, and how large
phylogas fetch-data --states va             # download + build demographics
```
Downloads the four synthetic-population files for the state (~0.10 GB for VA)
and joins them into the demographics table the painter reads. The 907 MB
EpiHiper contact network is **not** fetched by default — nothing in PhyloGAS
reads it. Add `--with-epihiper-inputs` if you intend to run the ABM yourself.

**4. Execute the pipeline:**
```bash
phylogas run --config config.yaml --cores all
phylogas run --config config.yaml --cores all --until paint_only
phylogas run --config config.yaml --profile slurm
phylogas run --config config.yaml --dry-run
```

**Or drive a single stage directly:**
```bash
phylogas train --config config.yaml        # Stage 1: entropy map from an MSA
phylogas paint --config config.yaml        # Stage 2: paint the network
phylogas subset-fasta -m samples.csv.xz -f full.fasta.xz -o subset.fasta.xz
```
Any config value can be overridden with a flag, and `--dry-run` prints the
underlying command:
```bash
phylogas paint --config config.yaml --num-ticks 30 --dry-run
```

Migrating from the older standalone scripts? See **[`command_map.md`](command_map.md)**.

---

## ⚙️ Configuration (`config.yaml`)

The pipeline is driven by one YAML file. Start from
[`config.template.yaml`](config.template.yaml), which is commented throughout and
is the authoritative reference for every key.

Strings may reference other keys with `{dotted.path}`; references are expanded
after loading, so ordering does not matter:

```yaml
project_name: "va_delta_wave"
data_dir: "data"
results_dir: "results"

population:
  state: "va"
  persontrait_file: "{data_dir}/va_2_4_0_demographics.csv"

genetic_painter:
  mutation_model: "rate_limit"      # rate_limit | simple | poor
  initial_viral_load: 10            # transmission bottleneck
  entropy_thresholds: "{data_dir}/run.03.threshold.file"
  probability_matrix: "{data_dir}/run.03.base.threshold.df.npy"
  output_prefix: "{results_dir}/01_synthetic_genomes/{project_name}"
  start_date: "auto"                # = ABM tick 0 + start_tick (abm.config)
  start_tick: 128
  num_ticks: 300

  # xz level 1 is ~37x faster than the lzma default (6) and still reaches
  # ~320x compression on this data. Level 6 adds hours to a 5M-record run.
  compression: "xz"
  compression_level: 1
  compression_threads: 8

sampling:
  algorithms: ["surs", "lasso_greedy"]   # recipe slugs; BeyondBaseline names also accepted
  batch_size: 400
  nextstrain_recipes: ["4S__surs"]       # which recipes get a tree; "all" for every one

nextstrain:
  enabled: true        # the ncov checkout is cloned for you; see above
```

**Recipes.** BeyondBaseline runs every scenario for each algorithm in
`sampling.algorithms`; `sampling.nextstrain_recipes` picks which of those
sample sets get subset and built. List the valid ids with `scenarios-recipes`.
SURS ignores the target, so its recipes are named by stride only (`1S__surs`,
`4S__surs`); the older `4S-4_LL-P__surs` spelling still resolves.

**Changing an upstream parameter.** Snakemake reruns a rule when its code or
params change, but not when you delete an intermediate whose downstream
targets already exist. To force a stage and everything after it, call
Snakemake directly (add `-n` first to see what would rerun):

```bash
snakemake --configfile config.yaml --cores all --forcerun simulate_linelist
```

---

## 📂 Output Structure

Upon successful completion, PhyloGAS generates a structured `results/` directory:

```text
results/<project_name>/
├── 00_mutation_model/       # entropy thresholds + substitution matrix (phylogas train)
├── 01_synthetic_genomes/    # painted FASTA + metadata for every infection (xz)
├── 02_simulated_linelists/  # linelist.csv.xz (ascertained), linelist_allevents.csv.xz (all infections)
├── 03_sampled_datasets/     # <recipe>_samples.csv.xz, <project>.<recipe>.fasta.xz
├── 04_nextstrain_builds/    # <recipe>/auspice.json, one per nextstrain_recipes entry
├── 05_benchmarks/           # clock_estimates.csv, Mugration_Metrics.csv,
│                            # AUC_truth_rankings.csv, sequence_divergence.csv, ...
└── ncov/                    # this project's ncov checkout (cloned, pinned);
                             # holds the staged inputs and ncov's own results/
```

`rule all` produces the painted set, the sampled subsets and (with
`nextstrain.enabled`) the trees and `clock_estimates.csv`. The other benchmark
CSVs are separate targets — name them, or use `phylogas benchmark ...`.

---

## 🗂️ Repository Layout

```text
PhyloGAS/
├── Snakefile                  # Pipeline DAG
├── config.template.yaml       # Commented configuration reference
├── command_map.md             # Old script invocations -> phylogas CLI
├── pyproject.toml             # Package metadata; defines the `phylogas` command
├── environment.yml            # Conda environment incl. satellite repos
├── scripts/legacy/            # Superseded scripts, kept for provenance
├── src/phylogas/
│   ├── cli.py                 # Unified command line interface
│   ├── config.py              # YAML loading + {placeholder} expansion
│   ├── dataverse.py           # Per-state DOIs + verified file-ID table
│   ├── painter/               # The Genetic Painter
│   │   ├── genetic_painter.py
│   │   ├── mutational_models/ # rate_limited | simple | poor
│   │   └── test/              # SLURM submission scripts
│   ├── seqprep/               # Seed acquisition, FASTA subsetting
│   ├── popprep/               # Demographics assembly (persontrait join)
│   └── benchmark/             # Scoring against ABM ground truth
│       ├── mugration.py       #   parsimony inference on the transmission tree
│       ├── scoring.py         #   cosine / F1 / pearson / masked MAE
│       ├── truth_metrics.py   #   weekly infections / variant counts from allevents
│       ├── clock.py           #   mu_truth, mu_sim, mu_operational
│       ├── lineage.py         #   within- vs between-chain clock decomposition
│       └── runner.py          #   batch drivers
├── docs/
│   ├── clock_modes.md         # what each clock rate measures, and how to read them
│   ├── design/                # open design questions (e.g. metadata_ownership.md)
│   └── staging_local_data.md  # using files already on a cluster
└── cfg/
    ├── nextstrain/ncov/       # ncov build overrides (clock, coalescent)
    └── ...                    # EpiHiper experiment configs
```

---

## ✅ Implementation Status

This project is mid-restructure. What is actually wired up today:

| Stage | Command | Status |
|---|---|---|
| — | `phylogas status` | Works — **start here** |
| — | `phylogas fetch-data` | Works — verified against live Dataverse |
| — | `phylogas build-demographics` | Works — verified on 5.2M MN rows |
| — | `phylogas ncov-checkout` | Works; clones `nextstrain.dir` at the pinned `ref`, no-op if present |
| 0. Seed acquisition | `phylogas prep-seeds` | Works |
| 1. Entropy training | `phylogas train` | Works, verified end-to-end |
| 2. Genetic painting | `phylogas paint` | Works, verified end-to-end |
| 3. Ascertainment | rule `simulate_linelist` (TwinSampler) | Works, verified at cluster scale |
| 4. Adaptive sampling | rule `sample_scenarios` (`beyond-baseline-sweep`) | Works; one sample set per recipe (`sampling.replicates` not yet wired) |
| 5. FASTA subsetting | `phylogas subset-fasta` | Works |
| 6. Nextstrain | rule `nextstrain_build` | Works, opt-in; verified with the Nextstrain CLI on a cluster, 4 states concurrently |
| 7. Clock benchmark | `phylogas benchmark clock` | Works for `truth` and `operational`; `inferred` (mu_sim) needs augur on `PATH` |
| 7. Other benchmarks | `phylogas benchmark truth / mugration / sequence` | Work — mugration scoring verified bit-identical to the pre-split implementation |
| — | `phylogas compare-strategies` | Works — runs a BeyondBaseline sweep, then ranks every strategy |
| — | Docker / Apptainer | **Not written** |

Verified end to end on the cluster: the Virginia Delta replicate (ticks
128–428, 5.35M painted infections) through painting, ascertainment, sampling,
a `4S__surs` Nextstrain build and `clock_estimates.csv`.

---

## 📚 Data Availability & Acknowledgements
Due to size constraints, the synthetic population data and EpiHiper outputs are not hosted in this repository. `phylogas fetch-data` pulls the population files from the UVA Dataverse ([doi:10.18130/V3/5LSDCY](https://dataverse.lib.virginia.edu/dataset.xhtml?persistentId=doi:10.18130/V3/5LSDCY)) and the simulation replicates from Zenodo ([doi:10.5281/zenodo.23067670](https://doi.org/10.5281/zenodo.23067670)).

**Core Components:**
* **Genetic Painter:** Integrated within `src/phylogas/painter/`.
* **Twin Sampler:** [github.com/NSSAC/TwinSampler](https://github.com/NSSAC/TwinSampler)
* **BeyondBaseline:** [github.com/NSSAC/BeyondBaseline](https://github.com/NSSAC/BeyondBaseline)

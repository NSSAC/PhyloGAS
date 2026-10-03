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

1. **Data Acquisition:** Automatically fetches multi-gigabyte synthetic population demographics and EpiHiper transmission dendrograms from the UVA Dataverse.
2. **The Genetic Painter (`src/genetic_painter/`):** Overlays biologically plausible SARS-CoV-2 viral genomes onto the transmission tree. It utilizes a two part architecture: a "Speedometer" (calibrated Poisson rate-limiting) dictates *how many* mutations occur per transmission, while a "Map" (MSA-derived entropy and substitution matrices) dictates *where* those mutations stably fixate.
3. **Ascertainment Simulation (Twin Sampler):** Simulates the real-world lag and demographic biases of infection reporting, generating a skewed, realistic line list of cases.
4. **Adaptive Sampling (Beyond Baseline):** Subsets the line list using Simple Uniform Random Sampling (SURS) and other sampling strategies (such as stratified sampling).
5. **Reconstruction & Benchmarking:** Runs the sampled sequences through standard public health phylodynamic pipelines (e.g., Nextstrain/TreeTime) and calculates Topological F1-Scores and Cosine Similarities against the absolute ground truth.

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
*(`environment.yml` also pulls the `TwinSampler` and `BeyondBaseline` libraries
directly from GitHub. Those repos each need a `pyproject.toml` for this to work
- see `docs/salvage_audit.md` if the pip step fails.)*

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
  start_date: "2021-04-07"
  start_tick: 128
  num_ticks: 300

  # xz level 1 is ~37x faster than the lzma default (6) and still reaches
  # ~320x compression on this data. Level 6 adds hours to a 5M-record run.
  compression: "xz"
  compression_level: 1
  compression_threads: 8

sampling:
  algorithms: ["surs", "lasso_greedy"]
  batch_size: 400
```

---

## 📂 Output Structure

Upon successful completion, PhyloGAS generates a structured `results/` directory:

```text
results/
├── 01_synthetic_genomes/        # Full ground-truth FASTA and Metadata (BGZF compressed)
├── 02_simulated_linelists/      # Ascertainment-biased case logs
├── 03_sampled_datasets/         # Subsets of FASTAs based on SURS, LASSO, etc.
├── 04_nextstrain_builds/        # Auspice JSONs and inferred trees from the pipeline
└── 05_benchmarks/               # Final CSVs containing Cosine Similarity & F1-Scores
```

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
│       └── runner.py          #   batch drivers
├── docs/salvage_audit.md      # What was rescued from synthetic_biosurveillance
└── cfg/                       # EpiHiper experiment configs
```

---

## ✅ Implementation Status

This project is mid-restructure. What is actually wired up today:

| Stage | Command | Status |
|---|---|---|
| — | `phylogas status` | Works — **start here** |
| — | `phylogas fetch-data` | Works — verified against live Dataverse |
| — | `phylogas build-demographics` | Works — verified on 5.2M MN rows |
| 0. Seed acquisition | `phylogas prep-seeds` | Works |
| 1. Entropy training | `phylogas train` | Works, verified end-to-end |
| 2. Genetic painting | `phylogas paint` | Works, verified end-to-end |
| 3. Ascertainment | `simulate_linelist` (TwinSampler) | External; Snakemake rule written, untested here |
| 4. Adaptive sampling | `scenarios-runner` (BeyondBaseline) | External; requires that repo's `pyproject.toml` |
| 5. FASTA subsetting | `phylogas subset-fasta` | Works |
| 6. Nextstrain | rule `nextstrain_build` | Opt-in, untested (needs an ncov checkout) |
| 7. Benchmarking | `phylogas benchmark` | Works — mugration scoring verified bit-identical to the pre-split implementation |
| — | `phylogas compare-strategies` | Works — runs a BeyondBaseline sweep, then ranks every strategy |
| — | Docker / Apptainer | **Not written** |

Verified working: `phylogas train` → `phylogas paint` via Snakemake produces
29,010 sequences with matching metadata, and a re-run is correctly a no-op.

---

## 📚 Data Availability & Acknowledgements
Due to size constraints, the heavy synthetic population data and raw EpiHiper transmission networks are not hosted in this repository. The pipeline automatically pulls necessary files from the UVA Dataverse ([doi:10.18130/V3/5LSDCY](https://dataverse.lib.virginia.edu/dataset.xhtml?persistentId=doi:10.18130/V3/5LSDCY)). 

**Core Components:**
* **Genetic Painter:** Integrated within `src/genetic_painter/`.
* **Twin Sampler:** [github.com/NSSAC/TwinSampler](https://github.com/NSSAC/TwinSampler)
* **BeyondBaseline:** [github.com/NSSAC/BeyondBaseline](https://github.com/NSSAC/BeyondBaseline)

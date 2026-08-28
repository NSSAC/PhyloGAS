# PhyloGAS
## PhyloGeographic Analysis Similars


[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.github.io)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11-blue.svg)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**PhyloGAS** is a comprehensive, digital-twin-based framework for benchmarking public health genomic surveillance strategies. 

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

## 🛠️ Installation

PhyloGAS relies on a master environment file that automatically installs its dependencies, including satellite repositories, via Git URLs. 

### Option A: Conda/Mamba (Recommended for Developers)
```bash
git clone https://github.com/NSSAC/PhyloGAS.git
cd PhyloGAS

# Create the environment and install all dependencies
mamba env create -f environment.yml

# Activate the environment
conda activate phylogas_env
```
*(Note: `environment.yml` automatically installs the `twin_sampler` and `BeyondBaseline` libraries directly from GitHub).*

### Option B: Docker / Apptainer (Recommended for Production/HPC)
For absolute reproducibility without environment conflicts on state or university clusters, use our pre-built container:
```bash
# Example Apptainer/Singularity execution
apptainer pull phylogas.sif docker://nssac/phylogas:latest
```

---

## 🚀 Quickstart & Usage

PhyloGAS is driven entirely by a configuration file: `config.yaml`.

**1. Configure your run:**
Copy the template and modify your desired parameters, such as viral load bottlenecks, target variant waves, and sampling algorithms.
```bash
cp config.template.yaml config.yaml
vim config.yaml
```

**2. Execute the pipeline:**
Let Snakemake handle the data dependencies, multi-threading, and script orchestration:
```bash
# Run locally using all available cores
snakemake --use-conda --cores all

# OR: Run on a SLURM cluster (Snakemake natively handles job submission)
snakemake --profile slurm_profile
```

---

## ⚙️ Configuration (`config.yaml`)

The pipeline relies on a unified YAML configuration. Here is an example of the configuration structure:

```yaml
# --- PhyloGAS Main Configuration ---

# 1. Project & Dataverse
project_name: "va_delta_wave_exp7"
dataverse_doi: "doi:10.18130/V3/5LSDCY"
data_dir: "data/"  # Dataverse downloads will populate here

# 2. Digital Twin Demographics
population:
  state: "va"
  persontrait_file: "data/va_persontrait_epihiper.csv"
  household_file: "data/va_household.csv"

# 3. The Genetic Engine (Genetic Painter)
genetic_painter:
  mutation_model: "rate_limit"     # Options: rate_limit, simple, poor
  initial_viral_load: 10           # Transmission bottleneck size
  peak_viral_load: 1e9             # Drives early-phase replication cycles
  reference_fasta: "data/reference.fasta"
  entropy_thresholds: "data/run.03.threshold.file"
  probability_matrix: "data/run.03.base.threshold.df.npy"

# 4. Surveillance & Sampling
surveillance:
  ascertainment_config: "data/ascertainment_parameters.yaml"
  batch_size: 400
  algorithms: 
    - "surs"
    - "stratified"
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

## 📚 Data Availability & Acknowledgements
Due to size constraints, the heavy synthetic population data and raw EpiHiper transmission networks are not hosted in this repository. The pipeline automatically pulls necessary files from the UVA Dataverse ([doi:10.18130/V3/5LSDCY](https://dataverse.lib.virginia.edu/dataset.xhtml?persistentId=doi:10.18130/V3/5LSDCY)). 

**Core Components:**
* **Genetic Painter:** Integrated within `src/genetic_painter/`.
* **Twin Sampler:** [github.com/aswarren/twin_sampler](https://github.com/aswarren/twin_sampler)
* **BeyondBaseline:** [github.com/aswarren/BeyondBaseline](https://github.com/aswarren/BeyondBaseline)

# PhyloGAS
## Sequence Acquisition / Preparation

### `seq_prep.py`

This script is a flexible tool for acquiring FASTA sequences from the Cov-Spectrum API. It operates in two main modes, controlled by the `--seed_mode` flag.

#### Mode 1: Seed Finding Mode (`--seed_mode`)

This mode replicates the logic of `seed_seq_prep.py`. It uses a cluster data file (from the UCSC SARS-CoV-2 Genome Browser) to identify the earliest "seed" sample from distinct importation clusters. This is useful for creating a phylogenetically representative, non-redundant set of sequences for a given variant wave. It also generates a consolidated importation schedule CSV, which is a crucial input for `label_components.py`.

**Usage:**
```bash
python seq_prep.py \
  --state "Virginia" \
  --pango "JN.1,XBB.1.5" \
  --output_folder ./schedules_and_seeds/ \
  --seed_mode \
  --outlier_method chaining \
  --rescue_cluster_days 365
```

This command will:
- Download the cluster tracker file (if not found in the output folder).
- Identify importation clusters for JN.1 and XBB.1.5 in Virginia.
- Perform outlier removal using the `chaining` method, rescuing multi-sample clusters that span up to a year.
- Generate `Virginia_JN_1_seed_strains.txt` and `Virginia_XBB_1_5_seed_strains.txt` containing sample IDs.
- Generate a consolidated `Virginia_schedule.csv` for use in other scripts.
- Download FASTA files for the identified seeds (unless `--no_download` is specified).

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



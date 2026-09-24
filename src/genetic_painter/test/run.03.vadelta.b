#!/bin/bash

set -euo pipefail

# NOTE: this path is relative to the directory you submit from (the `test/`
# directory), and it must point at the package layout, not the old flat one.
# The previous value "../src/genetic_painter.py" no longer exists.
code="../src/genetic_painter/genetic_painter.py"

## Inputs for both analyses.

### Analysis type.
# NOTE: the valid choice is "generate_sequence". The older
# "generate_sequence_analysis" spelling was removed in commit 8de18a5
# ("simplify mode") and argparse now rejects it outright.
analysis_type="generate_sequence"

random_number_seed=43 # Using the seed from your script

threshold_file="run.03.vadelta.output/run.03.threshold.file"
# Assuming base_threshold_df for the .npy file does not need .npy extension in variable
# as the script will add it if missing.
base_threshold_df="run.03.vadelta.output/run.03.base.threshold.df.npy" 

# Fasta file of genomic sequences.
# align_fasta="va_variant_BA.2.12.1_5_sequences.fasta" # Not used if seed_fasta is provided for gen_seq analysis
seed_fasta="/project/bii_nssac/biocomplexity/vdh_genomics/synthetic_biosurveillance/SARS-Cov2-Biosurveillance-Simulation/test/run.03.vadelta.output/Virginia_B_1_617_2_seed_sequences.fasta"

## Inputs for generate sequence analysis.

start_date="2021-04-07"
#start_date="2021-06-01"

input_graph_csv="/project/bii_nssac/epihiper-simulations/pipeline-jc/run/20250120_1/output_root/proj/20250120_1/batch_1/0.25/va/replicate_0/output.csv.gz"
#input_graph_csv="/project/biocomplexity/vdh_genomics/synthetic_biosurveillance/SARS-Cov2-Biosurveillance-Simulation/data/dendrogram/epihiper_exp7_dendrogram.csv"

#output file prefix - updated to reflect new parameters
output_prefix="/project/bii_nssac/biocomplexity/vdh_genomics/synthetic_biosurveillance/PhyloGAS/test/run.03.vadelta.output/run_03_vadelta_2026_09_10_128to428" # hold.filehanged to distinguish from previous runs

# Compress to xz format
compression_type="xz"

# --- Compression tuning -----------------------------------------------------
# This is the single biggest lever on runtime. Measured on 29,903 bp records:
#     preset 6 (the lzma default)  ~243 rec/s   <-- caused the original 8.68 h run
#     preset 1                   ~9,000 rec/s
# On this highly redundant data preset 1 still achieves ~320x compression
# (preset 6 reaches ~1500x, so files are a few x larger but written ~37x faster).
compression_level=1

# Pipe through the external multi-threaded `xz` binary instead of Python's
# single-threaded lzma module. Should match the -c value in the SBATCH header.
# Set to 0 to use in-process lzma.
compression_threads=8
# ---------------------------------------------------------------------------

# Persontrait file
#persontrait_file="/project/biocomplexity/vdh_genomics/synthetic_biosurveillance/SARS-Cov2-Biosurveillance-Simulation/data/merged_population_files/va_merged_person.csv"
persontrait_file="/project/biocomplexity/vdh_genomics/synthetic_biosurveillance/SARS-Cov2-Biosurveillance-Simulation/data/merged_population_files/va_2_4_0_demographics.csv"

# Add metadata to include from the persontrait file
add_metadata="gender,county,home_latitude,home_longitude,latino,race,smh_race,age_group"

# In code, this parameter's default value is False.
# This is used in the 'better' mutational model.
# The --proportional flag is added directly in the command string below.

# If true, then use the poor mutational model; otherwise
# use the better model.
# Default value is false.
# poor="False" # Not explicitly set, relies on script default, or add --poor False if needed

# Max number of values to process.
#limit=16521
# limit=1000 # Removing --limit as requested

# Add reference sequence to the output.
reference="/project/biocomplexity/vdh_genomics/synthetic_biosurveillance/SARS-Cov2-Biosurveillance-Simulation/data/training_sequences/reference.fasta"

# Either input_graph_painted_state or input_graph_painted_prefix
# can be set; defaults to input_graph_painted_state="var1E"
input_graph_painted_prefix="E2"

# New parameters for rate limiting and tick control
enable_rate_limit="--rate_limit" # Use this to toggle the feature easily
initial_viral_load_val=10 # Default or can be set as a variable
start_tick_val=128
#num_ticks_val=30
num_ticks_val=300

# Location and Reference Location (as strings, already handled well)
location_val='{"country":"USA","division":"Virginia","divisionAbbr":"VA","region":"North America"}'
reference_location_val='{"country":"China","division":"Wuhan","divisionAbbr":"Hu","region":"Asia","date":"2019-12-26"}'

# This script only paints; it does not train. Fail early and clearly if the
# entropy outputs from run.03.vadelta.a are missing.
for f in "${threshold_file}" "${seed_fasta}" "${input_graph_csv}" "${persontrait_file}" "${reference}"; do
    [[ -r "${f}" ]] || { echo "ERROR: cannot read required input: ${f}" >&2; exit 1; }
done
[[ -s "${base_threshold_df}" || -s "${base_threshold_df}.npy" ]] \
    || { echo "ERROR: probability matrix not found: ${base_threshold_df}[.npy]" >&2
         echo "       Run run.03.vadelta.a first to generate it." >&2; exit 1; }
mkdir -p "$(dirname "${output_prefix}")"

## Execute.
# Built as an array rather than a string so that paths containing spaces cannot
# be word-split and `eval` is not required.
command=(
    python -u "${code}"
    --analysis_type              "${analysis_type}"
    --random_number_seed         "${random_number_seed}"
    --threshold_file             "${threshold_file}"
    --base_threshold_df          "${base_threshold_df}"
    --start_date                 "${start_date}"
    --input_graph_csv            "${input_graph_csv}"
    --seed_fasta                 "${seed_fasta}"
    --output_prefix              "${output_prefix}"
    --compression                "${compression_type}"
    --compression_level          "${compression_level}"
    --compression_threads        "${compression_threads}"
    --persontrait_file           "${persontrait_file}"
    --add_metadata               "${add_metadata}"
    --input_graph_painted_prefix "${input_graph_painted_prefix}"
    --proportional
    --reference                  "${reference}"
    --location                   "${location_val}"
    --reference_location         "${reference_location_val}"
    ${enable_rate_limit}
    --initial_viral_load         "${initial_viral_load_val}"
    --start_tick                 "${start_tick_val}"
    --num_ticks                  "${num_ticks_val}"
)

echo "Executing command: ${command[*]}"
"${command[@]}"

# Verify the FASTA and metadata agree. A mismatch means the run was truncated
# (e.g. killed by the time limit) even though earlier steps looked fine.
fasta_out="${output_prefix}.sequences.fasta.xz"
meta_out="${output_prefix}.metadata.tsv.xz"
if [[ -s "${fasta_out}" && -s "${meta_out}" ]]; then
    n_fasta=$(xz -dc "${fasta_out}" | grep -c '^>' || true)
    n_meta=$(( $(xz -dc "${meta_out}" | wc -l) - 1 ))
    echo "Output check: ${n_fasta} FASTA records, ${n_meta} metadata rows"
    [[ "${n_fasta}" -eq "${n_meta}" ]] \
        || echo "WARNING: FASTA/metadata row count mismatch - output may be truncated." >&2
fi

# To disable rate limiting, you could set:
# enable_rate_limit=""
# And then ensure initial_viral_load, start_tick, num_ticks are not passed or handled by the script if enable_rate_limit is empty.
# However, the python script's argparse will handle defaults if these are not provided and --rate_limit is not set.

## Notes on commented out variables from original script:
## --align_fasta: Not needed if --seed_fasta is provided for generate_sequence_analysis.
## --proportional was directly added to the command string.
## --poor: Not set, so it will use the script's default (False). If you wanted to explicitly set it, you'd add --poor or --poor False.
## --limit: Removed as requested.

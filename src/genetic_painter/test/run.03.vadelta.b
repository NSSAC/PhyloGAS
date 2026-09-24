#!/bin/bash

code="../src/genetic_painter.py"

## Inputs for both analyses.

### Analysis type.
analysis_type="generate_sequence_analysis"

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
output_prefix="/project/bii_nssac/biocomplexity/vdh_genomics/synthetic_biosurveillance/SARS-Cov2-Biosurveillance-Simulation/test/run.03.vadelta.output/run_03_vadelta_2026_09_10_128to428" # hold.filehanged to distinguish from previous runs

# Compress to xz format
compression_type="xz"

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


## Execute.
command="python ${code}                                           \
    --analysis_type          ${analysis_type}            \
    --random_number_seed     ${random_number_seed}       \
    --threshold_file         \"${threshold_file}\"           \
    --base_threshold_df      \"${base_threshold_df}\"        \
    --start_date             ${start_date}               \
    --input_graph_csv        \"${input_graph_csv}\"          \
    --seed_fasta             \"${seed_fasta}\"              \
    --output_prefix          \"${output_prefix}\"            \
    --compression            ${compression_type}         \
    --persontrait_file       \"${persontrait_file}\"         \
    --add_metadata           \"${add_metadata}\"             \
    --input_graph_painted_prefix \"${input_graph_painted_prefix}\" \
    --proportional                                       \
    --reference              \"${reference}\"                \
    --location               '${location_val}'             \
    --reference_location     '${reference_location_val}'     \
    ${enable_rate_limit}                                 \
    --initial_viral_load     ${initial_viral_load_val}    \
    --start_tick             ${start_tick_val}            \
    --num_ticks              ${num_ticks_val}"

echo "Executing command: ${command}"
eval "${command}" # Using eval to correctly interpret quotes within the command string

# To disable rate limiting, you could set:
# enable_rate_limit=""
# And then ensure initial_viral_load, start_tick, num_ticks are not passed or handled by the script if enable_rate_limit is empty.
# However, the python script's argparse will handle defaults if these are not provided and --rate_limit is not set.

## Notes on commented out variables from original script:
## --align_fasta: Not needed if --seed_fasta is provided for generate_sequence_analysis.
## --proportional was directly added to the command string.
## --poor: Not set, so it will use the script's default (False). If you wanted to explicitly set it, you'd add --poor or --poor False.
## --limit: Removed as requested.

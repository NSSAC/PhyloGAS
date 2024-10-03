
code="../src/freq_v02.py"

## Inputs for both analyses.

### Analysis type.
analysis_type="generate_sequence_analysis"

random_number_seed=43

threshold_file="run.03.delta/run.03.threshold.file"
base_threshold_df="run.03.delta/run.03.base.threshold.df"

# Whether to write all threshold DFs into one file or
# each in a differnt file.
threshold_df_num_files="all_in_one_file"

# Fasta file of genomic sequences.
# align_fasta="va_variant_BA.2.12.1_5_sequences.fasta"
align_fasta="/project/biocomplexity/vdh_genomics/synthetic_biosurveillance/SARS-Cov2-Biosurveillance-Simulation/data/training_sequences/delta/delta.fa"

## Inputs for generate sequence analysis.

# start_date="2022-06-29"
start_date="2021-06-01"

input_graph_csv="epihiper_exp7_dendrogram.csv"

# Output file prefix.
output_prefix="run_03_with_delta_prop_ref_1000"

# Compress to xz format
compression_type="xz"

# Persontrait file
persontrait_file="/project/biocomplexity/vdh_genomics/synthetic_biosurveillance/SARS-Cov2-Biosurveillance-Simulation/data/merged_population_files/va_merged_person.csv"

# Add metadata to include from the persontrait file
# pid,hid,age,age_group,gender,county_fips,home_latitude,home_longitude,admin1,admin2,admin3,admin4,smh_race,latino,race
add_metadata="gender,county,home_latitude,home_longitude,latino,race,smh_race,age_group,pid"

# In code, this parameter's default value is False.
# This is used in the 'better' mutational model.

# If true, then use the poor mutational model; otherwise
# use the better model.
# Default value is false.
# poor="False"

# Max number of values to process.
# limit="16521"
limit="1000"

# Add reference sequence to the output.
reference=/project/biocomplexity/vdh_genomics/synthetic_biosurveillance/SARS-Cov2-Biosurveillance-Simulation/data/training_sequences/reference.fasta


## Execute.
python ${code}                                           \
    --analysis_type          ${analysis_type}            \
    --random_number_seed     ${random_number_seed}       \
    --threshold_file         ${threshold_file}           \
    --base_threshold_df      ${base_threshold_df}        \
    --threshold_df_num_files ${threshold_df_num_files}   \
    --align_fasta            ${align_fasta}              \
    --start_date             ${start_date}               \
    --input_graph_csv        ${input_graph_csv}          \
    --output_prefix          ${output_prefix}            \
    --compression            ${compression_type}         \
    --persontrait_file       ${persontrait_file}         \
    --add_metadata           ${add_metadata}             \
    --proportional                                       \
    --limit                  ${limit}                    \
    --reference              ${reference}

##     --proportional       ${proportional}        \
##     --poor               ${poor}                \
##     --reference          ${reference}

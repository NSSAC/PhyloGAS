
code="../src/freq_v02.py"

## Inputs for both analyses.

### Analysis type.
analysis_type="generate_sequence_analysis"


threshold_file="run.01.threshold.file"
base_threshold_df="run.01.base.threshold.df"


## Inputs for generate sequence analysis.

start_date="2022-06-29"

input_graph_csv="epihiper_exp7_dendrogram.csv"

# Output file prefix.
output_prefix="run_01_prefix"

# In code, this parameter's default value is False.
# This is used in the 'better' mutational model.
proportional=False

# If true, then use the poor mutational model; otherwise
# use the better model.
# Default value is false.
poor=False

# Max number of values to process.
limit="16521"

# Add reference sequence to the output.
reference=None


## Execute.
python ${code}                                  \
    --analysis_type      ${analysis_type}       \
    --threshold_file     ${threshold_file}      \
    --base_threshold_df  ${base_threshold_df}   \
    --start_date         ${start_date}          \
    --input_graph_csv    ${input_graph_csv}     \
    --output_prefix      ${output_prefix}       \
    --proportional       ${proportional}        \
    --poor               ${poor}                \
    --limit              ${limit}               \
    --reference          ${reference}



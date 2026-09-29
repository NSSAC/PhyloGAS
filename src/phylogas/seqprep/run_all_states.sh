#!/bin/bash
#
# run_all_states.sh
#
# This script runs seq_prep.py in bulk mode to download sequences
# for a specific Pango lineage across all US states and territories.

# --- Configuration ---
PYTHON_SCRIPT="seq_prep.py"
PANGO_LINEAGE="JN.1"
OUTPUT_BASE_DIR="./jn1_sequences_by_state" # A directory to store all the output files

# Optional: Specify a date range. Leave commented out to get all dates.
# DATE_FROM="2023-10-01"
# DATE_TO="2024-03-31"

# Array of US states and major territories
# Note: Cov-Spectrum's 'division' field may not cover all territories.
# This list is comprehensive; some may return no results.
STATES=(
    "Alabama" "Alaska" "Arizona" "Arkansas" "California" "Colorado"
    "Connecticut" "Delaware" "Florida" "Georgia" "Hawaii" "Idaho"
    "Illinois" "Indiana" "Iowa" "Kansas" "Kentucky" "Louisiana"
    "Maine" "Maryland" "Massachusetts" "Michigan" "Minnesota"
    "Mississippi" "Missouri" "Montana" "Nebraska" "Nevada"
    "New Hampshire" "New Jersey" "New Mexico" "New York"
    "North Carolina" "North Dakota" "Ohio" "Oklahoma" "Oregon"
    "Pennsylvania" "Rhode Island" "South Carolina" "South Dakota"
    "Tennessee" "Texas" "Utah" "Vermont" "Virginia" "Washington"
    "West Virginia" "Wisconsin" "Wyoming"
    "District of Columbia" "Puerto Rico" "Guam" "Virgin Islands" "American Samoa"
)

# --- Execution Logic ---

# Create the main output directory
mkdir -p "$OUTPUT_BASE_DIR"

# Construct date arguments if they are set
DATE_ARGS=""
if [[ -n "$DATE_FROM" && -n "$DATE_TO" ]]; then
    DATE_ARGS="--date_from ${DATE_FROM} --date_to ${DATE_TO}"
    echo "Using date range: $DATE_FROM to $DATE_TO"
else
    echo "No date range specified. Querying for all dates."
fi

# Loop through each state in the array
for STATE in "${STATES[@]}"; do
    echo "" # Add a blank line for readability
    echo "============================================================"
    echo "Processing State: ${STATE}"
    echo "============================================================"
    
    # Construct the full command in a bash array for safety
    # Note the use of `eval` to correctly handle the optional DATE_ARGS string
    cmd=(
        "python"
        "${PYTHON_SCRIPT}"
        "--state"
        "${STATE}"
        "--pango"
        "${PANGO_LINEAGE}"
        "--output_folder"
        "${OUTPUT_BASE_DIR}"
    )

    # --- Execute the command ---
    # The 'eval' is used here to correctly parse the optional date arguments string.
    # We pass the main command array and then the date arguments string.
    eval "${cmd[@]}" "$DATE_ARGS"
done

echo ""
echo "All states processed."

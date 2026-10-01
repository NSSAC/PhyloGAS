#!/usr/bin/env bash
# ==============================================================================
# Script: download_pgcoe_dataverse-v2.sh
# Description: Downloads and prepares synthetic population & EpiHiper files
#              from UVA Dataverse for PGCOE state simulations.
# Usage: ./download_pgcoe_dataverse-v2.sh [state_code ...]
# Example: ./download_pgcoe_dataverse-v2.sh va ca ga ma mn wa
# ==============================================================================

set -euo pipefail

# Default to all 6 PGCOE states if none specified
DEFAULT_STATES=("va" "ca" "ga" "ma" "mn" "wa")
STATES=("${@:-${DEFAULT_STATES[@]}}")

BASE_DIR="./data"

# Direct Dataverse File API ID Mapping
declare -A FILE_IDS=(
    # Virginia (VA)
    ["va_persontrait_epihiper.txt.xz"]="121039"
    ["va_contact_network_epihiper.txt.xz"]="121037"
    ["va_household.csv.xz"]="120566"
    ["va_person.csv.xz"]="120565"

    # California (CA)
    ["ca_persontrait_epihiper.txt.xz"]="121117"
    ["ca_contact_network_epihiper.txt.xz"]="121116"
    ["ca_household.csv.xz"]="120623"
    ["ca_person.csv.xz"]="120625"

    # Georgia (GA)
    ["ga_persontrait_epihiper.txt.xz"]="121106"
    ["ga_contact_network_epihiper.txt.xz"]="121107"
    ["ga_household.csv.xz"]="120595"
    ["ga_person.csv.xz"]="120598"

    # Massachusetts (MA)
    ["ma_persontrait_epihiper.txt.xz"]="121110"
    ["ma_contact_network_epihiper.txt.xz"]="121109"
    ["ma_household.csv.xz"]="120591"
    ["ma_person.csv.xz"]="120589"

    # Minnesota (MN)
    ["mn_persontrait_epihiper.txt.xz"]="121112"
    ["mn_contact_network_epihiper.txt.xz"]="121114"
    ["mn_household.csv.xz"]="120581"
    ["mn_person.csv.xz"]="120577"

    # Washington (WA)
    ["wa_persontrait_epihiper.txt.xz"]="121104"
    ["wa_contact_network_epihiper.txt.xz"]="121105"
    ["wa_household.csv.xz"]="120607"
    ["wa_person.csv.xz"]="120604"
)

# Required file keys per state
GET_STATE_FILES() {
    local st="$1"
    echo "${st}_persontrait_epihiper.txt.xz ${st}_contact_network_epihiper.txt.xz ${st}_household.csv.xz ${st}_person.csv.xz"
}

DATAVERSE_API="https://dataverse.lib.virginia.edu/api/access/datafile"

echo "=========================================================="
echo " Starting Dataverse Download for PGCOE States: ${STATES[*]}"
echo " Target Directory: ${BASE_DIR}"
echo "=========================================================="

for state in "${STATES[@]}"; do
    state_lc=$(echo "$state" | tr '[:upper:]' '[:lower:]')
    state_dir="${BASE_DIR}/${state_lc}"
    mkdir -p "$state_dir"

    echo -e "\n---> Processing State: ${state_lc^^}"
    files=$(GET_STATE_FILES "$state_lc")

    for file_name in $files; do
        file_id="${FILE_IDS[$file_name]:-}"
        if [[ -z "$file_id" ]]; then
            echo "Warning: No File ID registered for $file_name. Skipping."
            continue
        fi

        target_file="${state_dir}/${file_name}"
        decompressed_file="${target_file%.xz}"

        # Download if compressed/decompressed file does not exist
        if [[ ! -f "$decompressed_file" && ! -f "$target_file" ]]; then
            echo "Downloading ${file_name} (ID: ${file_id})..."
            curl -L -f --retry 3 --retry-delay 5 \
                "${DATAVERSE_API}/${file_id}" \
                -o "$target_file"
        elif [[ -f "$decompressed_file" ]]; then
            echo "Decompressed file ${decompressed_file} already exists. Skipping download."
            continue
        else
            echo "Compressed file ${target_file} already exists. Skipping download."
        fi

        # Decompress .xz files required by simulation tools
        if [[ -f "$target_file" && "$target_file" == *.xz ]]; then
            echo "Decompressing ${file_name}..."
            unxz -f -k "$target_file"
        fi
    done
done

echo -e "\n=========================================================="
echo " All downloads and preparations completed successfully!"
echo " Data stored in: ${BASE_DIR}/"
echo "=========================================================="

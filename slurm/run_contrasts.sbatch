#!/bin/bash
#SBATCH -p bii                   # Partition name
#SBATCH -A bii_nssac             # Account
#SBATCH -c 8                     # CPU cores per state task (no Nextstrain needed)
#SBATCH -t 12:00:00              # Time limit (hh:mm:ss)
#SBATCH -J phylogas_contrasts    # Job name
#SBATCH --array=0-4              # 5 states: 0=va, 1=ga, 2=ma, 3=mn, 4=wa
#SBATCH -o /dev/null             # Output redirected per state into slurm/logs/ below
#SBATCH -e /dev/null
#SBATCH --mem=64G

# Fire-and-forget SLURM array for the PhyloGAS contrast & mechanism benchmarks
# across all 5 states (va, ga, ma, mn, wa).
#
# Runs:
#   1) Contrast sampling arms in isolated sibling directories:
#      - 03_sampled_foil_s50   (Stratified + LASSO50-Stratified, max_supergroups=50)
#      - 03_sampled_nocounty   (Stratified + LASSO50-Greedy, --stratifiers age race sex)
#      - 03_sampled_countyonly (Stratified + LASSO50-Greedy, --stratifiers county)
#   2) Homophily decay, bridge enrichment, LASSO temporal weights & mechanism summary
#      -> results/<state>_delta_wave/05_benchmarks/mechanism/
#   3) Combined Mugration & ABM Ground-Truth benchmarks across all 44 strategies
#      -> results/<state>_delta_wave/05_benchmarks/Mugration_Metrics_contrasts.csv
#      -> results/<state>_delta_wave/05_benchmarks/AUC_truth_rankings_contrasts.csv
#
# Usage:
#   sbatch slurm/run_contrasts.sbatch                  # all 5 states (array 0-4)
#   sbatch --array=0 slurm/run_contrasts.sbatch        # Virginia only (index 0)

cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
[ -f slurm/prepare_inputs.sh ] || cd ..
[ -f slurm/prepare_inputs.sh ] || {
    echo "submit from the PhyloGAS checkout (or its slurm/ folder); job ran in $(pwd)" \
        > "${SLURM_SUBMIT_DIR:-.}/PHYLOGAS_CONTRASTS_${SLURM_JOB_ID:-local}_FAILED.txt"
    exit 1
}

STATES=(va ga ma mn wa)
if [ -z "${STATE:-}" ]; then
    IDX="${SLURM_ARRAY_TASK_ID:-0}"
    STATE="${STATES[$IDX]}"
fi

JOB_TAG="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}"
mkdir -p slurm/logs
exec >"slurm/logs/contrasts_${STATE}_${JOB_TAG}.out" 2>"slurm/logs/contrasts_${STATE}_${JOB_TAG}.err"

set -eo pipefail

source /project/biocomplexity/asw3xp/miniconda3/bin/activate
conda activate phylogas_env

set -u

trap 'echo "=== seff ${SLURM_JOB_ID:-} ==="; seff "${SLURM_JOB_ID:-}" 2>/dev/null || true' EXIT

CONFIG="configs/${STATE}.yaml"
echo "state: ${STATE}   config: ${CONFIG}   host: $(hostname)   phylogas $(phylogas --version 2>/dev/null)"

phylogas benchmark contrasts --config "${CONFIG}" "$@"

#!/bin/bash
#SBATCH -p bii                   # Partition name
#SBATCH -A bii_nssac             # Account
#SBATCH -c 8                     # CPU cores per state task (no Nextstrain needed)
#SBATCH -t 06:00:00              # Time limit (hh:mm:ss)
#SBATCH -J phylogas_canonical    # Job name
#SBATCH --array=0-4              # 5 states: 0=va, 1=ga, 2=ma, 3=mn, 4=wa
#SBATCH -o /dev/null             # Output redirected per state into slurm/logs/ below
#SBATCH -e /dev/null
#SBATCH --mem=64G

# Fire-and-forget SLURM array for the Canonical Public-Health Sampling Ladder
# across all 5 states (va, ga, ma, mn, wa) — no external population model needed!
#
# Strategies generated in `results/<state>_delta_wave/03_sampled_canonical/`:
#   - convenience         (SURS on K=2 academic hubs + inpatient/outbreak skew)
#   - sentinel            (SURS on K=6 regional sentinel lab hubs, Oltean et al. 2023)
#   - sentinel_caseprop   (Case-Proportional / Subsampler on sentinel swabs using statewide cases)
#   - sentinel_geocap     (Nextstrain `augur --group-by county` equal cap on sentinel swabs)
#   - statewide_surs      (Mandatory statewide submission + uniform random baseline)
#   - statewide_caseprop  (Mandatory statewide submission + county case-proportional quota)
#   - statewide_geocap    (Mandatory statewide submission + equal county cap)
#
# Outputs:
#   - results/<state>_delta_wave/05_benchmarks/Mugration_Metrics_canonical.csv
#   - results/<state>_delta_wave/05_benchmarks/AUC_truth_rankings_canonical.csv
#   - results/<state>_delta_wave/05_benchmarks/canonical_pools/<state>_sentinel_hubs.json
#   - results/<state>_delta_wave/05_benchmarks/canonical_pools/<state>_canonical_pool_summary.csv
#   - results/<state>_delta_wave/05_benchmarks/canonical_pools/<state>_canonical_strategy_mechanism_summary.csv
#
# Usage:
#   sbatch slurm/run_canonical.sbatch                  # all 5 states (array 0-4)
#   sbatch --array=4 slurm/run_canonical.sbatch        # Washington only (index 4)

cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
[ -f slurm/prepare_inputs.sh ] || cd ..
[ -f slurm/prepare_inputs.sh ] || {
    echo "submit from the PhyloGAS checkout (or its slurm/ folder); job ran in $(pwd)" \
        > "${SLURM_SUBMIT_DIR:-.}/PHYLOGAS_CANONICAL_${SLURM_JOB_ID:-local}_FAILED.txt"
    exit 1
}

STATES=(va ga ma mn wa)
if [ -z "${STATE:-}" ]; then
    IDX="${SLURM_ARRAY_TASK_ID:-0}"
    STATE="${STATES[$IDX]}"
fi

JOB_TAG="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}"
mkdir -p slurm/logs
exec >"slurm/logs/canonical_${STATE}_${JOB_TAG}.out" 2>"slurm/logs/canonical_${STATE}_${JOB_TAG}.err"

set -eo pipefail

source /project/biocomplexity/asw3xp/miniconda3/bin/activate
conda activate phylogas_env

set -u

trap 'echo "=== seff ${SLURM_JOB_ID:-} ==="; seff "${SLURM_JOB_ID:-}" 2>/dev/null || true' EXIT

CONFIG="configs/${STATE}.yaml"
echo "state: ${STATE}   config: ${CONFIG}   host: $(hostname)   phylogas $(phylogas --version 2>/dev/null)"

phylogas benchmark canonical --config "${CONFIG}" "$@"

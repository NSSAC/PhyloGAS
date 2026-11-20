#!/bin/bash
#SBATCH -p bii                   # Partition name
#SBATCH -A bii_nssac             # Account
#SBATCH -c 4                     # Single-pass streaming stdlib benchmark (~1-2 min/state)
#SBATCH -t 01:30:00              # Time limit (hh:mm:ss)
#SBATCH -J phylogas_detect       # Job name
#SBATCH --array=0-4              # 5 states: 0=va, 1=ga, 2=ma, 3=mn, 4=wa
#SBATCH -o /dev/null             # Output redirected per state into slurm/logs/ below
#SBATCH -e /dev/null
#SBATCH --mem=32G

# Fast 5-state SLURM array to compute Detection Lead-Time (T95 / UT95 / ΔT95)
# across ALL existing canonical + BeyondBaseline + contrast sample sets WITHOUT
# re-running sampling, homophily mechanism, mugration, or truth AUC.
#
# Writes per state in `results/<state>_delta_wave/05_benchmarks/`:
#   1) Detection_Summary_canonical.csv  & Detection_Events_canonical.csv   (component_id chains >=100)
#   2) Detection_Summary_contrasts.csv  & Detection_Events_contrasts.csv   (all 44 core+contrast arms)
#   3) Detection_Variants_canonical.csv & Detection_Variants_Events_canonical.csv (variant_benchmark sublineages)
#
# Usage:
#   sbatch slurm/run_detection.sbatch                  # all 5 states (array 0-4)
#   sbatch --array=4 slurm/run_detection.sbatch        # Washington only (index 4)

cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
[ -f slurm/prepare_inputs.sh ] || cd ..
[ -f slurm/prepare_inputs.sh ] || {
    echo "submit from the PhyloGAS checkout (or its slurm/ folder); job ran in $(pwd)" \
        > "${SLURM_SUBMIT_DIR:-.}/PHYLOGAS_DETECTION_${SLURM_JOB_ID:-local}_FAILED.txt"
    exit 1
}

STATES=(va ga ma mn wa)
if [ -z "${STATE:-}" ]; then
    IDX="${SLURM_ARRAY_TASK_ID:-0}"
    STATE="${STATES[$IDX]}"
fi

JOB_TAG="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}"
mkdir -p slurm/logs
exec >"slurm/logs/detection_${STATE}_${JOB_TAG}.out" 2>"slurm/logs/detection_${STATE}_${JOB_TAG}.err"

set -eo pipefail

source /project/biocomplexity/asw3xp/miniconda3/bin/activate
conda activate phylogas_env

set -u

trap 'echo "=== seff ${SLURM_JOB_ID:-} ==="; seff "${SLURM_JOB_ID:-}" 2>/dev/null || true' EXIT

CONFIG="configs/${STATE}.yaml"
RDIR="results/${STATE}_delta_wave"
echo "state: ${STATE}   config: ${CONFIG}   host: $(hostname)   phylogas $(phylogas --version 2>/dev/null)"

# 1. Canonical Public-Health Ladder (established importation chains, component_id >= 100)
phylogas benchmark canonical --config "${CONFIG}" \
    --skip-sampling --skip-mechanism --skip-mugration --skip-truth

# 2. Core + Contrast BeyondBaseline Suite (44 arms, component_id >= 100)
phylogas benchmark contrasts --config "${CONFIG}" \
    --skip-sampling --skip-mechanism --skip-mugration --skip-truth

# 3. Variant/Sublineage-level Detection on Canonical + Reference Arms (variant_benchmark >= 100)
START_ARG=()
if [ -f "${RDIR}/03_sampled_datasets/sampling_start_date.txt" ]; then
    START_ARG=(--start-date "$(cat "${RDIR}/03_sampled_datasets/sampling_start_date.txt")")
fi

SAMPLE_ARGS=("${RDIR}/03_sampled_canonical/"*_samples.csv.xz)
for ref in "1S__surs_samples.csv.xz" "4S__surs_samples.csv.xz" "1S-1_LL-P__lasso_greedy_samples.csv.xz"; do
    if [ -f "${RDIR}/03_sampled_datasets/${ref}" ]; then
        SAMPLE_ARGS+=("${RDIR}/03_sampled_datasets/${ref}")
    fi
done

phylogas benchmark detection \
    --samples "${SAMPLE_ARGS[@]}" \
    --infections "${RDIR}/05_benchmarks/allevents_variants.csv.xz" \
    --linelist "${RDIR}/02_simulated_linelists/linelist.csv.xz" \
    --target-col variant_benchmark \
    --min-chain-size 100 \
    "${START_ARG[@]}" \
    --out "${RDIR}/05_benchmarks/Detection_Variants_canonical.csv" \
    --out-events "${RDIR}/05_benchmarks/Detection_Variants_Events_canonical.csv"

echo "=== Completed all detection benchmarks for ${STATE} ==="

#!/bin/bash
#SBATCH -p bii                   # Partition name
#SBATCH -A bii_nssac             # Account
#SBATCH -c 8                     # Cores per outbreak task (JUNIPER subtrees)
#SBATCH -t 05:00:00              # Wall-clock limit (~30-60 min expected per task)
#SBATCH -J phylogas_outbreak     # Job name
#SBATCH --array=0-29             # 30 tasks = 5 states x 6 outbreaks/state (3 clusters x 2 pi)
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mem=32G

# Stage 2 (30-task wide SLURM array across all 5 states):
#   Runs `outbreaker2`, `TransPhylo` (on TreeTime dated tree), and `JUNIPER` (`juniper0`)
#   on one prepared focal outbreak directory (`N ~ 40..180` cases), then scores all
#   5 methods (`MicrobeTrace_1SNP`, `Nextstrain`, `outbreaker2`, `TransPhylo`, `JUNIPER`)
#   and updates `results/<state>_delta_wave/05_benchmarks/Outbreak_Reconstruction_Metrics.csv`.
#
# Usage:
#   JID=$(sbatch --parsable slurm/run_phylo_bench.sbatch)
#   sbatch --dependency=afterok:${JID} slurm/run_outbreak_bench.sbatch

cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
[ -f slurm/prepare_inputs.sh ] || cd ..
[ -f slurm/prepare_inputs.sh ] || {
    echo "submit from the PhyloGAS checkout (or its slurm/ folder); job ran in $(pwd)" \
        > "${SLURM_SUBMIT_DIR:-.}/PHYLOGAS_OB_${SLURM_JOB_ID:-local}_FAILED.txt"
    exit 1
}

REPO_ROOT="$(pwd)"
STATES=(va ga ma mn wa)
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
STATE_IDX=$(( TASK_ID / 6 ))
OB_IDX=$(( TASK_ID % 6 ))
STATE="${STATES[$STATE_IDX]}"

JOB_TAG="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}"
mkdir -p slurm/logs
exec >"slurm/logs/outbreak_${STATE}_idx${OB_IDX}_${JOB_TAG}.out" 2>"slurm/logs/outbreak_${STATE}_idx${OB_IDX}_${JOB_TAG}.err"

set -eo pipefail

CONFIG="configs/${STATE}.yaml"
RDIR="results/${STATE}_delta_wave"
INDEX_JSON="${RDIR}/07_outbreak_bench/outbreak_index_${STATE}.json"

if [ ! -f "${INDEX_JSON}" ]; then
    echo "ERROR: ${INDEX_JSON} not found. Run slurm/run_phylo_bench.sbatch first." >&2
    exit 1
fi

# Resolve this array task's outbreak directory via stdlib python3 before loading R
WORK_DIR="$(python3 -c '
import json, sys
data = json.load(open(sys.argv[1]))
idx = int(sys.argv[2])
if idx >= len(data):
    sys.exit(f"Index {idx} out of range for {len(data)} outbreaks")
print(data[idx]["run_dir"])
' "${INDEX_JSON}" "${OB_IDX}")"

echo "task: ${TASK_ID}   state: ${STATE}   ob_idx: ${OB_IDX}   workdir: ${WORK_DIR}   host: $(hostname)"

# ---- Step 0: Build (or refresh) TreeTime dated tree in decimal years ----
if [ ! -f "${WORK_DIR}/treetime/tree.nwk" ] || [ ! -f "${WORK_DIR}/treetime/branch_lengths.json" ]; then
    echo "=== [Step 0] Building Nextstrain (IQ-TREE + TreeTime) dated tree for ${WORK_DIR} ==="
    module purge
    if [ -d "/project/biocomplexity/${USER:-asw3xp}/nextstrain" ]; then
        export NEXTSTRAIN_HOME="${NEXTSTRAIN_HOME:-/project/biocomplexity/${USER:-asw3xp}/nextstrain}"
    fi
    source /project/biocomplexity/asw3xp/miniconda3/bin/activate
    conda activate phylogas_env
    unset PYTHONPATH
    phylogas benchmark outbreaks --config "${CONFIG}" --step nextstrain \
        --workdir "${WORK_DIR}" --threads "${SLURM_CPUS_PER_TASK:-8}"
    conda deactivate 2>/dev/null || true
else
    echo "=== [Step 0] Refreshing decimal-year dated tree from existing TreeTime output ==="
    module purge
    source /project/biocomplexity/asw3xp/miniconda3/bin/activate
    conda activate phylogas_env
    unset PYTHONPATH
    phylogas benchmark outbreaks --config "${CONFIG}" --step refresh-dated-tree \
        --workdir "${WORK_DIR}"
    conda deactivate 2>/dev/null || true
fi

# ---- Step 1: Run R Bayesian Outbreak Tools (outbreaker2, TransPhylo, JUNIPER) ----
module purge
module load gcc/14.2.0 openmpi/5.0.7 R/4.6.0

USER_RLIB="${HOME}/Rlibs/phylogas_R4.6"
SHARED_RLIB="/project/bii_nssac/people/bl4zc/Rlibs/juniper_R4.6"
if [ -d "${SHARED_RLIB}" ] && [ -r "${SHARED_RLIB}" ]; then
    export R_LIBS_USER="${USER_RLIB}:${SHARED_RLIB}"
else
    export R_LIBS_USER="${USER_RLIB}"
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

JUN_GLOBAL="${JUNIPER_GLOBAL:-1500}"
OB_TOOLS="${OUTBREAK_TOOLS:-outbreaker2,transphylo,juniper}"
Rscript "${REPO_ROOT}/scripts/run_outbreak_tools.R" \
    --workdir "${WORK_DIR}" \
    --tools "${OB_TOOLS}" \
    --juniper-global "${JUN_GLOBAL}"

# ---- Step 2: Score all 5 methods & safely aggregate per-state CSV ----
# Unload R modules & clear PYTHONPATH so R's scipy-bundle never shadows conda numpy/phylogas
module purge
unset PYTHONPATH

source /project/biocomplexity/asw3xp/miniconda3/bin/activate
conda activate phylogas_env
unset PYTHONPATH

set -u
trap 'echo "=== seff ${SLURM_JOB_ID:-} ==="; seff "${SLURM_JOB_ID:-}" 2>/dev/null || true' EXIT

phylogas benchmark outbreaks --config "${CONFIG}" --step score --workdir "${WORK_DIR}"

(
    flock -x 9
    phylogas benchmark outbreaks --config "${CONFIG}" --step aggregate
) 9>"${RDIR}/07_outbreak_bench/.aggregate.lock"

echo "=== Completed outbreak reconstruction benchmark for ${WORK_DIR} ==="

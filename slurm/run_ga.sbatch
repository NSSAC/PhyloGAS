#!/bin/bash
#SBATCH -p bii                   # Partition name
#SBATCH -A bii_nssac             # Account
#SBATCH -c 32                    # Number of CPU cores
#SBATCH -t 48:00:00              # Time limit (hh:mm:ss)
#SBATCH -J gadelta_job           # Job name
#SBATCH -o /dev/null              # Output goes to slurm/logs/ -- redirected below,
#SBATCH -e /dev/null              # since -o paths depend on where you submit from
#SBATCH --mem=128G
### SBATCH --mem=1TB                # Memory request; some systems may require --mem=1000G

# PhyloGAS end to end for Georgia: stage inputs, then training -> painting ->
# line list -> sampling -> Nextstrain -> benchmarks.
#
#   sbatch slurm/run_ga.sbatch     (from the PhyloGAS checkout, or from slurm/)

# sbatch COPIES this file into a spool directory before executing it, so
# ${BASH_SOURCE[0]} points at the spool copy, not at the repo. SLURM_SUBMIT_DIR
# is the reliable anchor; the BASH_SOURCE fallback only applies when this
# script is run directly (outside SLURM).
cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# Accept submission from the PhyloGAS checkout or from its slurm/ folder.
[ -f slurm/prepare_inputs.sh ] || cd ..
[ -f slurm/prepare_inputs.sh ] || {  # logs/ is unknown here, so note it where submitted
    echo "submit from the PhyloGAS checkout (or its slurm/ folder); job ran in $(pwd)" \
        > "${SLURM_SUBMIT_DIR:-.}/PHYLOGAS_JOB_${SLURM_JOB_ID:-local}_FAILED.txt"; exit 1; }

# All output to slurm/logs/, whichever folder the job was submitted from.
# Done before anything else so conda errors are captured too. SLURM's own
# messages (time limit, out-of-memory kill) do not pass through here: check
# those with  sacct -j <jobid>  or  seff <jobid>.
mkdir -p slurm/logs
exec >"slurm/logs/gadelta_${SLURM_JOB_ID:-local}.out" 2>"slurm/logs/gadelta_${SLURM_JOB_ID:-local}.err"

# Fail fast so a broken step cannot silently produce a truncated dataset.
# NOTE: `-u` (error on unset variable) is deliberately NOT enabled yet. Conda's
# activation hooks read variables before assigning them -- e.g.
#   etc/conda/activate.d/gdal-activate.sh:  export _CONDA_SET_GDAL_DATA=$GDAL_DATA
# which aborts the job under `-u` with "GDAL_DATA: unbound variable".
# Those hooks are third-party, so we activate first and tighten afterwards.
set -eo pipefail

source /project/biocomplexity/asw3xp/miniconda3/bin/activate
conda activate phylogas_env

# Now that conda's hooks have run, catch our own typos.
set -u


# Report the resources actually consumed so the next submission can be tuned.
trap 'echo "=== seff ${SLURM_JOB_ID:-} ==="; seff "${SLURM_JOB_ID:-}" 2>/dev/null || true' EXIT

STATE=ga
CONFIG=configs/${STATE}.yaml
echo "state: ${STATE}   config: ${CONFIG}   host: $(hostname)   phylogas $(phylogas --version 2>/dev/null)"

# Inputs (downloads, seeds) -- skipped where already present and correct.
bash slurm/prepare_inputs.sh "${STATE}"

# `full` = the default target plus the ground-truth benchmarks.
# --rerun-incomplete: a job killed mid-write (time limit) leaves partial
# outputs that would otherwise block the resubmission.
phylogas run --config "${CONFIG}" --cores "${SLURM_CPUS_PER_TASK:-8}" full \
    --rerun-incomplete --keep-going

#!/bin/bash
#SBATCH -p bii          # Partition name
#SBATCH -A bii_nssac             # Account
#SBATCH -c 30                    # Number of CPU cores
#SBATCH -t 24:00:00              # Time limit (hh:mm:ss)
#SBATCH -J vadelta_job           # (Optional) Job name
#SBATCH -o vadelta_%j.out        # (Optional) Standard output file (%j expands to job ID)
#SBATCH -e vadelta_%j.err        # (Optional) Standard error file
### SBATCH --mem=1TB                # Memory request; some systems may require --mem=1000G

# Fail fast so a broken step cannot silently produce a truncated dataset.
# NOTE: `-u` (error on unset variable) is deliberately NOT enabled yet. Conda's
# activation hooks read variables before assigning them -- e.g.
#   etc/conda/activate.d/gdal-activate.sh:  export _CONDA_SET_GDAL_DATA=$GDAL_DATA
# which aborts the job under `-u` with "GDAL_DATA: unbound variable".
# Those hooks are third-party, so we activate first and tighten afterwards.
set -eo pipefail

# Run your command/script
source /project/biocomplexity/asw3xp/miniconda3/bin/activate

# Now that conda's hooks have run, catch our own typos.
set -u

# Report the resources actually consumed so the next submission can be tuned.
trap 'echo "=== seff ${SLURM_JOB_ID:-} ==="; seff "${SLURM_JOB_ID:-}" 2>/dev/null || true' EXIT

#bash run.03.vadelta.b
bash run.03.vadelta.both

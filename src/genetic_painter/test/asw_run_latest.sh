#!/bin/bash
#SBATCH -p bii          # Partition name
#SBATCH -A bii_nssac             # Account
#SBATCH -c 30                    # Number of CPU cores
#SBATCH -t 24:00:00              # Time limit (hh:mm:ss)
#SBATCH -J vadelta_job           # (Optional) Job name
#SBATCH -o vadelta_%j.out        # (Optional) Standard output file (%j expands to job ID)
#SBATCH -e vadelta_%j.err        # (Optional) Standard error file
### SBATCH --mem=1TB                # Memory request; some systems may require --mem=1000G

set -euo pipefail                # Fail fast: an unset var or failed step should notsilently produce a truncated dataset.
# Run your command/script
source /project/biocomplexity/asw3xp/miniconda3/bin/activate

# Report the resources actually consumed so the next submission can be tuned.
trap 'echo "=== seff ${SLURM_JOB_ID:-} ==="; seff "${SLURM_JOB_ID:-}" 2>/dev/null || true' EXIT

#bash run.03.vadelta.b
bash run.03.vadelta.both

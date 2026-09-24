#!/bin/bash
#SBATCH -p bii-largemem          # Partition name
#SBATCH -A bii_nssac             # Account
#SBATCH -c 30                    # Number of CPU cores
#SBATCH --mem=1TB                # Memory request; some systems may require --mem=1000G
#SBATCH -t 24:00:00              # Time limit (hh:mm:ss)
#SBATCH -J vadelta_job           # (Optional) Job name
#SBATCH -o vadelta_%j.out        # (Optional) Standard output file (%j expands to job ID)
#SBATCH -e vadelta_%j.err        # (Optional) Standard error file

# Run your command/script
source /project/biocomplexity/asw3xp/miniconda3/bin/activate
#bash run.03.vadelta.b
bash run.03.vadelta.both

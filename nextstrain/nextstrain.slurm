#!/bin/bash
#SBATCH --partition=bii
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=62:00:00
#SBATCH --account=bii_nssac
#SBATCH --mem=0

module purge
module load singularity/3.7.1
#module load bioconda/py3.8
#module load anaconda/2020.11-py3.8
#source activate /project/biocomplexity/dbg2pu/python38-env  
cd /project/biocomplexity/dbg2pu/nextstrain/ncov
singularity /project/biocomplexity/asw3xp/nextstrain.sif exec nextstrain check-setup --set-default && nextstrain build . --configfile $1

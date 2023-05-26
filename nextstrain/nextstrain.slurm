#!/bin/bash
#SBATCH --partition=bii
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=72:00:00
#SBATCH --account=nssac_students
#SBATCH --mem=0

module purge
module load bioconda/py3.8
module load anaconda/2020.11-py3.8
source activate /project/biocomplexity/dbg2pu/python38-env  

nextstrain build . --configfile $1

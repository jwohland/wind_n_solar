#!/bin/bash
#SBATCH --q overrun
#SBATCH --time-min=02:00:00
#SBATCH --ntasks=1
#SBATCH --job-name=retrieval
#SBATCH --output=retrieval.txt

srun bash retrieve_20CRv3_NERSC.sh
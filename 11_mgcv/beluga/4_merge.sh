#!/bin/bash
#SBATCH --time=00:50:00
#SBATCH --mem=7G
#SBATCH --account=def-claprise
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err




echo 'Start 4_merge all CpGs for AA!'

date  # Print the current system time after starting

source /home/laj773/Venv/MethEnv/bin/activate
module load r/4.4.0

export R_LIBS_USER=~/scratch/R/library

Rscript  4_merge.R 

echo 'Finished 4_merge all CpGs for AA!'

date  # Print the current system time after starting

sleep 30
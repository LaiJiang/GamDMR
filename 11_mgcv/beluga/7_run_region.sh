#!/bin/bash
#SBATCH --time=02:59:00
#SBATCH --mem=7G
#SBATCH --account=def-claprise
#SBATCH --array=1-10
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err




echo 'Start 7_run_region all regions for AA!'

date  # Print the current system time after starting

source /home/laj773/Venv/MethEnv/bin/activate
module load r/4.4.0

export R_LIBS_USER=~/scratch/R/library

Rscript  7_run_region.R $SLURM_ARRAY_TASK_ID

echo 'Finished 7_run_region all regions for AA!'

date  # Print the current system time after starting

sleep 30
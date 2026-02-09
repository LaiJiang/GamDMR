#!/bin/bash
#SBATCH --time=02:55:00
#SBATCH --mem=7G
#SBATCH --account=def-claprise
#SBATCH --array=1-20
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err




echo 'Start 8_rerun all CpGs for AA!'

date  # Print the current system time after starting

source /home/laj773/Venv/MethEnv/bin/activate
module load r/4.4.0

export R_LIBS_USER=~/scratch/R/library

Rscript  1_test_run.R $SLURM_ARRAY_TASK_ID

echo 'Finished 9_somnibus all CpGs for AA!'

date  # Print the current system time after starting

sleep 30
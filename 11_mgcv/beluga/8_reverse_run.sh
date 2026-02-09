#!/bin/bash
#SBATCH --time=02:59:00
#SBATCH --mem=7G
#SBATCH --account=def-claprise
#SBATCH --array=1-170
#SBATCH --output=logsB2/slurm-%A_%a.out
#SBATCH --error=logsB2/slurm-%A_%a.err




echo 'Start 8_reverse_run all regions for AA!'

date  # Print the current system time after starting

source /home/laj773/Venv/MethEnv/bin/activate
module load r/4.4.0

export R_LIBS_USER=~/scratch/R/library

Rscript  8_reverse_run.R $SLURM_ARRAY_TASK_ID

echo 'Finished 8_reverse_run all regions for AA!'

date  # Print the current system time after starting

sleep 30
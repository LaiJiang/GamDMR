#!/bin/bash
#SBATCH --time=02:59:00
#SBATCH --mem=7G
#SBATCH --account=def-claprise
#SBATCH --array=1-2
#SBATCH --output=logsB3/slurm-%A_%a.out
#SBATCH --error=logsB3/slurm-%A_%a.err




echo 'Start 10_batch3 all regions for AA!'

date  # Print the current system time after starting

source /home/laj773/Venv/MethEnv/bin/activate
module load r/4.4.0

export R_LIBS_USER=~/scratch/R/library

Rscript  10_batch3.R $SLURM_ARRAY_TASK_ID

echo 'Finished 10_batch3 all regions for AA!'

date  # Print the current system time after starting

sleep 30
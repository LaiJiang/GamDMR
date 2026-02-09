#!/bin/bash
#SBATCH --time=02:50:00
#SBATCH --mem=10G
#SBATCH --account=def-claprise
#SBATCH --array=1-500

echo 'Start 6_run all CpGs for AA!'

date  # Print the current system time after starting

source /home/laj773/Venv/MethEnv/bin/activate
module load r/4.4.0

export R_LIBS_USER=~/scratch/R/library

Rscript  6_3_job.R $SLURM_ARRAY_TASK_ID

echo 'Finished 6_run all CpGs for AA!'

date  # Print the current system time after starting

sleep 30
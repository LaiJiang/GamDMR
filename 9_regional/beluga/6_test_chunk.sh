#!/bin/bash
#SBATCH --time=02:55:00
#SBATCH --mem=15G
#SBATCH --account=def-claprise
#SBATCH --array=1-20
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err




echo 'Start 9_regional somnibus for AA!: 6_test_chunk.R'

date  # Print the current system time after starting

source /home/laj773/Venv/MethEnv/bin/activate
module load r/4.4.0

export R_LIBS_USER=~/scratch/R/library

Rscript  6_test_chunk.R $SLURM_ARRAY_TASK_ID

echo 'Finished 6_test_chunk.R in this chunk for AA!'

date  # Print the current system time after starting

sleep 30
#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --mem=24G
#SBATCH --account=def-claprise
#SBATCH --array=1-863
#SBATCH --output=StageIIlogs/slurm-%A_%a.out
#SBATCH --error=StageIIlogs/slurm-%A_%a.err

echo "Start 9_stage2.R for AA with task IDs failed in stage 1"

date

source /home/laj773/Venv/MethEnv/bin/activate
module load r/4.4.0

export R_LIBS_USER=~/scratch/R/library

# Read the real ID from the txt file based on SLURM_ARRAY_TASK_ID (line number)
REAL_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /home/laj773/scratch/UQAC/meth/results/9_regional/8_col_stage1.txt)

echo "Mapped SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID to real data chunk ID $REAL_ID"

Rscript 9_stage2.R $REAL_ID

echo "Finished 6_test_chunk.R for real data chunk ID $REAL_ID"

date

sleep 30

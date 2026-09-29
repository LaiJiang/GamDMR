#!/bin/bash
#SBATCH --account=def-claprise
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-990%100
#SBATCH --output=logs/mgcv_leakfree_%A_%a.out
#SBATCH --error=logs/mgcv_leakfree_%A_%a.err

module load StdEnv/2023
module load gcc/12.3
module load r/4.4.0

SPLIT_ID=${SPLIT_ID}
JOB_ID=${SLURM_ARRAY_TASK_ID}
N_JOBS=990

Rscript \
  01_run_mgcv_leakfree.R \
  ${SPLIT_ID} \
  ${JOB_ID} \
  ${N_JOBS}

#!/bin/bash
#SBATCH --job-name=cpg_train_assoc
#SBATCH --account=def-claprise
#SBATCH --array=1-990%100
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=24:00:00
#SBATCH --output=logs/cpg_train_assoc_%A_%a.out
#SBATCH --error=logs/cpg_train_assoc_%A_%a.err

set -euo pipefail

module purge
module load StdEnv/2023
module load r/4.4.0

# Use the same personal R library that contains glmmLasso/lme4/lmerTest.
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"

# ---------------------------------------------------------------------------
# Required preparation:
# Save the final_table object once before submission, for example in R:
#
# saveRDS(
#   final_table,
#   "~/scratch/UQAC/meth/results/15_revision/5_cv/6_single_cpg/final_table_for_cv.rds"
# )
# ---------------------------------------------------------------------------

export PATH_WK="${PATH_WK:-$HOME/scratch/UQAC/meth/}"

export FINAL_TABLE_FILE="${FINAL_TABLE_FILE:-$PATH_WK/results/15_revision/5_cv/6_single_cpg/final_table_for_cv.rds}"
export TRAIN_SPLIT_FILE="${TRAIN_SPLIT_FILE:-$PATH_WK/scr/11_mgcv/dat/100_family_training_splits.csv}"
export PHENO_FILE="${PHENO_FILE:-$PATH_WK/scr/11_mgcv/dat/18_pheno_BMI.RData}"
export METH_CHUNK_DIR="${METH_CHUNK_DIR:-$PATH_WK/data/meth_split}"
export SINGLE_CPG_ASSOC_OUTPUT="${SINGLE_CPG_ASSOC_OUTPUT:-$PATH_WK/results/15_revision/5_cv/6_single_cpg/1_training_associations}"

export N_JOBS="${N_JOBS:-990}"
export MIN_COMPLETE_N="${MIN_COMPLETE_N:-31}"
export CACHE_CHUNKS="${CACHE_CHUNKS:-3}"
export OVERWRITE="${OVERWRITE:-0}"

mkdir -p logs
mkdir -p "$SINGLE_CPG_ASSOC_OUTPUT"

echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "FINAL_TABLE_FILE=$FINAL_TABLE_FILE"
echo "OUTPUT=$SINGLE_CPG_ASSOC_OUTPUT"

Rscript --vanilla 12_run_single_cpg_training_associations.R "$SLURM_ARRAY_TASK_ID"

echo "Finished: $(date)"

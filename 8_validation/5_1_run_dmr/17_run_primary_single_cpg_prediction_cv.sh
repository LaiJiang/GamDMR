#!/bin/bash
#SBATCH --job-name=cpg_primary_pred
#SBATCH --account=def-claprise
#SBATCH --array=1-100%50
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/cpg_primary_pred_%A_%a.out
#SBATCH --error=logs/cpg_primary_pred_%A_%a.err

set -euo pipefail

module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/scratch/R/library:$HOME/R/x86_64-pc-linux-gnu-library/4.4"

export PATH_WK="${PATH_WK:-$HOME/scratch/UQAC/meth/}"
export SINGLE_CPG_SELECTED_DIR="${SINGLE_CPG_SELECTED_DIR:-$PATH_WK/results/15_revision/5_cv/6_single_cpg/2_collected_significant/selected_by_split}"
export TRAIN_SPLIT_FILE="${TRAIN_SPLIT_FILE:-$PATH_WK/scr/11_mgcv/dat/100_family_training_splits.csv}"
export PHENO_FILE="${PHENO_FILE:-$PATH_WK/scr/11_mgcv/dat/18_pheno_BMI.RData}"
export METH_CHUNK_DIR="${METH_CHUNK_DIR:-$PATH_WK/data/meth_split}"
export SINGLE_CPG_PREDICTION_OUTPUT="${SINGLE_CPG_PREDICTION_OUTPUT:-$PATH_WK/results/15_revision/5_cv/5_prediction/6_single_cpg_primary}"

export INNER_FOLDS="${INNER_FOLDS:-5}"
export SEED_BASE="${SEED_BASE:-20260803}"
export MIN_TRAIN_CPG_OBSERVED_FRACTION="${MIN_TRAIN_CPG_OBSERVED_FRACTION:-0.50}"
export VARIANCE_EPSILON="${VARIANCE_EPSILON:-1e-8}"
export MAX_CPGS_PER_MODEL="${MAX_CPGS_PER_MODEL:-0}"
export OVERWRITE="${OVERWRITE:-0}"

mkdir -p logs
mkdir -p "$SINGLE_CPG_PREDICTION_OUTPUT"

echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Started: $(date)"
echo "Selected CpGs: $SINGLE_CPG_SELECTED_DIR"
echo "Output: $SINGLE_CPG_PREDICTION_OUTPUT"

Rscript --vanilla \
  16_run_primary_single_cpg_prediction_cv.R \
  "$SLURM_ARRAY_TASK_ID"

echo "Finished: $(date)"

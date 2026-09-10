#!/bin/bash
#SBATCH --job-name=cpg_primary_assoc
#SBATCH --account=def-claprise
#SBATCH --array=1-990%100
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=logs/cpg_primary_assoc_%A_%a.out
#SBATCH --error=logs/cpg_primary_assoc_%A_%a.err

set -euo pipefail

module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/scratch/R/library:$HOME/R/x86_64-pc-linux-gnu-library/4.4"

echo "R library paths:"
Rscript -e 'print(.libPaths())'

echo "Checking required packages:"
Rscript -e '
pkgs <- c("data.table", "lme4", "lmerTest", "glmmLasso")
for (p in pkgs) {
  cat(p, ": ")
  if (requireNamespace(p, quietly = TRUE)) {
    cat(find.package(p), "\n")
  } else {
    stop("Required package not found: ", p)
  }
}'

export PATH_WK="${PATH_WK:-$HOME/scratch/UQAC/meth/}"

export PRIMARY_CPG_MANIFEST_ROOT="${PRIMARY_CPG_MANIFEST_ROOT:-$PATH_WK/results/15_revision/5_cv/6_single_cpg/0_primary_all_cpg_manifest}"

export TRAIN_SPLIT_FILE="${TRAIN_SPLIT_FILE:-$PATH_WK/scr/11_mgcv/dat/100_family_training_splits.csv}"
export PHENO_FILE="${PHENO_FILE:-$PATH_WK/scr/11_mgcv/dat/18_pheno_BMI.RData}"
export METH_CHUNK_DIR="${METH_CHUNK_DIR:-$PATH_WK/data/meth_split}"

# EXACTLY the same result directory used by the pilot.
export SINGLE_CPG_ASSOC_OUTPUT="${SINGLE_CPG_ASSOC_OUTPUT:-$PATH_WK/results/15_revision/5_cv/6_single_cpg/1_training_associations}"

export N_JOBS="${N_JOBS:-990}"
export MIN_COMPLETE_N="${MIN_COMPLETE_N:-31}"

# Primary default: all three single-CpG models.
export MODELS="${MODELS:-M1,M2,M3}"

# Default is intentionally 1 because these primary jobs are meant to replace
# the pilot files with the same filenames.
export OVERWRITE="${OVERWRITE:-1}"

export PROGRESS_EVERY="${PROGRESS_EVERY:-25}"
export WRITE_BUFFER_CPGS="${WRITE_BUFFER_CPGS:-10}"

mkdir -p logs
mkdir -p "$SINGLE_CPG_ASSOC_OUTPUT"

echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "MODELS=$MODELS"
echo "PRIMARY_CPG_MANIFEST_ROOT=$PRIMARY_CPG_MANIFEST_ROOT"
echo "OUTPUT=$SINGLE_CPG_ASSOC_OUTPUT"

Rscript --vanilla \
  13_run_primary_single_cpg_training_associations.R \
  "$SLURM_ARRAY_TASK_ID"

echo "Finished: $(date)"

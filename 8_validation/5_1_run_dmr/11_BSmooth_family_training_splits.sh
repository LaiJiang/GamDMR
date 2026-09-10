#!/bin/bash
#SBATCH --time=02:59:00
#SBATCH --mem=4G
#SBATCH --account=def-claprise
#SBATCH --array=1-990%100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=BSmoothTrainLogs/slurm-%A_%a.out
#SBATCH --error=BSmoothTrainLogs/slurm-%A_%a.err

set -euo pipefail

echo "Start BSmooth family-level training-split analysis"
date
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/%p-library/%v}"
mkdir -p "$R_LIBS_USER"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Production defaults; sbatch --export values can override them.
export N_JOBS="${N_JOBS:-990}"
export SPLIT_START="${SPLIT_START:-1}"
export SPLIT_END="${SPLIT_END:-100}"
export MAX_REGIONS="${MAX_REGIONS:-0}"
export MAX_RUNTIME_MINUTES="${MAX_RUNTIME_MINUTES:-165}"
export OVERWRITE="${OVERWRITE:-0}"
export RETRY_ERRORS="${RETRY_ERRORS:-0}"

# BSmooth settings retained from the original analysis.
export BSMOOTH_MIN_CPGS="${BSMOOTH_MIN_CPGS:-10}"
export BSMOOTH_NS="${BSMOOTH_NS:-25}"
export BSMOOTH_H="${BSMOOTH_H:-400}"
export BSMOOTH_MAX_GAP="${BSMOOTH_MAX_GAP:-50000000}"
export BSMOOTH_TSTAT_CUTOFF="${BSMOOTH_TSTAT_CUTOFF:-4.417}"
export BSMOOTH_MIN_GROUP_SAMPLES="${BSMOOTH_MIN_GROUP_SAMPLES:-2}"
export BSMOOTH_ESTIMATE_VAR="${BSMOOTH_ESTIMATE_VAR:-group2}"
export BSMOOTH_LOCAL_CORRECT="${BSMOOTH_LOCAL_CORRECT:-1}"

Rscript --vanilla 10_BSmooth_family_training_splits.R \
  "${SLURM_ARRAY_TASK_ID}"

echo "Finished BSmooth family-level training-split analysis"
date

#!/bin/bash
#SBATCH --time=02:59:00
#SBATCH --mem=7G
#SBATCH --account=def-claprise
#SBATCH --array=1-100%20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=MGCVPredictLogs/slurm-%A_%a.out
#SBATCH --error=MGCVPredictLogs/slurm-%A_%a.err

set -euo pipefail

echo "Start MGCV primary family-split prediction analysis"
date
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

module --force purge
module load StdEnv/2023
module load gcc/12.3
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/%p-library/%v}"
mkdir -p "$R_LIBS_USER"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Primary-analysis settings.
export INNER_FOLDS="${INNER_FOLDS:-5}"
export N_DMR_PCS="${N_DMR_PCS:-2}"
export MAX_DMRS="${MAX_DMRS:-0}"
export MIN_TRAIN_CPG_OBSERVED_FRACTION="${MIN_TRAIN_CPG_OBSERVED_FRACTION:-0.50}"
export MIN_VARIABLE_CPGS="${MIN_VARIABLE_CPGS:-1}"
export VARIANCE_EPSILON="${VARIANCE_EPSILON:-1e-8}"
export SEED_BASE="${SEED_BASE:-20260803}"
export OVERWRITE="${OVERWRITE:-0}"

Rscript --vanilla 1_run_mgcv_primary_prediction_cv.R \
  "${SLURM_ARRAY_TASK_ID}"

echo "Finished MGCV primary family-split prediction analysis"
date

#!/bin/bash
#SBATCH --time=02:50:00
#SBATCH --mem=4G
#SBATCH --account=def-claprise
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=MGCVTestCollectLogs/slurm-%j.out
#SBATCH --error=MGCVTestCollectLogs/slurm-%j.err

set -euo pipefail

echo "Start MGCV all-subject PCA prediction result collection"
date
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"

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

export EXPECTED_SPLITS="${EXPECTED_SPLITS:-100}"
export OVERWRITE="${OVERWRITE:-0}"

# Directory containing performance/, predictions/, coefficients/, feature_qc/,
# and split_qc/ from 10_run_mgcv_test_pca_all_subjects.R.
export MGCV_PREDICTION_OUTPUT="${MGCV_PREDICTION_OUTPUT:-$HOME/scratch/UQAC/meth/results/15_revision/5_cv/5_prediction/1_mgcv_pca_all_subjects_test}"

# Optional comparison with the valid training-only PCA analysis. The collector
# performs this comparison only when the lambda rules match, unless
# ALLOW_LAMBDA_MISMATCH=1 is explicitly supplied.
export COMPARE_WITH_TRAINING_ONLY="${COMPARE_WITH_TRAINING_ONLY:-1}"
export TRAINING_ONLY_PERFORMANCE_FILE="${TRAINING_ONLY_PERFORMANCE_FILE:-$HOME/scratch/UQAC/meth/results/15_revision/5_cv/5_prediction/1_mgcv_primary/collected/mgcv_primary_performance_all_splits.tsv}"
export ALLOW_LAMBDA_MISMATCH="${ALLOW_LAMBDA_MISMATCH:-0}"

echo "MGCV_PREDICTION_OUTPUT=${MGCV_PREDICTION_OUTPUT}"
echo "EXPECTED_SPLITS=${EXPECTED_SPLITS}"
echo "COMPARE_WITH_TRAINING_ONLY=${COMPARE_WITH_TRAINING_ONLY}"
echo "TRAINING_ONLY_PERFORMANCE_FILE=${TRAINING_ONLY_PERFORMANCE_FILE}"

Rscript --vanilla 11_col_mgcv_prediction_cv.R

echo "Finished MGCV all-subject PCA prediction result collection"
date

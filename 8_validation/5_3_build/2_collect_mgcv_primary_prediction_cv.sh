#!/bin/bash
#SBATCH --time=02:50:00
#SBATCH --mem=4G
#SBATCH --account=def-claprise
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=MGCVPredictCollectLogs/slurm-%j.out
#SBATCH --error=MGCVPredictCollectLogs/slurm-%j.err

set -euo pipefail

echo "Start MGCV primary prediction result collection"
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

Rscript --vanilla 2_collect_mgcv_primary_prediction_cv.R

echo "Finished MGCV primary prediction result collection"
date

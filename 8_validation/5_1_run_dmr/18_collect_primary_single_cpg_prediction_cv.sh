#!/bin/bash
#SBATCH --job-name=collect_cpg_pred
#SBATCH --account=def-claprise
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/collect_cpg_pred_%j.out
#SBATCH --error=logs/collect_cpg_pred_%j.err

set -euo pipefail
module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20
export R_LIBS_USER="$HOME/scratch/R/library:$HOME/R/x86_64-pc-linux-gnu-library/4.4"
export PATH_WK="${PATH_WK:-$HOME/scratch/UQAC/meth/}"
export SINGLE_CPG_PREDICTION_OUTPUT="${SINGLE_CPG_PREDICTION_OUTPUT:-$PATH_WK/results/15_revision/5_cv/5_prediction/6_single_cpg_primary}"
export EXPECTED_SPLITS="${EXPECTED_SPLITS:-100}"
export ALLOW_INCOMPLETE="${ALLOW_INCOMPLETE:-0}"
mkdir -p logs
Rscript --vanilla 17_collect_primary_single_cpg_prediction_cv.R

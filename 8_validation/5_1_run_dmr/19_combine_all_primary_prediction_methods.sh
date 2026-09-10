#!/bin/bash
#SBATCH --job-name=combine_pred_methods
#SBATCH --account=def-claprise
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH --output=logs/combine_pred_methods_%j.out
#SBATCH --error=logs/combine_pred_methods_%j.err

set -euo pipefail

module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/scratch/R/library:$HOME/R/x86_64-pc-linux-gnu-library/4.4"

export PATH_WK="${PATH_WK:-$HOME/scratch/UQAC/meth/}"
export ALL_METHOD_COMPARISON_OUTPUT="${ALL_METHOD_COMPARISON_OUTPUT:-$PATH_WK/results/15_revision/5_cv/5_prediction/7_all_method_comparison}"
export EXPECTED_SPLITS="${EXPECTED_SPLITS:-100}"
export OVERWRITE="${OVERWRITE:-1}"

mkdir -p logs

echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Output: $ALL_METHOD_COMPARISON_OUTPUT"

Rscript -e '
pkgs <- c("data.table","pROC","PRROC","ggplot2")
for (p in pkgs) {
  cat(p, ": ")
  if (requireNamespace(p, quietly=TRUE)) {
    cat(find.package(p), "\n")
  } else {
    stop("Required package not found: ", p)
  }
}'

Rscript --vanilla \
  18_combine_all_primary_prediction_methods.R

echo "Finished: $(date)"

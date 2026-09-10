#!/bin/bash
#SBATCH --job-name=GAMdepthCollect
#SBATCH --account=def-claprise
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs_weighted_GAM/GAMdepthCollect_%j.out
#SBATCH --error=logs_weighted_GAM/GAMdepthCollect_%j.err

set -euo pipefail

module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/scratch/R/library:$HOME/R/x86_64-pc-linux-gnu-library/4.4"

export PATH_WK="${PATH_WK:-$HOME/scratch/UQAC/meth/}"

export WEIGHTED_GAM_OUTPUT="${WEIGHTED_GAM_OUTPUT:-$HOME/scratch/UQAC/meth/results/15_revision/7_read_depth_weighted_GAM}"

# Original unweighted outputs from the uploaded 20_run_region_BMI.R.
export UNWEIGHTED_GAM_OUTPUT="${UNWEIGHTED_GAM_OUTPUT:-$HOME/scratch/UQAC/meth/results/11_mgcv/BMI_B1}"

export WEIGHTED_GAM_COMPARISON_OUTPUT="${WEIGHTED_GAM_COMPARISON_OUTPUT:-$HOME/scratch/UQAC/meth/results/15_revision/7_read_depth_weighted_GAM/comparison}"

# Final-DMR sensitivity thresholds. Change here only if your final manuscript
# uses different primary cutoffs.


mkdir -p logs_weighted_GAM

Rscript 23_collect_compare_read_depth_weighted_GAM.R

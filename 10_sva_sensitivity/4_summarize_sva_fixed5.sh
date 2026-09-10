#!/bin/bash
#SBATCH --job-name=sum_sva_fix5
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/sum_sva_fix5_%j.out
#SBATCH --error=logs/sum_sva_fix5_%j.err

set -euo pipefail
module purge
module load StdEnv/2023
module load gcc/12.3
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/scratch/R/library"
export METH_BASE_DIR="$HOME/scratch/UQAC/meth"
export SVA_SENS_OUT="$METH_BASE_DIR/results/15_revision/7_sva_threshold_sensitivity"

mkdir -p logs
Rscript 4_summarize_sva_fixed5.R

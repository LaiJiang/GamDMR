#!/bin/bash
#SBATCH --job-name=prep_sva_sens
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/prep_sva_sens_%j.out
#SBATCH --error=logs/prep_sva_sens_%j.err

set -euo pipefail

module purge
module load StdEnv/2023
module load gcc/12.3
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/scratch/R/library"
export METH_BASE_DIR="$HOME/scratch/UQAC/meth"


# Override these only if your files are elsewhere.
export SVA_METH_FILE="$METH_BASE_DIR/dat/peek_methylation_sva_0035.csv"
export SVA_PHENO_RDATA="$METH_BASE_DIR/scr/11_mgcv/dat/bsmooth_only_pheno_bmi.RData"
export SVA_SENS_OUT="$METH_BASE_DIR/results/15_revision/7_sva_threshold_sensitivity"
export SVA_EXPECTED_N="349"

mkdir -p logs
mkdir -p "$SVA_SENS_OUT"

echo "Start: $(date)"
echo "Host:  $(hostname)"

Rscript 0_prepare_sva_sensitivity_input.R

echo "End:   $(date)"

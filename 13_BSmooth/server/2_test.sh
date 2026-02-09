#!/bin/bash
#SBATCH --time=02:59:00
#SBATCH --mem=24G
#SBATCH --account=def-claprise
#SBATCH --array=1-10%2
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err

set -euo pipefail

echo "Start 2_test.R all regions with BSmooth for AA!"
date
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

# Clean environment, then load R + Bioc bundle
module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20

# Ensure your user library (where bsseq is installed) is active
export R_LIBS_USER="$HOME/R/%p-library/%v"
mkdir -p "$R_LIBS_USER"

# (Optional but nice) keep threaded libs from oversubscribing
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Make sure log directory exists (in case sbatch runs from a fresh dir)
mkdir -p logs

# Run
Rscript 2_test.R "$SLURM_ARRAY_TASK_ID"

echo "Finished 2_test.R all regions with BSmooth for AA!"
date
sleep 30

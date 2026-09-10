#!/bin/bash
#SBATCH --job-name=GAMdepthW
#SBATCH --account=def-claprise
#SBATCH --time=08:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-900%60
#SBATCH --output=logs_weighted_GAM/GAMdepthW_%A_%a.out
#SBATCH --error=logs_weighted_GAM/GAMdepthW_%A_%a.err

set -euo pipefail

module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/scratch/R/library:$HOME/R/x86_64-pc-linux-gnu-library/4.4"

export PATH_WK="${PATH_WK:-$HOME/scratch/UQAC/meth/}"

# This must agree with the SLURM array size.
export N_GAM_JOBS="${N_GAM_JOBS:-900}"

# Explicitly enforce the same minimum read-depth threshold described in the
# manuscript/reviewer response.
export MIN_READ_DEPTH="${MIN_READ_DEPTH:-5}"
export MIN_CPGS="${MIN_CPGS:-10}"

export WEIGHTED_GAM_OUTPUT="${WEIGHTED_GAM_OUTPUT:-$HOME/scratch/UQAC/meth/results/15_revision/7_read_depth_weighted_GAM}"

export OVERWRITE="${OVERWRITE:-0}"

mkdir -p logs_weighted_GAM

Rscript 22_run_region_BMI_read_depth_weighted.R "${SLURM_ARRAY_TASK_ID}"

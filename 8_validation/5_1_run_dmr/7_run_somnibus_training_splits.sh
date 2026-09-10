#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=7G
#SBATCH --account=def-claprise
#SBATCH --array=1-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=SomnibusFitLogs/slurm-%A_%a.out
#SBATCH --error=SomnibusFitLogs/slurm-%A_%a.err

set -euo pipefail

echo "Start SOMNiBUS family-level training-split models"
date
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"

module --force purge
module load StdEnv/2023
module load gcc/12.3
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/scratch/R/library"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# One array task reads one of the 900 prepared bundles.
export N_BUNDLES=900

# Match the previous SOMNiBUS analysis.
export SOMNIBUS_MIN_CPGS=51
export SOMNIBUS_P0=0.003
export SOMNIBUS_P1=0.9
export SOMNIBUS_EPSILON=1e-6
export SOMNIBUS_EPSILON_LAMBDA=1e-3
export SOMNIBUS_MAX_STEP=200

# Resume existing task-specific output. Set OVERWRITE=1 only intentionally.
export OVERWRITE=0
export RETRY_ERRORS=0
export FLUSH_EVERY=10

# Production defaults: all regions and all 100 splits.
export MAX_REGIONS=0
export SPLIT_START=1
export SPLIT_END=100

# Stop cleanly before the 24-hour SLURM limit and resume after resubmission.
export MAX_RUNTIME_MINUTES=1380

Rscript 7_run_somnibus_training_splits.R "${SLURM_ARRAY_TASK_ID}"

echo "Finished SOMNiBUS family-level training-split models"
date

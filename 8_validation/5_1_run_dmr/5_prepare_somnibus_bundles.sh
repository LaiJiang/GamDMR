#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --mem=7G
#SBATCH --account=def-claprise
#SBATCH --array=1-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=PrepLogs/slurm-%A_%a.out
#SBATCH --error=PrepLogs/slurm-%A_%a.err

set -euo pipefail

echo "Start SOMNiBUS bundle preparation"
date
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"

module purge
module load StdEnv/2023
module load gcc/12.3
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/scratch/R/library"


export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Must match the upper bound of --array.
export N_JOBS=900

# Matches the previous splitDataByRegion() stage.
export SOMNIBUS_GAP_BP=250
export SOMNIBUS_MIN_CPGS=51
export SOMNIBUS_MAX_CPGS=2000

# Change to 1 only when intentionally replacing existing outputs.
export OVERWRITE=0

Rscript 5_prepare_somnibus_bundles.R "${SLURM_ARRAY_TASK_ID}"

echo "Finished SOMNiBUS bundle preparation"
date

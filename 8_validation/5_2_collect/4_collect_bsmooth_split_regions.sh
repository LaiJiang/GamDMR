#!/bin/bash
#SBATCH --time=02:50:00
#SBATCH --mem=7G
#SBATCH --account=def-claprise
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=BSmoothCollectLogs/slurm-%j.out
#SBATCH --error=BSmoothCollectLogs/slurm-%j.err

set -euo pipefail

echo "Start BSmooth family-split result collection"
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

export SPLIT_START="${SPLIT_START:-1}"
export SPLIT_END="${SPLIT_END:-100}"
export EXPECTED_JOBS="${EXPECTED_JOBS:-990}"
export OVERWRITE="${OVERWRITE:-0}"

Rscript --vanilla 4_collect_bsmooth_split_regions.R

echo "Finished BSmooth family-split result collection"
date

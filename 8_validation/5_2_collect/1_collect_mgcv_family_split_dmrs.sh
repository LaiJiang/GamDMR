#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --account=def-claprise
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=MGCVCollectLogs/slurm-%j.out
#SBATCH --error=MGCVCollectLogs/slurm-%j.err

set -euo pipefail

echo "Start collecting MGCV family-split DMR results"
date
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"

module --force purge
module load StdEnv/2023
module load gcc/12.3
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

Rscript --vanilla 4_collect_mgcv_family_split_dmrs.R

echo "Finished collecting MGCV family-split DMR results"
date

#!/bin/bash
#SBATCH --job-name=readDepthCollect
#SBATCH --account=def-claprise
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/read_depth_collect_%j.out
#SBATCH --error=logs/read_depth_collect_%j.err

set -euo pipefail

module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/scratch/R/library:$HOME/R/x86_64-pc-linux-gnu-library/4.4"

export PATH_WK="${PATH_WK:-$HOME/scratch/UQAC/meth/}"
export N_READ_DEPTH_JOBS="${N_READ_DEPTH_JOBS:-900}"
export ALLOW_INCOMPLETE="${ALLOW_INCOMPLETE:-0}"

mkdir -p logs

Rscript 21_collect_read_depth_variability.R

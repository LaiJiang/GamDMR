#!/bin/bash
#SBATCH --job-name=readDepth
#SBATCH --account=def-claprise
#SBATCH --time=02:50:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-900%60
#SBATCH --output=logs/read_depth_%A_%a.out
#SBATCH --error=logs/read_depth_%A_%a.err

set -euo pipefail

module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/scratch/R/library:$HOME/R/x86_64-pc-linux-gnu-library/4.4"

export PATH_WK="${PATH_WK:-$HOME/scratch/UQAC/meth/}"
export N_READ_DEPTH_JOBS="${N_READ_DEPTH_JOBS:-900}"
export MAX_HIST_DEPTH="${MAX_HIST_DEPTH:-200}"
export ANALYTIC_ONLY="${ANALYTIC_ONLY:-1}"
export OVERWRITE="${OVERWRITE:-0}"

mkdir -p logs

Rscript 1_run_read_depth_variability_chunks.R "${SLURM_ARRAY_TASK_ID}"

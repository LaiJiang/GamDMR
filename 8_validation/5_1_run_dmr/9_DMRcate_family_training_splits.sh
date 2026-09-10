#!/bin/bash
#SBATCH --time=02:59:00
#SBATCH --mem=4G
#SBATCH --account=def-claprise
#SBATCH --array=1-900
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=DMRcateTrainLogs/slurm-%A_%a.out
#SBATCH --error=DMRcateTrainLogs/slurm-%A_%a.err

set -euo pipefail

echo "Start DMRcate family-level training-split analysis"
date
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/%p-library/%v}"
mkdir -p "$R_LIBS_USER"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Production defaults. Command-line --export values can override these.
export N_JOBS="${N_JOBS:-900}"
export SPLIT_START="${SPLIT_START:-1}"
export SPLIT_END="${SPLIT_END:-100}"
export MAX_REGIONS="${MAX_REGIONS:-0}"
export MAX_RUNTIME_MINUTES="${MAX_RUNTIME_MINUTES:-1380}"
export OVERWRITE="${OVERWRITE:-0}"
export RETRY_ERRORS="${RETRY_ERRORS:-0}"

# DMRcate settings matching the previous stringent analysis.
export DMRCATE_MIN_CPGS="${DMRCATE_MIN_CPGS:-10}"
export DMRCATE_MIN_COV="${DMRCATE_MIN_COV:-5}"
export DMRCATE_COVERAGE_FRACTION="${DMRCATE_COVERAGE_FRACTION:-0.70}"
export DMRCATE_SEED_FDR="${DMRCATE_SEED_FDR:-0.05}"
export DMRCATE_MIN_SEED_CPGS="${DMRCATE_MIN_SEED_CPGS:-5}"
export DMRCATE_EXTRACT_FDR="${DMRCATE_EXTRACT_FDR:-0.01}"
export DMRCATE_MIN_ABS_MD="${DMRCATE_MIN_ABS_MD:-0.05}"
export DMRCATE_MIN_WIDTH="${DMRCATE_MIN_WIDTH:-250}"
export DMRCATE_LAMBDA="${DMRCATE_LAMBDA:-1000}"
export DMRCATE_C="${DMRCATE_C:-2}"

Rscript --vanilla 9_DMRcate_family_training_splits.R \
  "${SLURM_ARRAY_TASK_ID}"

echo "Finished DMRcate family-level training-split analysis"
date

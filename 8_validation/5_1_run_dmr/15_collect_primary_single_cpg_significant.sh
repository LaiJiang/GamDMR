#!/bin/bash
#SBATCH --job-name=collect_cpg_sig
#SBATCH --account=def-claprise
#SBATCH --array=1-99%50
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=06:00:00
#SBATCH --output=logs/collect_cpg_sig_%A_%a.out
#SBATCH --error=logs/collect_cpg_sig_%A_%a.err

set -euo pipefail

module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/scratch/R/library:$HOME/R/x86_64-pc-linux-gnu-library/4.4"

export PATH_WK="${PATH_WK:-$HOME/scratch/UQAC/meth/}"
export SINGLE_CPG_ASSOC_OUTPUT="${SINGLE_CPG_ASSOC_OUTPUT:-$PATH_WK/results/15_revision/5_cv/6_single_cpg/1_training_associations}"
export SINGLE_CPG_COLLECT_ROOT="${SINGLE_CPG_COLLECT_ROOT:-$PATH_WK/results/15_revision/5_cv/6_single_cpg/2_collected_significant}"

export N_ASSOC_JOBS="${N_ASSOC_JOBS:-990}"
export N_COLLECT_JOBS="${N_COLLECT_JOBS:-99}"
export CANDIDATE_P_MAX="${CANDIDATE_P_MAX:-0.001}"
export ALLOW_INCOMPLETE="${ALLOW_INCOMPLETE:-0}"
export OVERWRITE_COLLECT="${OVERWRITE_COLLECT:-1}"

mkdir -p logs
mkdir -p "$SINGLE_CPG_COLLECT_ROOT"

Rscript --vanilla \
  15_collect_primary_single_cpg_significant.R \
  "$SLURM_ARRAY_TASK_ID"

#!/bin/bash
#SBATCH --job-name=finalize_cpg_sig
#SBATCH --account=def-claprise
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/finalize_cpg_sig_%j.out
#SBATCH --error=logs/finalize_cpg_sig_%j.err

set -euo pipefail

module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/scratch/R/library:$HOME/R/x86_64-pc-linux-gnu-library/4.4"

export PATH_WK="${PATH_WK:-$HOME/scratch/UQAC/meth/}"
export SINGLE_CPG_COLLECT_ROOT="${SINGLE_CPG_COLLECT_ROOT:-$PATH_WK/results/15_revision/5_cv/6_single_cpg/2_collected_significant}"

export N_COLLECT_JOBS="${N_COLLECT_JOBS:-99}"
export SIGNIFICANCE_P="${SIGNIFICANCE_P:-1e-5}"
export CANDIDATE_P_MAX="${CANDIDATE_P_MAX:-0.001}"
export BONF_ALPHA="${BONF_ALPHA:-0.05}"
export EXPECTED_SPLITS="${EXPECTED_SPLITS:-100}"
export ALLOW_INCOMPLETE="${ALLOW_INCOMPLETE:-0}"

mkdir -p logs

Rscript --vanilla \
  16_finalize_primary_single_cpg_selection.R

#!/bin/bash
#SBATCH --job-name=prep_cpg_primary
#SBATCH --account=def-claprise
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=06:00:00
#SBATCH --output=logs/prep_cpg_primary_%j.out
#SBATCH --error=logs/prep_cpg_primary_%j.err

set -euo pipefail

module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/scratch/R/library:$HOME/R/x86_64-pc-linux-gnu-library/4.4"

export PATH_WK="${PATH_WK:-$HOME/scratch/UQAC/meth/}"
export METH_CHUNK_DIR="${METH_CHUNK_DIR:-$PATH_WK/data/meth_split}"

export PRIMARY_CPG_MANIFEST_ROOT="${PRIMARY_CPG_MANIFEST_ROOT:-$PATH_WK/results/15_revision/5_cv/6_single_cpg/0_primary_all_cpg_manifest}"

export N_JOBS="${N_JOBS:-990}"
export EXPECTED_CPGS="${EXPECTED_CPGS:-4609564}"
export OVERWRITE_MANIFEST="${OVERWRITE_MANIFEST:-0}"

mkdir -p logs

echo "Started: $(date)"
echo "Host: $(hostname)"
echo "METH_CHUNK_DIR=$METH_CHUNK_DIR"
echo "PRIMARY_CPG_MANIFEST_ROOT=$PRIMARY_CPG_MANIFEST_ROOT"
echo "N_JOBS=$N_JOBS"

Rscript --vanilla 12_prepare_primary_single_cpg_manifest.R

echo "Finished: $(date)"

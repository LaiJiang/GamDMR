#!/bin/bash
#SBATCH --time=02:50:00
#SBATCH --mem=7G
#SBATCH --account=def-claprise
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=SomnibusCollectLogs/slurm-%j.out
#SBATCH --error=SomnibusCollectLogs/slurm-%j.err

set -euo pipefail

echo "Start SOMNiBUS family-split result collection"
date
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"

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

export SPLIT_START="${SPLIT_START:-1}"
export SPLIT_END="${SPLIT_END:-100}"

# Match the previous 16_somnibus_regions.csv rule.
export SOMNIBUS_REGION_P_THRESHOLD="${SOMNIBUS_REGION_P_THRESHOLD:-0.01}"

# Number of task-specific SOMNiBUS result files expected.
export EXPECTED_BUNDLES="${EXPECTED_BUNDLES:-900}"

# Set to 1 only when intentionally replacing previous collector outputs.
export OVERWRITE="${OVERWRITE:-0}"

Rscript --vanilla 2_collect_somnibus_split_regions.R

echo "Finished SOMNiBUS family-split result collection"
date

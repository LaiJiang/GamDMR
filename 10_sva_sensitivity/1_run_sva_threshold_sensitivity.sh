#!/bin/bash
#SBATCH --job-name=sva_vfilter
#SBATCH --array=1-3
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/sva_vfilter_%A_%a.out
#SBATCH --error=logs/sva_vfilter_%A_%a.err

set -euo pipefail

module purge
module load StdEnv/2023
module load gcc/12.3
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/scratch/R/library"
export METH_BASE_DIR="$HOME/scratch/UQAC/meth"
export SVA_SENS_OUT="$METH_BASE_DIR/results/15_revision/7_sva_threshold_sensitivity"
export SVA_SEED="20260816"

mkdir -p logs
mkdir -p "$SVA_SENS_OUT"

case "${SLURM_ARRAY_TASK_ID}" in
  1) VF=10000 ;;
  2) VF=25000 ;;
  3) VF=50000 ;;
  *) echo "Invalid array task"; exit 1 ;;
esac

echo "Start:   $(date)"
echo "Host:    $(hostname)"
echo "vfilter: ${VF}"

Rscript 1_run_sva_threshold_sensitivity.R "${VF}"

echo "End:     $(date)"

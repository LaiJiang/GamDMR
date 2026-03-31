#!/bin/bash
#SBATCH --time=02:59:00
#SBATCH --mem=7G
#SBATCH --account=XXXXXXX
#SBATCH --array=1-300
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err

set -euo pipefail

echo "Starting regional GAM-DMR analysis"
date

module purge
module load StdEnv/2023
module load gcc/12.3
module load r/4.4.0
module load r-bundle-bioconductor/3.20

Rscript 1_run_region_clean.R "${SLURM_ARRAY_TASK_ID}"

echo "Finished regional GAM-DMR analysis"
date

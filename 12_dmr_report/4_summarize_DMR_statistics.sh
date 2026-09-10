#!/usr/bin/env bash
#SBATCH --job-name=DMRsummary
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/%u/scratch/UQAC/meth/results/15_revision/9_dmr_report/DMR_statistics/logs/DMRsummary_%j.out
#SBATCH --error=/home/%u/scratch/UQAC/meth/results/15_revision/9_dmr_report/DMR_statistics/logs/DMRsummary_%j.err

set -euo pipefail

module purge
module load StdEnv/2023
module load gcc/12.3
module load r/4.4.0e-bioconductor/3.20

export R_LIBS_USER="${HOME}/scratch/R/library"

BASE="${HOME}/scratch/UQAC/meth"
SCRIPT_DIR="${BASE}/scr/15_revision/9_dmr_report"
OUT_ROOT="${BASE}/results/15_revision/9_dmr_report/DMR_statistics"

mkdir -p "${OUT_ROOT}/logs"

cd "${SCRIPT_DIR}"

echo "============================================================"
echo "Summarizing DMR characteristics"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Host:   $(hostname)"
echo "Start:  $(date)"
echo "============================================================"

Rscript 4_summarize_DMR_statistics.R "${OUT_ROOT}"

echo
echo "Finished: $(date)"

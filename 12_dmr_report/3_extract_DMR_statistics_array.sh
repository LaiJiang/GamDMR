#!/usr/bin/env bash
#SBATCH --job-name=DMRchunk
#SBATCH --time=04:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/%u/scratch/UQAC/meth/results/15_revision/9_dmr_report/DMR_statistics/logs/DMRchunk_%A_%a.out
#SBATCH --error=/home/%u/scratch/UQAC/meth/results/15_revision/9_dmr_report/DMR_statistics/logs/DMRchunk_%A_%a.err

set -euo pipefail

module purge
module load StdEnv/2023
module load gcc/12.3
module load r/4.4.0e-bioconductor/3.20

export R_LIBS_USER="${HOME}/scratch/R/library"

BASE="${HOME}/scratch/UQAC/meth"
SCRIPT_DIR="${BASE}/scr/15_revision/9_dmr_report"
OUT_ROOT="${BASE}/results/15_revision/9_dmr_report/DMR_statistics"
CHUNK_DIR="${BASE}/data/meth_split"

# Passed by 2_submit_DMR_pipeline.sh.
CHUNKS_PER_TASK="${CHUNKS_PER_TASK:-10}"

mkdir -p "${OUT_ROOT}/logs"
mkdir -p "${OUT_ROOT}/partials"

cd "${SCRIPT_DIR}"

echo "============================================================"
echo "Bundled DMR chunk array task"
echo "Job ID:          ${SLURM_JOB_ID:-NA}"
echo "Array job:       ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array task:      ${SLURM_ARRAY_TASK_ID:-NA}"
echo "Chunks per task: ${CHUNKS_PER_TASK}"
echo "Host:            $(hostname)"
echo "Start:           $(date)"
echo "============================================================"

Rscript 3_extract_DMR_statistics_array.R \
    "${SLURM_ARRAY_TASK_ID}" \
    "${OUT_ROOT}" \
    "${CHUNK_DIR}" \
    "${CHUNKS_PER_TASK}"

echo
echo "Finished: $(date)"

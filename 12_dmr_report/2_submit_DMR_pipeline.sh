#!/usr/bin/env bash

# =============================================================================
# 2_submit_DMR_pipeline.sh
#
# Master launcher for chunk-BUNDLED SLURM array.
#
# Why bundle chunks?
# ------------------
# Rorqual counts every array task toward the user's submitted-job limit.
# The DMR manifest currently requires 1733 source chunks, so submitting
# --array=1-1733%100 violates the normal-account MaxSubmitJobs limit even
# though only 100 would run concurrently.
#
# This launcher therefore groups multiple source chunks into each array task.
#
# Defaults:
#   CHUNKS_PER_TASK=10
#   MAX_CONCURRENT=100
#
# For 1733 chunks:
#   ceil(1733 / 10) = 174 submitted array tasks
#
# Usage:
#   ./2_submit_DMR_pipeline.sh
#
# Optional:
#   CHUNKS_PER_TASK=20 MAX_CONCURRENT=50 ./2_submit_DMR_pipeline.sh
# =============================================================================

set -euo pipefail

BASE="${HOME}/scratch/UQAC/meth"
SCRIPT_DIR="${BASE}/scr/15_revision/9_dmr_report"
OUT_ROOT="${BASE}/results/15_revision/9_dmr_report/DMR_statistics"

DMR_FILE="${BASE}/results/15_revision/9_dmr_report/DMR_input/all_DMRs_long.csv"
REGION_FILE="${BASE}/scr/11_mgcv/dat/region_file_1_chunk.csv"
PHENO_RDATA="${BASE}/scr/11_mgcv/dat/18_pheno_BMI.RData"

MAX_CONCURRENT="${MAX_CONCURRENT:-100}"
CHUNKS_PER_TASK="${CHUNKS_PER_TASK:-10}"

if ! [[ "${MAX_CONCURRENT}" =~ ^[0-9]+$ ]] || [[ "${MAX_CONCURRENT}" -lt 1 ]]; then
    echo "ERROR: MAX_CONCURRENT must be a positive integer." >&2
    exit 1
fi

if ! [[ "${CHUNKS_PER_TASK}" =~ ^[0-9]+$ ]] || [[ "${CHUNKS_PER_TASK}" -lt 1 ]]; then
    echo "ERROR: CHUNKS_PER_TASK must be a positive integer." >&2
    exit 1
fi

mkdir -p "${OUT_ROOT}/work"
mkdir -p "${OUT_ROOT}/partials"
mkdir -p "${OUT_ROOT}/logs"

cd "${SCRIPT_DIR}"

module purge
module load StdEnv/2023
module load gcc/12.3
module load r/4.4.0e-bioconductor/3.20

export R_LIBS_USER="${HOME}/scratch/R/library"

echo "============================================================"
echo "Preparing chunk manifest"
echo "============================================================"

Rscript 2_prepare_DMR_chunk_manifest.R \
    "${DMR_FILE}" \
    "${REGION_FILE}" \
    "${PHENO_RDATA}" \
    "${OUT_ROOT}"

MANIFEST="${OUT_ROOT}/work/chunk_manifest.csv"

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: manifest was not created: ${MANIFEST}" >&2
    exit 1
fi

# CSV has one header line.
N_CHUNKS=$(( $(wc -l < "${MANIFEST}") - 1 ))

if [[ "${N_CHUNKS}" -lt 1 ]]; then
    echo "ERROR: no source chunks found in manifest." >&2
    exit 1
fi

# Ceiling division.
N_ARRAY_TASKS=$(( (N_CHUNKS + CHUNKS_PER_TASK - 1) / CHUNKS_PER_TASK ))

# There is no benefit in allowing more simultaneous tasks than exist.
if [[ "${MAX_CONCURRENT}" -gt "${N_ARRAY_TASKS}" ]]; then
    ARRAY_THROTTLE="${N_ARRAY_TASKS}"
else
    ARRAY_THROTTLE="${MAX_CONCURRENT}"
fi

echo
echo "Source methylation chunks:    ${N_CHUNKS}"
echo "Chunks per array task:        ${CHUNKS_PER_TASK}"
echo "Total submitted array tasks:  ${N_ARRAY_TASKS}"
echo "Maximum concurrently running: ${ARRAY_THROTTLE}"
echo

if [[ "${N_ARRAY_TASKS}" -ge 1000 ]]; then
    echo "WARNING: ${N_ARRAY_TASKS} array tasks is still near/above the normal"
    echo "Rorqual submitted-job limit. Increase CHUNKS_PER_TASK."
    echo
fi

ARRAY_JOB_ID=$(sbatch \
    --parsable \
    --array="1-${N_ARRAY_TASKS}%${ARRAY_THROTTLE}" \
    --export="ALL,CHUNKS_PER_TASK=${CHUNKS_PER_TASK}" \
    3_extract_DMR_statistics_array.sh)

echo "Submitted bundled chunk-array job: ${ARRAY_JOB_ID}"

SUMMARY_JOB_ID=$(sbatch \
    --parsable \
    --dependency="afterok:${ARRAY_JOB_ID}" \
    4_summarize_DMR_statistics.sh)

echo "Submitted dependent summary job: ${SUMMARY_JOB_ID}"

echo
echo "============================================================"
echo "Pipeline submitted"
echo "============================================================"
echo "Source chunks:      ${N_CHUNKS}"
echo "Chunks/task:        ${CHUNKS_PER_TASK}"
echo "Array tasks:        ${N_ARRAY_TASKS}"
echo "Array job:          ${ARRAY_JOB_ID}"
echo "Summary job:        ${SUMMARY_JOB_ID}"
echo
echo "Monitor:"
echo "  squeue -u ${USER}"
echo
echo "Array logs:"
echo "  ${OUT_ROOT}/logs/DMRchunk_${ARRAY_JOB_ID}_*.out"
echo
echo "When everything finishes:"
echo "  cat ${OUT_ROOT}/DMR_characteristics_report.txt"
echo "============================================================"

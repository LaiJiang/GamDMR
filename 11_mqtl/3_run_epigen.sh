#!/usr/bin/env bash

# =============================================================================
# 3_run_epigen.sh
#
# Run EPIGEN meQTL overlap analysis for the four DMR methods:
#   1. GAM-DMR
#   2. SOMNiBUS
#   3. BSmooth
#   4. DMRcate
#
# Then automatically create:
#   - one method-level supplementary comparison table
#   - one combined DMR-level supplementary table
#
# The four EPIGEN scans are run SEQUENTIALLY on purpose. meQTL_full.txt.gz is
# ~1.8 GB compressed; running four scans simultaneously on the same local disk
# usually causes unnecessary I/O contention.
# =============================================================================

set -euo pipefail

# ------------------------------- Paths ---------------------------------------

BASE="/mnt/c/Per/LaiJiang/Project/UQAC/meth"

SCRIPT_DIR="${BASE}/scr/15_revision/8_mqtl"
DMR_INPUT_DIR="${BASE}/results/15_revision/8_mqtl/DMR_input"
EPIGEN_FILE="${BASE}/dat/EPIGEN/meQTL_full.txt.gz"

OVERLAP_SCRIPT="${SCRIPT_DIR}/2_epigen_dmr_overlap.py"
SUMMARY_SCRIPT="${SCRIPT_DIR}/4_summarize_epigen.py"

OUT_ROOT="${BASE}/results/15_revision/8_mqtl/EPIGEN"

mkdir -p "${OUT_ROOT}"

# ----------------------------- Input files -----------------------------------

GAM_INPUT="${DMR_INPUT_DIR}/gam_dmr_input.csv"
SOMNIBUS_INPUT="${DMR_INPUT_DIR}/somnibus_dmr_input.csv"
BSMOOTH_INPUT="${DMR_INPUT_DIR}/BSmooth_dmr_input.csv"
DMRCATE_INPUT="${DMR_INPUT_DIR}/DMRcate_dmr_input.csv"

# ----------------------------- Sanity checks ---------------------------------

echo "============================================================"
echo "EPIGEN meQTL overlap analysis: four DMR methods"
echo "============================================================"
echo "Started: $(date)"
echo "EPIGEN: ${EPIGEN_FILE}"
echo "Output: ${OUT_ROOT}"
echo

for f in \
    "${OVERLAP_SCRIPT}" \
    "${SUMMARY_SCRIPT}" \
    "${EPIGEN_FILE}" \
    "${GAM_INPUT}" \
    "${SOMNIBUS_INPUT}" \
    "${BSMOOTH_INPUT}" \
    "${DMRCATE_INPUT}"
do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required file does not exist:"
        echo "  ${f}"
        exit 1
    fi
done

# Validate the gzip file before spending time on four analyses.
echo "Checking EPIGEN gzip integrity..."
gzip -t "${EPIGEN_FILE}"
echo "EPIGEN gzip file is valid."
echo

# ----------------------------- Run helper ------------------------------------

run_method () {
    local method_label="$1"
    local method_slug="$2"
    local dmr_file="$3"
    local out_dir="${OUT_ROOT}/${method_slug}"

    mkdir -p "${out_dir}"

    echo
    echo "============================================================"
    echo "Running ${method_label}"
    echo "Input:  ${dmr_file}"
    echo "Output: ${out_dir}"
    echo "Start:  $(date)"
    echo "============================================================"

    python3 "${OVERLAP_SCRIPT}" \
        --dmr "${dmr_file}" \
        --meqtl "${EPIGEN_FILE}" \
        --out "${out_dir}" \
        2>&1 | tee "${out_dir}/run.log"

    echo "Finished ${method_label}: $(date)"
}

# ----------------------------- Four methods ----------------------------------

run_method "GAM-DMR"  "GAM_DMR"  "${GAM_INPUT}"
run_method "SOMNiBUS" "SOMNiBUS" "${SOMNIBUS_INPUT}"
run_method "BSmooth"  "BSmooth"  "${BSMOOTH_INPUT}"
run_method "DMRcate"  "DMRcate"  "${DMRCATE_INPUT}"

# ----------------------------- Summarize -------------------------------------

echo
echo "============================================================"
echo "All four analyses completed."
echo "Creating supplementary comparison tables..."
echo "============================================================"

python3 "${SUMMARY_SCRIPT}" \
    --root "${OUT_ROOT}"

echo
echo "============================================================"
echo "EPIGEN four-method analysis completed successfully."
echo "Finished: $(date)"
echo
echo "Main paper/supplement outputs:"
echo "  ${OUT_ROOT}/Supplementary_Table_EPIGEN_meQTL_overlap_by_method.csv"
echo "  ${OUT_ROOT}/Supplementary_Table_EPIGEN_meQTL_overlap_by_DMR.csv"
echo
echo "Method-specific detailed results:"
echo "  ${OUT_ROOT}/GAM_DMR/"
echo "  ${OUT_ROOT}/SOMNiBUS/"
echo "  ${OUT_ROOT}/BSmooth/"
echo "  ${OUT_ROOT}/DMRcate/"
echo "============================================================"

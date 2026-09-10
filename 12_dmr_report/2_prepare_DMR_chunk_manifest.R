#!/usr/bin/env Rscript

# =============================================================================
# 2_prepare_DMR_chunk_manifest.R
#
# Prepare a chunk-based work manifest for descriptive DMR methylation statistics.
#
# Input DMR file:
#   ~/scratch/UQAC/meth/results/15_revision/9_dmr_report/DMR_input/all_DMRs_long.csv
#
# Expected columns:
#   method, DMR_chr, DMR_region_start, DMR_region_end, DMR_width_bp
#
# Existing project inputs:
#   ~/scratch/UQAC/meth/scr/11_mgcv/dat/region_file_1_chunk.csv
#   ~/scratch/UQAC/meth/scr/11_mgcv/dat/18_pheno_BMI.RData
#
# Outputs under:
#   ~/scratch/UQAC/meth/results/15_revision/9_dmr_report/DMR_statistics/work/
#
#   dmr_table.rds
#   dmr_chunk_map.csv
#   chunk_manifest.csv
#   analysis_subjects.csv
#   preparation_report.txt
#
# Each unique methylation chunk becomes one SLURM-array task.
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

BASE <- path.expand("~/scratch/UQAC/meth")

DMR_FILE <- if (length(args) >= 1L) args[1L] else file.path(
    BASE,
    "results/15_revision/9_dmr_report/DMR_input/all_DMRs_long.csv"
)

REGION_FILE <- if (length(args) >= 2L) args[2L] else file.path(
    BASE,
    "scr/11_mgcv/dat/region_file_1_chunk.csv"
)

PHENO_RDATA <- if (length(args) >= 3L) args[3L] else file.path(
    BASE,
    "scr/11_mgcv/dat/18_pheno_BMI.RData"
)

OUT_ROOT <- if (length(args) >= 4L) args[4L] else file.path(
    BASE,
    "results/15_revision/9_dmr_report/DMR_statistics"
)

WORK_DIR <- file.path(OUT_ROOT, "work")
PARTIAL_DIR <- file.path(OUT_ROOT, "partials")
LOG_DIR <- file.path(OUT_ROOT, "logs")

dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PARTIAL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

norm_chr <- function(x) {
    x <- as.character(x)
    x <- sub("^chr", "", x, ignore.case = TRUE)
    toupper(trimws(x))
}

cat("============================================================\n")
cat("Preparing DMR -> methylation-chunk manifest\n")
cat("============================================================\n")
cat("DMR file:     ", DMR_FILE, "\n", sep = "")
cat("Region index: ", REGION_FILE, "\n", sep = "")
cat("Phenotype:    ", PHENO_RDATA, "\n", sep = "")
cat("Output root:  ", OUT_ROOT, "\n\n", sep = "")

for (f in c(DMR_FILE, REGION_FILE, PHENO_RDATA)) {
    if (!file.exists(f)) stop("Required file does not exist: ", f)
}

# -------------------------------------------------------------------------
# DMR table
# -------------------------------------------------------------------------

dmr <- fread(DMR_FILE)

required_dmr <- c(
    "method",
    "DMR_chr",
    "DMR_region_start",
    "DMR_region_end",
    "DMR_width_bp"
)

missing_dmr <- setdiff(required_dmr, names(dmr))
if (length(missing_dmr)) {
    stop("DMR file missing columns: ", paste(missing_dmr, collapse = ", "))
}

dmr[, DMR_chr := norm_chr(DMR_chr)]
dmr[, DMR_region_start := as.integer(DMR_region_start)]
dmr[, DMR_region_end := as.integer(DMR_region_end)]
dmr[, DMR_width_bp := as.integer(DMR_width_bp)]

if (anyNA(dmr[, ..required_dmr])) {
    stop("DMR input contains missing values in required columns.")
}

if (any(dmr$DMR_region_start > dmr$DMR_region_end)) {
    stop("At least one DMR has start > end.")
}

# Recompute width from coordinates so the final table is internally consistent.
dmr[, DMR_width_bp := DMR_region_end - DMR_region_start + 1L]

# Stable unique identifiers for the 12,030 method-specific DMR rows.
dmr[, DMR_row_id := .I]
dmr[, DMR_method_index := seq_len(.N), by = method]
dmr[, DMR_method_id := sprintf(
    "%s_%05d",
    gsub("[^A-Za-z0-9]+", "_", method),
    DMR_method_index
)]

cat("DMR rows loaded: ", nrow(dmr), "\n", sep = "")
print(dmr[, .N, by = method])
cat("\n")

# -------------------------------------------------------------------------
# Analytic subjects
# -------------------------------------------------------------------------

env <- new.env(parent = emptyenv())
load(PHENO_RDATA, envir = env)

if (!exists("pheno_file", envir = env, inherits = FALSE)) {
    stop("Object 'pheno_file' not found in ", PHENO_RDATA)
}

pheno <- as.data.table(get("pheno_file", envir = env))

if (!all(c("ID", "AA_only") %in% names(pheno))) {
    stop("pheno_file must contain ID and AA_only.")
}

pheno[, ID := as.character(ID)]
pheno[, AA_only := suppressWarnings(as.integer(as.character(AA_only)))]
pheno <- pheno[AA_only %in% c(0L, 1L), .(ID, AA_only)]
pheno <- unique(pheno, by = "ID")

if (anyDuplicated(pheno$ID)) stop("Duplicate IDs remain in phenotype table.")

pheno[, subject_index := .I]
setcolorder(pheno, c("subject_index", "ID", "AA_only"))

fwrite(pheno, file.path(WORK_DIR, "analysis_subjects.csv"))

cat("Analytic subjects: ", nrow(pheno), "\n", sep = "")
cat("AA cases:          ", pheno[AA_only == 1L, .N], "\n", sep = "")
cat("Controls:          ", pheno[AA_only == 0L, .N], "\n\n", sep = "")

# -------------------------------------------------------------------------
# Region -> chunk index
# -------------------------------------------------------------------------

region_header <- names(fread(REGION_FILE, nrows = 0L))
required_region <- c(
    "data_chunk_id", "chr", "region_start", "region_end"
)

missing_region <- setdiff(required_region, region_header)
if (length(missing_region)) {
    stop(
        "region_file_1_chunk.csv missing columns: ",
        paste(missing_region, collapse = ", ")
    )
}

reg <- fread(
    REGION_FILE,
    select = required_region,
    showProgress = FALSE
)

reg[, chr := norm_chr(chr)]
reg[, data_chunk_id := as.integer(data_chunk_id)]
reg[, region_start := as.integer(region_start)]
reg[, region_end := as.integer(region_end)]

reg <- reg[
    !is.na(data_chunk_id) &
    !is.na(chr) &
    !is.na(region_start) &
    !is.na(region_end)
]

if (!nrow(reg)) stop("No valid rows in region_file_1_chunk.csv.")

# The same chunk may own many predefined regions. Keep all region intervals
# for accurate "any overlap" mapping, then de-duplicate DMR/chunk pairs.
setkey(reg, chr, region_start, region_end)

dmr_overlap <- dmr[, .(
    DMR_row_id,
    chr = DMR_chr,
    dmr_start = DMR_region_start,
    dmr_end = DMR_region_end
)]

mapped <- foverlaps(
    dmr_overlap,
    reg,
    by.x = c("chr", "dmr_start", "dmr_end"),
    by.y = c("chr", "region_start", "region_end"),
    type = "any",
    nomatch = 0L
)

chunk_map <- unique(mapped[, .(
    DMR_row_id,
    data_chunk_id
)])

setorder(chunk_map, data_chunk_id, DMR_row_id)

unmapped_ids <- setdiff(dmr$DMR_row_id, chunk_map$DMR_row_id)

if (length(unmapped_ids)) {
    unmapped_file <- file.path(WORK_DIR, "DMRs_without_chunk_mapping.csv")
    fwrite(dmr[DMR_row_id %in% unmapped_ids], unmapped_file)

    stop(
        length(unmapped_ids),
        " DMR rows could not be mapped to a methylation chunk. ",
        "See: ", unmapped_file
    )
}

# Number of chunks intersected by each DMR is useful QC.
nchunks <- chunk_map[, .(n_source_chunks = uniqueN(data_chunk_id)), by = DMR_row_id]
dmr <- merge(dmr, nchunks, by = "DMR_row_id", all.x = TRUE, sort = FALSE)
setorder(dmr, DMR_row_id)

# -------------------------------------------------------------------------
# Manifest: exactly one array task per required chunk
# -------------------------------------------------------------------------

manifest <- unique(chunk_map[, .(data_chunk_id)])
setorder(manifest, data_chunk_id)
manifest[, array_task_id := .I]
manifest[, n_DMR_rows := chunk_map[.SD, on = "data_chunk_id", .N, by = .EACHI]$N]

setcolorder(manifest, c("array_task_id", "data_chunk_id", "n_DMR_rows"))

fwrite(chunk_map, file.path(WORK_DIR, "dmr_chunk_map.csv"))
fwrite(manifest, file.path(WORK_DIR, "chunk_manifest.csv"))
saveRDS(dmr, file.path(WORK_DIR, "dmr_table.rds"))

# Remove stale partials from a previous run so the final collector can require
# exactly the current manifest.
old_partials <- list.files(
    PARTIAL_DIR,
    pattern = "^chunk_[0-9]+\\.rds$",
    full.names = TRUE
)
if (length(old_partials)) file.remove(old_partials)

# -------------------------------------------------------------------------
# Report
# -------------------------------------------------------------------------

report <- c(
    "DMR CHUNK-MANIFEST PREPARATION REPORT",
    "=====================================",
    paste0("DMR rows: ", format(nrow(dmr), big.mark = ",")),
    paste0("Unique source chunks required: ", nrow(manifest)),
    paste0("DMR-chunk pairs: ", format(nrow(chunk_map), big.mark = ",")),
    paste0(
        "DMRs spanning >1 source chunk: ",
        dmr[n_source_chunks > 1L, .N]
    ),
    paste0("Analytic subjects: ", nrow(pheno)),
    paste0("AA cases: ", pheno[AA_only == 1L, .N]),
    paste0("Controls: ", pheno[AA_only == 0L, .N]),
    "",
    "DMR counts by method:"
)

method_lines <- dmr[, paste0(
    "  ", method, ": ", format(.N, big.mark = ",")
), by = method]$V1

report <- c(
    report,
    method_lines,
    "",
    paste0("Manifest: ", file.path(WORK_DIR, "chunk_manifest.csv")),
    paste0("DMR-chunk map: ", file.path(WORK_DIR, "dmr_chunk_map.csv")),
    paste0("DMR table RDS: ", file.path(WORK_DIR, "dmr_table.rds")),
    paste0("Subjects: ", file.path(WORK_DIR, "analysis_subjects.csv"))
)

writeLines(report, file.path(WORK_DIR, "preparation_report.txt"))
cat(paste(report, collapse = "\n"), "\n")

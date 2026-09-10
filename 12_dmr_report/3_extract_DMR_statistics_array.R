#!/usr/bin/env Rscript

# =============================================================================
# 3_extract_DMR_statistics_array.R
#
# Bundled chunk-array worker.
#
# One SLURM array task processes CHUNKS_PER_TASK source methylation chunks
# sequentially.  Each source chunk is still read exactly once by its assigned
# task, and a separate partial RDS is written for every source chunk so the
# existing final collector remains unchanged.
#
# Example with 1733 source chunks and CHUNKS_PER_TASK=10:
#   array task 1   -> manifest rows 1:10
#   array task 2   -> manifest rows 11:20
#   ...
#   array task 174 -> manifest rows 1731:1733
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (!length(args)) stop("Pass SLURM_ARRAY_TASK_ID as argument 1.")

task_id <- suppressWarnings(as.integer(args[1L]))
if (is.na(task_id) || task_id < 1L) stop("Invalid array task ID.")

BASE <- path.expand("~/scratch/UQAC/meth")

OUT_ROOT <- if (length(args) >= 2L) args[2L] else file.path(
    BASE,
    "results/15_revision/9_dmr_report/DMR_statistics"
)

CHUNK_DIR <- if (length(args) >= 3L) args[3L] else file.path(
    BASE,
    "data/meth_split"
)

chunks_per_task <- if (length(args) >= 4L) {
    suppressWarnings(as.integer(args[4L]))
} else {
    10L
}

if (is.na(chunks_per_task) || chunks_per_task < 1L) {
    stop("CHUNKS_PER_TASK must be a positive integer.")
}

WORK_DIR <- file.path(OUT_ROOT, "work")
PARTIAL_DIR <- file.path(OUT_ROOT, "partials")

MANIFEST_FILE <- file.path(WORK_DIR, "chunk_manifest.csv")
MAP_FILE <- file.path(WORK_DIR, "dmr_chunk_map.csv")
DMR_RDS <- file.path(WORK_DIR, "dmr_table.rds")
SUBJECT_FILE <- file.path(WORK_DIR, "analysis_subjects.csv")

for (f in c(MANIFEST_FILE, MAP_FILE, DMR_RDS, SUBJECT_FILE)) {
    if (!file.exists(f)) stop("Required prepared file missing: ", f)
}

manifest <- fread(MANIFEST_FILE)
dmr <- readRDS(DMR_RDS)
chunk_map <- fread(MAP_FILE)
subjects <- fread(SUBJECT_FILE)

# -------------------------------------------------------------------------
# Determine the manifest rows assigned to this bundled array task.
# -------------------------------------------------------------------------

first_row <- (task_id - 1L) * chunks_per_task + 1L
last_row <- min(task_id * chunks_per_task, nrow(manifest))

if (first_row > nrow(manifest)) {
    stop(
        "Array task ", task_id,
        " has no manifest rows. Manifest contains ",
        nrow(manifest), " source chunks."
    )
}

manifest_task <- manifest[first_row:last_row]

cat("============================================================\n")
cat("Bundled DMR methylation chunk task\n")
cat("============================================================\n")
cat("Array task:        ", task_id, "\n", sep = "")
cat("Manifest rows:     ", first_row, "-", last_row, "\n", sep = "")
cat("Source chunks:     ", nrow(manifest_task), "\n", sep = "")
cat("Subjects:          ", nrow(subjects), "\n", sep = "")
cat("Chunk IDs:         ", paste(manifest_task$data_chunk_id, collapse = ", "), "\n\n", sep = "")

# Expected analytic methylation-column names.  The GAM pipeline uses ID_meth.
expected_meth_cols <- paste0(subjects$ID, "_meth")
n_subjects <- nrow(subjects)

# -------------------------------------------------------------------------
# Process one source chunk.
# -------------------------------------------------------------------------

process_one_chunk <- function(chunk_id) {

    chunk_file <- file.path(
        CHUNK_DIR,
        sprintf("chunk_%04d.csv", chunk_id)
    )

    if (!file.exists(chunk_file)) {
        stop("Methylation chunk does not exist: ", chunk_file)
    }

    dmr_ids <- chunk_map[
        data_chunk_id == chunk_id,
        unique(DMR_row_id)
    ]

    if (!length(dmr_ids)) {
        stop(
            "Manifest assigned chunk ", chunk_id,
            " but no DMR rows were mapped."
        )
    }

    dmr_here <- dmr[DMR_row_id %in% dmr_ids]
    setorder(dmr_here, DMR_row_id)

    cat(
        "\n--- chunk_", sprintf("%04d", chunk_id),
        " | DMR rows=", nrow(dmr_here), " ---\n",
        sep = ""
    )

    # ---------------------------------------------------------------------
    # Read only coordinates + analytic *_meth columns that exist.
    # ---------------------------------------------------------------------

    header <- names(fread(
        chunk_file,
        nrows = 0L,
        showProgress = FALSE
    ))

    required_coord <- c("chr", "start", "end")
    missing_coord <- setdiff(required_coord, header)

    if (length(missing_coord)) {
        stop(
            "Chunk ", chunk_id,
            " missing coordinate columns: ",
            paste(missing_coord, collapse = ", ")
        )
    }

    present <- expected_meth_cols %in% header

    cat(
        "Analytic methylation columns present: ",
        sum(present), "/", length(present), "\n", sep = ""
    )

    if (!any(present)) {
        stop(
            "No analytic *_meth columns found in chunk ",
            chunk_id
        )
    }

    meth_cols_present <- expected_meth_cols[present]
    sel_cols <- c(required_coord, meth_cols_present)

    dat <- fread(
        chunk_file,
        select = sel_cols,
        showProgress = FALSE
    )

    dat[, chr := as.character(chr)]
    dat[, start := as.integer(start)]
    dat[, end := as.integer(end)]

    # Keep explicit global subject ordering.
    global_index <- match(
        sub("_meth$", "", meth_cols_present),
        subjects$ID
    )

    if (anyNA(global_index)) {
        stop("Internal subject-column matching error.")
    }

    entries <- vector("list", nrow(dmr_here))

    for (j in seq_len(nrow(dmr_here))) {

        d <- dmr_here[j]

        idx <- which(
            dat$chr == as.character(d$DMR_chr) &
            dat$start >= d$DMR_region_start &
            dat$start <= d$DMR_region_end
        )

        sum_vec <- numeric(n_subjects)
        n_vec <- integer(n_subjects)

        if (!length(idx)) {

            entries[[j]] <- list(
                DMR_row_id = as.integer(d$DMR_row_id),
                cpg_keys = character(0),
                sum_beta = sum_vec,
                n_beta = n_vec
            )

            next
        }

        # Defensive de-duplication within a source chunk.
        cpg_keys <- paste(
            dat$chr[idx],
            dat$start[idx],
            dat$end[idx],
            sep = ":"
        )

        keep <- !duplicated(cpg_keys)
        idx <- idx[keep]
        cpg_keys <- cpg_keys[keep]

        x <- as.matrix(dat[idx, ..meth_cols_present])
        storage.mode(x) <- "double"

        # Validate raw methylation-proportion scale only on values actually
        # entering the DMR.
        finite_x <- x[is.finite(x)]
        if (length(finite_x) &&
            (min(finite_x) < -1e-8 || max(finite_x) > 1 + 1e-8)) {
            stop(
                "*_meth values outside [0,1] in chunk ",
                chunk_id,
                ", DMR_row_id ", d$DMR_row_id
            )
        }

        sum_present <- colSums(x, na.rm = TRUE)
        n_present <- colSums(!is.na(x))

        sum_vec[global_index] <- sum_present
        n_vec[global_index] <- as.integer(n_present)

        entries[[j]] <- list(
            DMR_row_id = as.integer(d$DMR_row_id),
            cpg_keys = cpg_keys,
            sum_beta = sum_vec,
            n_beta = n_vec
        )
    }

    partial <- list(
        bundle_task_id = task_id,
        data_chunk_id = as.integer(chunk_id),
        subject_ids = subjects$ID,
        AA_only = subjects$AA_only,
        DMR_row_ids = as.integer(dmr_here$DMR_row_id),
        entries = entries,
        created = as.character(Sys.time())
    )

    partial_file <- file.path(
        PARTIAL_DIR,
        sprintf("chunk_%04d.rds", chunk_id)
    )

    # Write atomically so the collector never sees a half-written partial.
    tmp_file <- paste0(
        partial_file,
        ".tmp_",
        Sys.getpid()
    )

    saveRDS(
        partial,
        tmp_file,
        compress = "xz"
    )

    if (!file.rename(tmp_file, partial_file)) {
        stop(
            "Could not atomically rename partial result to: ",
            partial_file
        )
    }

    cat("Saved: ", partial_file, "\n", sep = "")

    rm(dat, partial, entries)
    invisible(gc(FALSE))
}

# -------------------------------------------------------------------------
# Sequentially process every source chunk assigned to this array task.
# -------------------------------------------------------------------------

for (chunk_id in manifest_task$data_chunk_id) {
    process_one_chunk(as.integer(chunk_id))
}

cat("\n============================================================\n")
cat("Bundled array task completed successfully\n")
cat("Array task: ", task_id, "\n", sep = "")
cat("Chunks processed: ", nrow(manifest_task), "\n", sep = "")
cat("Completed: ", format(Sys.time()), "\n", sep = "")
cat("============================================================\n")

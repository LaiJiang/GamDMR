#!/usr/bin/env Rscript

# =============================================================================
# 12_prepare_primary_single_cpg_manifest.R
#
# Build an OUTCOME-INDEPENDENT primary single-CpG universe directly from the
# raw/QC methylation chunk files. No p-values, coefficients, AA labels, or
# pre-filtered CpG lists are used.
#
# The script creates:
#   1) one global CpG manifest;
#   2) N_JOBS balanced per-job manifests.
#
# Each per-job manifest contains CpG coordinates + source chunk/row only.
# The downstream primary association runner will evaluate M1, M2, and M3
# separately within each of the 100 OUTER TRAINING-family splits.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

data.table::setDTthreads(1L)

# ------------------------------- configuration -------------------------------

PATH_WK <- path.expand(Sys.getenv("PATH_WK", "~/scratch/UQAC/meth/"))

METH_CHUNK_DIR <- path.expand(Sys.getenv(
  "METH_CHUNK_DIR",
  file.path(PATH_WK, "data/meth_split")
))

MANIFEST_ROOT <- path.expand(Sys.getenv(
  "PRIMARY_CPG_MANIFEST_ROOT",
  file.path(
    PATH_WK,
    "results/15_revision/5_cv/6_single_cpg/0_primary_all_cpg_manifest"
  )
))

N_JOBS <- suppressWarnings(as.integer(Sys.getenv("N_JOBS", "990")))
OVERWRITE_MANIFEST <- Sys.getenv("OVERWRITE_MANIFEST", "0") %in%
  c("1", "TRUE", "true", "T", "yes", "YES")

EXPECTED_CPGS <- suppressWarnings(as.numeric(
  Sys.getenv("EXPECTED_CPGS", "4609564")
))

if (is.na(N_JOBS) || N_JOBS < 1L) stop("N_JOBS must be >= 1.")

dir.create(MANIFEST_ROOT, recursive = TRUE, showWarnings = FALSE)

job_manifest_dir <- file.path(MANIFEST_ROOT, "job_manifests")
dir.create(job_manifest_dir, recursive = TRUE, showWarnings = FALSE)

global_manifest_file <- file.path(
  MANIFEST_ROOT,
  "primary_all_cpg_manifest.tsv"
)

manifest_qc_file <- file.path(
  MANIFEST_ROOT,
  "primary_manifest_qc.tsv"
)

job_counts_file <- file.path(
  MANIFEST_ROOT,
  "primary_manifest_job_counts.tsv"
)

if (!OVERWRITE_MANIFEST) {
  existing_jobs <- list.files(
    job_manifest_dir,
    pattern = "^primary_cpg_manifest_job_[0-9]{4}\\.tsv$",
    full.names = TRUE
  )

  if (file.exists(global_manifest_file) ||
      file.exists(manifest_qc_file) ||
      length(existing_jobs) > 0L) {
    stop(
      "Primary manifest outputs already exist. ",
      "Set OVERWRITE_MANIFEST=1 to recreate them."
    )
  }
}

if (OVERWRITE_MANIFEST) {
  unlink(global_manifest_file, force = TRUE)
  unlink(manifest_qc_file, force = TRUE)
  unlink(job_counts_file, force = TRUE)

  old_jobs <- list.files(
    job_manifest_dir,
    pattern = "^primary_cpg_manifest_job_[0-9]{4}\\.tsv$",
    full.names = TRUE
  )
  if (length(old_jobs) > 0L) unlink(old_jobs, force = TRUE)
}

# -------------------------------- utilities ----------------------------------

normalize_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  sub("\\.0$", "", x)
}

extract_chunk_id <- function(path) {
  b <- basename(path)
  out <- sub("^chunk_([0-9]+)\\.csv$", "\\1", b)
  suppressWarnings(as.integer(out))
}

make_cpg_label <- function(chr, start, end) {
  chr <- as.character(chr)
  start <- suppressWarnings(as.integer(start))
  end <- suppressWarnings(as.integer(end))

  out <- paste0(chr, ":", start)

  has_end <- !is.na(end)
  out[has_end] <- paste0(
    chr[has_end], ":",
    start[has_end], "-",
    end[has_end]
  )

  out
}

# ----------------------------- discover chunks -------------------------------

chunk_files <- list.files(
  METH_CHUNK_DIR,
  pattern = "^chunk_[0-9]+\\.csv$",
  full.names = TRUE
)

if (length(chunk_files) == 0L) {
  stop("No chunk_XXXX.csv files found under: ", METH_CHUNK_DIR)
}

chunk_ids <- vapply(chunk_files, extract_chunk_id, integer(1L))

if (anyNA(chunk_ids)) {
  stop("Could not parse one or more chunk IDs.")
}

ord <- order(chunk_ids)
chunk_files <- chunk_files[ord]
chunk_ids <- chunk_ids[ord]

if (anyDuplicated(chunk_ids)) {
  stop("Duplicated chunk IDs detected.")
}

cat("Primary single-CpG manifest preparation\n")
cat("Methylation chunk directory:", METH_CHUNK_DIR, "\n")
cat("Number of chunk files:", length(chunk_files), "\n")
cat("N_JOBS:", N_JOBS, "\n")

# --------------------------- build global manifest ---------------------------

global_offset <- 0
chunk_qc_list <- vector("list", length(chunk_files))

for (ii in seq_along(chunk_files)) {
  f <- chunk_files[ii]
  chunk_id <- chunk_ids[ii]

  cat(
    sprintf(
      "[%d/%d] Reading CpG coordinates from chunk %d: %s\n",
      ii, length(chunk_files), chunk_id, basename(f)
    )
  )

  hdr <- fread(f, nrows = 0L, showProgress = FALSE)
  hnames <- names(hdr)

  if (!all(c("chr", "start") %in% hnames)) {
    stop("Chunk ", chunk_id, " lacks chr/start columns: ", f)
  }

  select_cols <- c("chr", "start")
  if ("end" %in% hnames) select_cols <- c(select_cols, "end")

  coord <- fread(
    f,
    select = select_cols,
    showProgress = FALSE
  )

  if (!"end" %in% names(coord)) {
    coord[, end := NA_integer_]
  }

  coord[, chr := normalize_chr(chr)]
  coord[, start := suppressWarnings(as.integer(start))]
  coord[, end := suppressWarnings(as.integer(end))]

  coord[, cpg_row_in_chunk := .I]
  coord[, data_chunk_id := as.integer(chunk_id)]

  n_chunk <- nrow(coord)

  if (n_chunk > 0L) {
    coord[, primary_cpg_id := global_offset + .I]
    coord[, CpG := make_cpg_label(chr, start, end)]

    invalid_coord <- is.na(coord$chr) |
      !nzchar(coord$chr) |
      is.na(coord$start)

    chunk_qc_list[[ii]] <- data.table(
      data_chunk_id = as.integer(chunk_id),
      chunk_file = f,
      n_cpgs = n_chunk,
      n_invalid_coordinates = sum(invalid_coord)
    )

    write_dt <- coord[
      ,
      .(
        primary_cpg_id,
        CpG,
        cpg_chr = chr,
        cpg_start = start,
        cpg_end = end,
        data_chunk_id,
        cpg_row_in_chunk
      )
    ]

    file_exists_before <- file.exists(global_manifest_file)

    fwrite(
      write_dt,
      global_manifest_file,
      sep = "\t",
      quote = FALSE,
      na = "NA",
      append = file_exists_before,
      col.names = !file_exists_before
    )

    global_offset <- global_offset + n_chunk
  } else {
    chunk_qc_list[[ii]] <- data.table(
      data_chunk_id = as.integer(chunk_id),
      chunk_file = f,
      n_cpgs = 0L,
      n_invalid_coordinates = 0L
    )
  }

  rm(coord)
  if (ii %% 25L == 0L) invisible(gc(FALSE))
}

n_total <- global_offset

if (n_total == 0L) stop("No CpGs were found.")

cat("\nGlobal primary CpG count:", format(n_total, big.mark = ","), "\n")

# ----------------------- balanced per-job manifests --------------------------

# The global coordinate-only manifest is modest compared with the methylation
# matrix, so read it once and partition deterministically by primary_cpg_id.
all_cpg <- fread(global_manifest_file, showProgress = TRUE)

if (nrow(all_cpg) != n_total) {
  stop(
    "Global manifest row count mismatch: expected ",
    n_total, ", observed ", nrow(all_cpg)
  )
}

# Balanced contiguous partition:
# job = floor((id - 1) * N_JOBS / N) + 1
all_cpg[
  ,
  job_id := pmin(
    N_JOBS,
    floor((primary_cpg_id - 1) * N_JOBS / n_total) + 1L
  )
]

setkey(all_cpg, job_id)

job_counts <- all_cpg[, .N, by = job_id]
setnames(job_counts, "N", "n_cpgs")
setorder(job_counts, job_id)

# Write every job file, including an empty header-only file if N_JOBS exceeds
# the number of CpGs.
manifest_cols <- c(
  "primary_cpg_id",
  "CpG",
  "cpg_chr",
  "cpg_start",
  "cpg_end",
  "data_chunk_id",
  "cpg_row_in_chunk"
)

empty_manifest <- all_cpg[0, ..manifest_cols]

for (j in seq_len(N_JOBS)) {
  out_file <- file.path(
    job_manifest_dir,
    sprintf("primary_cpg_manifest_job_%04d.tsv", j)
  )

  dtj <- all_cpg[.(j), ..manifest_cols]

  if (nrow(dtj) == 0L) dtj <- empty_manifest

  fwrite(
    dtj,
    out_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  if (j %% 100L == 0L || j == N_JOBS) {
    cat("Wrote job manifests:", j, "/", N_JOBS, "\n")
  }
}

fwrite(
  job_counts,
  job_counts_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

chunk_qc <- rbindlist(chunk_qc_list, use.names = TRUE, fill = TRUE)

manifest_qc <- data.table(
  n_chunk_files = length(chunk_files),
  min_chunk_id = min(chunk_ids),
  max_chunk_id = max(chunk_ids),
  n_primary_cpgs = n_total,
  expected_cpgs_reference = EXPECTED_CPGS,
  difference_from_expected = if (is.finite(EXPECTED_CPGS)) {
    n_total - EXPECTED_CPGS
  } else {
    NA_real_
  },
  n_invalid_coordinates = sum(chunk_qc$n_invalid_coordinates),
  n_jobs = N_JOBS,
  min_cpgs_per_job = min(job_counts$n_cpgs),
  max_cpgs_per_job = max(job_counts$n_cpgs),
  mean_cpgs_per_job = mean(job_counts$n_cpgs),
  global_manifest_file = global_manifest_file,
  job_manifest_dir = job_manifest_dir,
  created_at = format(
    Sys.time(),
    tz = "America/Toronto",
    usetz = TRUE
  )
)

fwrite(
  manifest_qc,
  manifest_qc_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  chunk_qc,
  file.path(MANIFEST_ROOT, "primary_manifest_chunk_qc.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\nManifest preparation completed.\n")
print(manifest_qc)

if (is.finite(EXPECTED_CPGS) && n_total != EXPECTED_CPGS) {
  cat(
    "\nNOTE: observed CpG count differs from EXPECTED_CPGS by ",
    n_total - EXPECTED_CPGS,
    ". This is a QC note only; no outcome-based filtering was applied.\n",
    sep = ""
  )
}

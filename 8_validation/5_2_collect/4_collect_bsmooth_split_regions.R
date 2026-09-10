#!/usr/bin/env Rscript

# Collect BSmooth family-training results across all task-specific files.
# Per-split files reproduce the previous 6_BSmooth_DMRs.csv logic:
# merge summary files and retain rows with n_dmrs >= 1.

suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)

PATH_wk <- path.expand("~/scratch/UQAC/meth/")
input_dir <- file.path(
  PATH_wk, "results/15_revision/5_cv/4_bsmooth/training_splits"
)
output_root <- file.path(
  PATH_wk, "results/15_revision/5_cv/4_bsmooth/collected_by_split"
)
summary_split_dir <- file.path(output_root, "6_BSmooth_DMRs_by_split")
detail_split_dir <- file.path(output_root, "6_BSmooth_DMR_details_by_split")
qc_dir <- file.path(output_root, "qc")

for (d in c(summary_split_dir, detail_split_dir, qc_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

env_int <- function(name, default) {
  x <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (is.na(x)) stop(name, " must be an integer.")
  x
}
env_bool <- function(name, default = FALSE) {
  Sys.getenv(name, if (default) "1" else "0") %in%
    c("1", "TRUE", "true", "T", "yes", "YES")
}

split_start <- env_int("SPLIT_START", 1L)
split_end <- env_int("SPLIT_END", 100L)
expected_jobs <- env_int("EXPECTED_JOBS", 990L)
overwrite <- env_bool("OVERWRITE", FALSE)

if (split_start < 1L || split_end < split_start) {
  stop("Invalid SPLIT_START/SPLIT_END.")
}
if (expected_jobs < 1L) stop("EXPECTED_JOBS must be positive.")

split_ids <- seq.int(split_start, split_end)

combined_summary_file <- file.path(output_root, "6_BSmooth_DMRs_all_splits.csv")
combined_detail_file <- file.path(output_root, "6_BSmooth_DMR_details_all_splits.tsv")
all_summary_file <- file.path(output_root, "bsmooth_training_summary_all_rows.tsv")
all_detail_file <- file.path(output_root, "bsmooth_training_dmr_details_all_rows.tsv")
status_file <- file.path(qc_dir, "bsmooth_status_counts_by_split.tsv")
split_qc_file <- file.path(qc_dir, "bsmooth_collection_qc_by_split.tsv")
duplicate_file <- file.path(qc_dir, "bsmooth_duplicate_region_split_rows.tsv")
missing_file <- file.path(qc_dir, "bsmooth_missing_job_files.tsv")
count_check_file <- file.path(qc_dir, "bsmooth_summary_detail_count_check.tsv")
summary_manifest_file <- file.path(qc_dir, "bsmooth_summary_input_manifest.tsv")
detail_manifest_file <- file.path(qc_dir, "bsmooth_detail_input_manifest.tsv")

old_outputs <- c(
  combined_summary_file, combined_detail_file, all_summary_file, all_detail_file,
  status_file, split_qc_file, duplicate_file, missing_file, count_check_file,
  summary_manifest_file, detail_manifest_file,
  list.files(summary_split_dir, full.names = TRUE),
  list.files(detail_split_dir, full.names = TRUE)
)
old_outputs <- old_outputs[file.exists(old_outputs)]

if (length(old_outputs) && !overwrite) {
  stop(
    "Collector output already exists. Set OVERWRITE=1 to replace it. ",
    "First file: ", old_outputs[1L]
  )
}
if (overwrite && length(old_outputs)) unlink(old_outputs, force = TRUE)

if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)

summary_files <- list.files(
  input_dir,
  pattern = "^bsmooth_training_summary_job_[0-9]{4}\\.tsv$",
  full.names = TRUE
)
detail_files <- list.files(
  input_dir,
  pattern = "^bsmooth_training_dmrs_job_[0-9]{4}\\.tsv$",
  full.names = TRUE
)

if (!length(summary_files)) stop("No BSmooth summary files found.")

extract_job <- function(paths, prefix) {
  as.integer(sub(
    paste0("^", prefix, "_job_([0-9]{4})\\.tsv$"),
    "\\1", basename(paths)
  ))
}

summary_jobs <- extract_job(summary_files, "bsmooth_training_summary")
detail_jobs <- extract_job(detail_files, "bsmooth_training_dmrs")

summary_manifest <- data.table(
  source_file = basename(summary_files),
  source_path = summary_files,
  job_id = summary_jobs,
  file_size_bytes = file.info(summary_files)$size
)
detail_manifest <- data.table(
  source_file = basename(detail_files),
  source_path = detail_files,
  job_id = detail_jobs,
  file_size_bytes = file.info(detail_files)$size
)
setorder(summary_manifest, job_id)
setorder(detail_manifest, job_id)
fwrite(summary_manifest, summary_manifest_file, sep = "\t", quote = FALSE)
fwrite(detail_manifest, detail_manifest_file, sep = "\t", quote = FALSE)

missing_summary_jobs <- setdiff(seq_len(expected_jobs), summary_jobs)
missing_detail_jobs <- setdiff(seq_len(expected_jobs), detail_jobs)

fwrite(
  rbindlist(list(
    data.table(file_type = "summary", job_id = missing_summary_jobs),
    data.table(file_type = "detail", job_id = missing_detail_jobs)
  )),
  missing_file, sep = "\t", quote = FALSE
)

cat("Summary files found:", length(summary_files), "\n")
cat("Detailed DMR files found:", length(detail_files), "\n")
cat("Expected jobs:", expected_jobs, "\n")

if (length(missing_summary_jobs)) {
  warning(
    "Missing ", length(missing_summary_jobs), " summary files. First IDs: ",
    paste(head(missing_summary_jobs, 20L), collapse = ", ")
  )
}

summary_columns <- c(
  "prep_job_id", "region_id", "data_chunk_id", "chr",
  "region_start", "region_end", "splitID", "status", "message",
  "n_train_FIDs_requested", "n_train_FIDs_observed",
  "n_train_samples", "n_train_cases", "n_train_controls",
  "n_cpgs", "n_finite_stats", "tstat_column",
  "abs_t_q95", "max_abs_t", "n_dmrs", "fit_seconds"
)

read_summary <- function(f) {
  header <- names(fread(f, nrows = 0L, showProgress = FALSE))
  missing <- setdiff(summary_columns, header)
  if (length(missing)) {
    stop(basename(f), " is missing: ", paste(missing, collapse = ", "))
  }
  x <- fread(f, select = summary_columns, showProgress = FALSE)
  x[, source_file := basename(f)]
  x
}

summary_list <- vector("list", length(summary_files))
for (i in seq_along(summary_files)) {
  if (i == 1L || i == length(summary_files) || i %% 50L == 0L) {
    cat("Reading summary file", i, "of", length(summary_files), "\n")
  }
  summary_list[[i]] <- read_summary(summary_files[i])
}
summary_all <- rbindlist(summary_list, use.names = TRUE, fill = TRUE)
rm(summary_list)
gc(FALSE)

summary_all[, `:=`(
  splitID = as.integer(splitID),
  region_id = as.integer(region_id),
  n_dmrs = as.integer(n_dmrs),
  chr = as.character(chr),
  region_start = as.integer(region_start),
  region_end = as.integer(region_end)
)]
summary_all <- summary_all[
  !is.na(splitID) & splitID >= split_start & splitID <= split_end
]
if (!nrow(summary_all)) stop("No summary rows remain.")

# Resolve resumed/retried duplicate region-split rows.
summary_all[, key := paste(region_id, splitID, sep = "::")]
dups <- summary_all[duplicated(key) | duplicated(key, fromLast = TRUE)]

if (nrow(dups)) {
  fwrite(
    dups[, .(
      n_rows = .N,
      statuses = paste(sort(unique(status)), collapse = " | "),
      n_dmrs_values = paste(sort(unique(n_dmrs)), collapse = " | "),
      source_files = paste(sort(unique(source_file)), collapse = " | ")
    ), by = .(region_id, splitID)],
    duplicate_file, sep = "\t", quote = FALSE
  )

  summary_all[, priority :=
    1000L * as.integer(status == "FIT_OK") +
    100L * as.integer(!is.na(n_dmrs)) +
    pmax(fifelse(is.na(n_dmrs), 0L, n_dmrs), 0L)
  ]
  setorder(summary_all, key, -priority, source_file)
  summary_all <- summary_all[!duplicated(key)]
  summary_all[, priority := NULL]
} else {
  fwrite(
    data.table(
      region_id = integer(), splitID = integer(), n_rows = integer(),
      statuses = character(), n_dmrs_values = character(),
      source_files = character()
    ),
    duplicate_file, sep = "\t", quote = FALSE
  )
}

setorder(summary_all, splitID, region_id)
fwrite(summary_all, all_summary_file, sep = "\t", quote = FALSE, na = "NA")

# Collect actual dmrFinder output rows.
detail_all <- data.table()

if (length(detail_files)) {
  read_detail <- function(f) {
    x <- fread(f, showProgress = FALSE)
    missing <- setdiff(c("region_id", "splitID"), names(x))
    if (length(missing)) {
      stop(basename(f), " is missing: ", paste(missing, collapse = ", "))
    }
    x[, source_file := basename(f)]
    x
  }

  detail_list <- vector("list", length(detail_files))
  for (i in seq_along(detail_files)) {
    if (i == 1L || i == length(detail_files) || i %% 50L == 0L) {
      cat("Reading detailed DMR file", i, "of", length(detail_files), "\n")
    }
    detail_list[[i]] <- read_detail(detail_files[i])
  }

  detail_all <- rbindlist(detail_list, use.names = TRUE, fill = TRUE)
  rm(detail_list)
  gc(FALSE)

  detail_all[, `:=`(
    splitID = as.integer(splitID),
    region_id = as.integer(region_id)
  )]
  detail_all <- detail_all[
    !is.na(splitID) & splitID >= split_start & splitID <= split_end
  ]

  if (nrow(detail_all)) {
    ord <- intersect(
      c("splitID", "region_id", "chr", "start", "end"),
      names(detail_all)
    )
    setorderv(detail_all, ord)
    detail_all[, dmr_index_within_region := seq_len(.N),
               by = .(splitID, region_id)]
    fwrite(detail_all, all_detail_file, sep = "\t", quote = FALSE, na = "NA")
  }
}

summary_selected_list <- vector("list", length(split_ids))
detail_selected_list <- vector("list", length(split_ids))
qc_list <- vector("list", length(split_ids))
names(summary_selected_list) <- names(detail_selected_list) <-
  names(qc_list) <- as.character(split_ids)

for (sid in split_ids) {
  cat("Writing split", sid, "\n")

  split_all <- summary_all[splitID == sid]
  # Same selection rule as the old 6_col_BSmooth.R.
  selected <- split_all[!is.na(n_dmrs) & n_dmrs >= 1L]
  setorder(selected, region_id)

  selected_out <- copy(selected)
  selected_out[, c("key", "source_file") := NULL]

  fwrite(
    selected_out,
    file.path(
      summary_split_dir,
      sprintf("6_BSmooth_DMRs_split_%03d.csv", sid)
    ),
    sep = ",", quote = FALSE, na = "NA"
  )
  summary_selected_list[[as.character(sid)]] <- selected_out

  if (nrow(detail_all)) {
    split_detail <- detail_all[splitID == sid]
  } else {
    split_detail <- data.table(region_id = integer(), splitID = integer())
  }

  fwrite(
    split_detail,
    file.path(
      detail_split_dir,
      sprintf("6_BSmooth_DMR_details_split_%03d.tsv", sid)
    ),
    sep = "\t", quote = FALSE, na = "NA"
  )
  detail_selected_list[[as.character(sid)]] <- split_detail

  qc_list[[as.character(sid)]] <- data.table(
    splitID = sid,
    n_summary_rows = nrow(split_all),
    n_unique_regions = uniqueN(split_all$region_id),
    n_fit_ok = split_all[status == "FIT_OK", .N],
    n_parent_regions_with_dmrs = nrow(selected),
    n_reported_dmrs_from_summary = sum(selected$n_dmrs, na.rm = TRUE),
    n_detailed_dmr_rows = nrow(split_detail)
  )
}

summary_selected_all <- rbindlist(
  summary_selected_list, use.names = TRUE, fill = TRUE
)
setorder(summary_selected_all, splitID, region_id)
fwrite(
  summary_selected_all,
  combined_summary_file,
  sep = ",", quote = FALSE, na = "NA"
)

detail_selected_all <- rbindlist(
  detail_selected_list, use.names = TRUE, fill = TRUE
)
if (nrow(detail_selected_all) &&
    "dmr_index_within_region" %in% names(detail_selected_all)) {
  setorder(
    detail_selected_all, splitID, region_id, dmr_index_within_region
  )
}
fwrite(
  detail_selected_all,
  combined_detail_file,
  sep = "\t", quote = FALSE, na = "NA"
)

status_counts <- summary_all[, .N, by = .(splitID, status)]
setorder(status_counts, splitID, status)
fwrite(status_counts, status_file, sep = "\t", quote = FALSE)

split_qc <- rbindlist(qc_list, use.names = TRUE, fill = TRUE)
setorder(split_qc, splitID)
fwrite(split_qc, split_qc_file, sep = "\t", quote = FALSE)

count_check <- split_qc[, .(
  splitID,
  n_parent_regions_with_dmrs,
  n_reported_dmrs_from_summary,
  n_detailed_dmr_rows,
  detail_count_matches_summary =
    n_reported_dmrs_from_summary == n_detailed_dmr_rows
)]
fwrite(count_check, count_check_file, sep = "\t", quote = FALSE)

cat("\nBSmooth collection finished.\n")
cat("Per-split summary directory:", summary_split_dir, "\n")
cat("Combined summary file:", combined_summary_file, "\n")
cat("Per-split detail directory:", detail_split_dir, "\n")
cat("Combined detailed DMR file:", combined_detail_file, "\n")
cat("All summary rows:", all_summary_file, "\n")
if (file.exists(all_detail_file)) cat("All detail rows:", all_detail_file, "\n")
cat("QC directory:", qc_dir, "\n")
cat("Collected summary rows:", nrow(summary_all), "\n")
cat("Parent region-split rows with n_dmrs >= 1:",
    nrow(summary_selected_all), "\n")
cat("Detailed DMR rows:", nrow(detail_selected_all), "\n")
cat("Missing summary files:", length(missing_summary_jobs), "\n")

print(split_qc[, .(
  n_splits = .N,
  median_summary_rows = median(n_summary_rows),
  median_parent_regions_with_dmrs =
    median(n_parent_regions_with_dmrs),
  min_parent_regions_with_dmrs =
    min(n_parent_regions_with_dmrs),
  max_parent_regions_with_dmrs =
    max(n_parent_regions_with_dmrs),
  total_detailed_dmrs = sum(n_detailed_dmr_rows)
)])

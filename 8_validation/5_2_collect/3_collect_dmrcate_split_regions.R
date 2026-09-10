#!/usr/bin/env Rscript

# Collect DMRcate family-training results.
#
# Main outputs:
#   - one 7_DMRcate_DMRs_split_XXX.csv for each split
#   - one 7_DMRcate_DMRs_all_splits.csv binding all splits
#
# The per-split files reproduce the previous collector logic:
# merge summary files and retain rows with n_dmrs >= 1.
#
# Detailed DMR rows from dmrcate_training_dmrs_job_*.tsv are also collected
# for auditing and downstream use.

suppressPackageStartupMessages({
  library(data.table)
})

data.table::setDTthreads(1L)

PATH_wk <- path.expand("~/scratch/UQAC/meth/")

input_dir <- file.path(
  PATH_wk,
  "results/15_revision/5_cv/3_dmrcate/training_splits"
)

output_root <- file.path(
  PATH_wk,
  "results/15_revision/5_cv/3_dmrcate/collected_by_split"
)

summary_split_dir <- file.path(
  output_root,
  "7_DMRcate_DMRs_by_split"
)

detail_split_dir <- file.path(
  output_root,
  "7_DMRcate_DMR_details_by_split"
)

qc_dir <- file.path(output_root, "qc")

dir.create(summary_split_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(detail_split_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

split_start <- as.integer(Sys.getenv("SPLIT_START", "1"))
split_end <- as.integer(Sys.getenv("SPLIT_END", "100"))
expected_jobs <- as.integer(Sys.getenv("EXPECTED_JOBS", "900"))
overwrite <- identical(Sys.getenv("OVERWRITE", "0"), "1")

if (
  is.na(split_start) ||
    is.na(split_end) ||
    split_start < 1L ||
    split_end < split_start
) {
  stop("Invalid SPLIT_START/SPLIT_END values.")
}

if (is.na(expected_jobs) || expected_jobs < 1L) {
  stop("EXPECTED_JOBS must be a positive integer.")
}

split_ids <- seq.int(split_start, split_end)

combined_summary_file <- file.path(
  output_root,
  "7_DMRcate_DMRs_all_splits.csv"
)

combined_detail_file <- file.path(
  output_root,
  "7_DMRcate_DMR_details_all_splits.tsv"
)

all_summary_file <- file.path(
  output_root,
  "dmrcate_training_summary_all_rows.tsv"
)

all_detail_file <- file.path(
  output_root,
  "dmrcate_training_dmr_details_all_rows.tsv"
)

status_qc_file <- file.path(
  qc_dir,
  "dmrcate_status_counts_by_split.tsv"
)

split_qc_file <- file.path(
  qc_dir,
  "dmrcate_collection_qc_by_split.tsv"
)

duplicate_qc_file <- file.path(
  qc_dir,
  "dmrcate_duplicate_region_split_rows.tsv"
)

missing_jobs_file <- file.path(
  qc_dir,
  "dmrcate_missing_jobs.tsv"
)

summary_detail_check_file <- file.path(
  qc_dir,
  "dmrcate_summary_detail_count_check.tsv"
)

existing_outputs <- c(
  combined_summary_file,
  combined_detail_file,
  all_summary_file,
  all_detail_file,
  status_qc_file,
  split_qc_file,
  duplicate_qc_file,
  missing_jobs_file,
  summary_detail_check_file,
  list.files(
    summary_split_dir,
    pattern = "^7_DMRcate_DMRs_split_[0-9]{3}\\.csv$",
    full.names = TRUE
  ),
  list.files(
    detail_split_dir,
    pattern = "^7_DMRcate_DMR_details_split_[0-9]{3}\\.tsv$",
    full.names = TRUE
  )
)

existing_outputs <- existing_outputs[file.exists(existing_outputs)]

if (length(existing_outputs) > 0L && !overwrite) {
  stop(
    "Collector output already exists. Set OVERWRITE=1 to replace it. ",
    "First existing file: ",
    existing_outputs[1L]
  )
}

if (overwrite && length(existing_outputs) > 0L) {
  unlink(existing_outputs, force = TRUE)
}

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir)
}

summary_files <- list.files(
  input_dir,
  pattern = "^dmrcate_training_summary_job_[0-9]{4}\\.tsv$",
  full.names = TRUE
)

detail_files <- list.files(
  input_dir,
  pattern = "^dmrcate_training_dmrs_job_[0-9]{4}\\.tsv$",
  full.names = TRUE
)

if (length(summary_files) == 0L) {
  stop("No DMRcate training summary files found in: ", input_dir)
}

extract_job_id <- function(paths, prefix) {
  as.integer(
    sub(
      paste0("^", prefix, "_job_([0-9]{4})\\.tsv$"),
      "\\1",
      basename(paths)
    )
  )
}

summary_job_ids <- extract_job_id(
  summary_files,
  "dmrcate_training_summary"
)

detail_job_ids <- extract_job_id(
  detail_files,
  "dmrcate_training_dmrs"
)

missing_summary_jobs <- setdiff(
  seq_len(expected_jobs),
  summary_job_ids
)

missing_detail_jobs <- setdiff(
  seq_len(expected_jobs),
  detail_job_ids
)

missing_jobs_dt <- rbindlist(
  list(
    data.table(
      file_type = "summary",
      job_id = missing_summary_jobs
    ),
    data.table(
      file_type = "detail",
      job_id = missing_detail_jobs
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  missing_jobs_dt,
  missing_jobs_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Summary files found:", length(summary_files), "\n")
cat("Detailed DMR files found:", length(detail_files), "\n")
cat("Expected jobs:", expected_jobs, "\n")

if (length(missing_summary_jobs) > 0L) {
  warning(
    "Missing ",
    length(missing_summary_jobs),
    " summary files. First missing job IDs: ",
    paste(head(missing_summary_jobs, 20L), collapse = ", ")
  )
}

required_summary_columns <- c(
  "prep_job_id",
  "region_id",
  "data_chunk_id",
  "chr",
  "region_start",
  "region_end",
  "splitID",
  "status",
  "message",
  "n_train_FIDs_requested",
  "n_train_FIDs_observed",
  "n_train_samples",
  "n_train_cases",
  "n_train_controls",
  "n_cpgs_before_coverage",
  "n_cpgs_after_coverage",
  "n_seed_cpgs",
  "n_dmrs",
  "fit_seconds"
)

read_summary_file <- function(file_path) {
  header <- names(
    fread(file_path, nrows = 0L, showProgress = FALSE)
  )

  missing_columns <- setdiff(
    required_summary_columns,
    header
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Summary file is missing required columns: ",
      basename(file_path),
      " [",
      paste(missing_columns, collapse = ", "),
      "]"
    )
  }

  dt <- fread(
    file_path,
    select = required_summary_columns,
    showProgress = FALSE
  )

  dt[, source_file := basename(file_path)]
  dt
}

summary_list <- vector("list", length(summary_files))

for (i in seq_along(summary_files)) {
  if (i == 1L || i == length(summary_files) || i %% 50L == 0L) {
    cat("Reading summary file", i, "of", length(summary_files), "\n")
  }

  summary_list[[i]] <- read_summary_file(summary_files[i])
}

summary_all <- rbindlist(
  summary_list,
  use.names = TRUE,
  fill = TRUE
)

rm(summary_list)
invisible(gc(FALSE))

if (nrow(summary_all) == 0L) {
  stop("The collected DMRcate summary table is empty.")
}

summary_all[, splitID := as.integer(splitID)]
summary_all[, region_id := as.integer(region_id)]
summary_all[, n_dmrs := as.integer(n_dmrs)]
summary_all[, chr := as.character(chr)]
summary_all[, region_start := as.integer(region_start)]
summary_all[, region_end := as.integer(region_end)]

summary_all <- summary_all[
  !is.na(splitID) &
    splitID >= split_start &
    splitID <= split_end
]

if (nrow(summary_all) == 0L) {
  stop("No summary rows remain in the requested split range.")
}


# ---------------------------- duplicate handling -----------------------------

summary_all[
  ,
  region_split_key := paste(region_id, splitID, sep = "::")
]

duplicate_rows <- summary_all[
  duplicated(region_split_key) |
    duplicated(region_split_key, fromLast = TRUE)
]

if (nrow(duplicate_rows) > 0L) {
  duplicate_qc <- duplicate_rows[
    ,
    .(
      n_rows = .N,
      statuses = paste(sort(unique(status)), collapse = " | "),
      n_dmrs_values = paste(sort(unique(n_dmrs)), collapse = " | "),
      source_files = paste(sort(unique(source_file)), collapse = " | ")
    ),
    by = .(
      region_id,
      splitID
    )
  ]

  fwrite(
    duplicate_qc,
    duplicate_qc_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  summary_all[
    ,
    duplicate_priority := (
      1000L * as.integer(status == "FIT_OK") +
        100L * as.integer(!is.na(n_dmrs)) +
        pmax(fifelse(is.na(n_dmrs), 0L, n_dmrs), 0L)
    )
  ]

  setorder(
    summary_all,
    region_split_key,
    -duplicate_priority,
    source_file
  )

  summary_all <- summary_all[
    !duplicated(region_split_key)
  ]

  summary_all[, duplicate_priority := NULL]
} else {
  fwrite(
    data.table(
      region_id = integer(),
      splitID = integer(),
      n_rows = integer(),
      statuses = character(),
      n_dmrs_values = character(),
      source_files = character()
    ),
    duplicate_qc_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

setorder(summary_all, splitID, region_id)

fwrite(
  summary_all,
  all_summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ----------------------------- detailed DMR rows ------------------------------

detail_all <- data.table()

if (length(detail_files) > 0L) {
  read_detail_file <- function(file_path) {
    dt <- fread(file_path, showProgress = FALSE)

    required_detail <- c(
      "region_id",
      "splitID"
    )

    missing_columns <- setdiff(
      required_detail,
      names(dt)
    )

    if (length(missing_columns) > 0L) {
      stop(
        "Detailed DMR file is missing required columns: ",
        basename(file_path),
        " [",
        paste(missing_columns, collapse = ", "),
        "]"
      )
    }

    dt[, source_file := basename(file_path)]
    dt
  }

  detail_list <- vector("list", length(detail_files))

  for (i in seq_along(detail_files)) {
    if (i == 1L || i == length(detail_files) || i %% 50L == 0L) {
      cat("Reading detailed DMR file", i, "of", length(detail_files), "\n")
    }

    detail_list[[i]] <- read_detail_file(detail_files[i])
  }

  detail_all <- rbindlist(
    detail_list,
    use.names = TRUE,
    fill = TRUE
  )

  rm(detail_list)
  invisible(gc(FALSE))

  detail_all[, splitID := as.integer(splitID)]
  detail_all[, region_id := as.integer(region_id)]

  detail_all <- detail_all[
    !is.na(splitID) &
      splitID >= split_start &
      splitID <= split_end
  ]

  if (nrow(detail_all) > 0L) {
    order_cols <- intersect(
      c("splitID", "region_id", "seqnames", "start", "end"),
      names(detail_all)
    )

    setorderv(detail_all, order_cols)

    detail_all[
      ,
      dmr_index_within_region := seq_len(.N),
      by = .(
        splitID,
        region_id
      )
    ]

    fwrite(
      detail_all,
      all_detail_file,
      sep = "\t",
      quote = FALSE,
      na = "NA"
    )
  }
}

# -------------------------- split-specific outputs ----------------------------

selected_summary_list <- vector(
  "list",
  length(split_ids)
)

selected_detail_list <- vector(
  "list",
  length(split_ids)
)

split_qc_list <- vector(
  "list",
  length(split_ids)
)

names(selected_summary_list) <- as.character(split_ids)
names(selected_detail_list) <- as.character(split_ids)
names(split_qc_list) <- as.character(split_ids)

for (split_id in split_ids) {
  cat("Writing split", split_id, "\n")

  split_all <- summary_all[
    splitID == split_id
  ]

  # Match the previous 7_col_DMRcate.R rule.
  split_selected <- split_all[
    !is.na(n_dmrs) &
      n_dmrs >= 1L
  ]

  setorder(split_selected, region_id)

  split_summary_file <- file.path(
    summary_split_dir,
    sprintf(
      "7_DMRcate_DMRs_split_%03d.csv",
      split_id
    )
  )

  split_selected_output <- copy(split_selected)

  split_selected_output[
    ,
    c(
      "region_split_key",
      "source_file"
    ) := NULL
  ]

  fwrite(
    split_selected_output,
    split_summary_file,
    sep = ",",
    quote = FALSE,
    na = "NA"
  )

  selected_summary_list[[as.character(split_id)]] <-
    split_selected_output

  if (nrow(detail_all) > 0L) {
    split_detail <- detail_all[
      splitID == split_id
    ]
  } else {
    split_detail <- data.table(
      region_id = integer(),
      splitID = integer()
    )
  }

  split_detail_file <- file.path(
    detail_split_dir,
    sprintf(
      "7_DMRcate_DMR_details_split_%03d.tsv",
      split_id
    )
  )

  fwrite(
    split_detail,
    split_detail_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  selected_detail_list[[as.character(split_id)]] <- split_detail

  split_qc_list[[as.character(split_id)]] <- data.table(
    splitID = split_id,
    n_summary_rows = nrow(split_all),
    n_unique_regions = uniqueN(split_all$region_id),
    n_fit_ok = split_all[status == "FIT_OK", .N],
    n_parent_regions_with_dmrs = nrow(split_selected),
    n_reported_dmrs_from_summary = sum(
      split_selected$n_dmrs,
      na.rm = TRUE
    ),
    n_detailed_dmr_rows = nrow(split_detail)
  )
}


# ----------------------------- combined outputs ------------------------------

selected_summary_all <- rbindlist(
  selected_summary_list,
  use.names = TRUE,
  fill = TRUE
)

setorder(
  selected_summary_all,
  splitID,
  region_id
)

fwrite(
  selected_summary_all,
  combined_summary_file,
  sep = ",",
  quote = FALSE,
  na = "NA"
)

selected_detail_all <- rbindlist(
  selected_detail_list,
  use.names = TRUE,
  fill = TRUE
)

if (
  nrow(selected_detail_all) > 0L &&
    "dmr_index_within_region" %in% names(selected_detail_all)
) {
  setorder(
    selected_detail_all,
    splitID,
    region_id,
    dmr_index_within_region
  )
}

fwrite(
  selected_detail_all,
  combined_detail_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ----------------------------------- QC --------------------------------------

status_counts <- summary_all[
  ,
  .N,
  by = .(
    splitID,
    status
  )
]

setorder(
  status_counts,
  splitID,
  status
)

fwrite(
  status_counts,
  status_qc_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

split_qc <- rbindlist(
  split_qc_list,
  use.names = TRUE,
  fill = TRUE
)

setorder(split_qc, splitID)

fwrite(
  split_qc,
  split_qc_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

summary_detail_check <- split_qc[
  ,
  .(
    splitID,
    n_parent_regions_with_dmrs,
    n_reported_dmrs_from_summary,
    n_detailed_dmr_rows,
    detail_count_matches_summary = (
      n_reported_dmrs_from_summary ==
        n_detailed_dmr_rows
    )
  )
]

fwrite(
  summary_detail_check,
  summary_detail_check_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------- final report --------------------------------

cat("\nDMRcate collection finished.\n")
cat("Per-split summary directory:", summary_split_dir, "\n")
cat("Combined summary file:", combined_summary_file, "\n")
cat("Per-split detail directory:", detail_split_dir, "\n")
cat("Combined detailed DMR file:", combined_detail_file, "\n")
cat("All summary rows:", all_summary_file, "\n")

if (file.exists(all_detail_file)) {
  cat("All detailed DMR rows:", all_detail_file, "\n")
}

cat("QC directory:", qc_dir, "\n")
cat("Collected summary rows:", nrow(summary_all), "\n")
cat(
  "Parent region-split rows with n_dmrs >= 1:",
  nrow(selected_summary_all),
  "\n"
)
cat("Detailed DMR rows:", nrow(selected_detail_all), "\n")
cat("Missing summary job files:", length(missing_summary_jobs), "\n")

print(
  split_qc[
    ,
    .(
      n_splits = .N,
      median_summary_rows = median(n_summary_rows),
      median_parent_regions_with_dmrs = median(
        n_parent_regions_with_dmrs
      ),
      min_parent_regions_with_dmrs = min(
        n_parent_regions_with_dmrs
      ),
      max_parent_regions_with_dmrs = max(
        n_parent_regions_with_dmrs
      ),
      total_detailed_dmrs = sum(n_detailed_dmr_rows)
    )
  ]
)

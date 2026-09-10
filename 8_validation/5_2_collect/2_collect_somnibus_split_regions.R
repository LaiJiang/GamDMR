#!/usr/bin/env Rscript

# Collect SOMNiBUS family-training results across all bundle tasks.
# Produces one legacy-format region file per split and one combined file.

suppressPackageStartupMessages({
  library(data.table)
})

data.table::setDTthreads(1L)

PATH_wk <- path.expand("~/scratch/UQAC/meth/")

input_dir <- file.path(
  PATH_wk,
  "results/15_revision/5_cv/2_somnibus/training_splits"
)

output_root <- file.path(
  PATH_wk,
  "results/15_revision/5_cv/2_somnibus/collected_by_split"
)

split_output_dir <- file.path(
  output_root,
  "16_somnibus_regions_by_split"
)

qc_dir <- file.path(output_root, "qc")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(split_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

split_start <- suppressWarnings(
  as.integer(Sys.getenv("SPLIT_START", "1"))
)

split_end <- suppressWarnings(
  as.integer(Sys.getenv("SPLIT_END", "100"))
)

region_p_threshold <- suppressWarnings(
  as.numeric(Sys.getenv("SOMNIBUS_REGION_P_THRESHOLD", "0.01"))
)

expected_bundles <- suppressWarnings(
  as.integer(Sys.getenv("EXPECTED_BUNDLES", "900"))
)

overwrite <- identical(Sys.getenv("OVERWRITE", "0"), "1")

if (
  is.na(split_start) ||
    is.na(split_end) ||
    split_start < 1L ||
    split_end < split_start
) {
  stop("Invalid SPLIT_START/SPLIT_END values.")
}

if (
  !is.finite(region_p_threshold) ||
    region_p_threshold <= 0 ||
    region_p_threshold >= 1
) {
  stop("SOMNIBUS_REGION_P_THRESHOLD must be in (0, 1).")
}

if (is.na(expected_bundles) || expected_bundles < 1L) {
  stop("EXPECTED_BUNDLES must be a positive integer.")
}

split_ids_expected <- seq.int(split_start, split_end)

all_adjusted_file <- file.path(
  output_root,
  "somnibus_adjusted_results_all_splits.tsv"
)

all_regions_file <- file.path(
  output_root,
  "16_somnibus_regions_all_splits.csv"
)

all_regions_detailed_file <- file.path(
  output_root,
  "16_somnibus_regions_all_splits_detailed.tsv"
)

split_qc_file <- file.path(
  qc_dir,
  "somnibus_collection_qc_by_split.tsv"
)

status_qc_file <- file.path(
  qc_dir,
  "somnibus_status_counts_by_split.tsv"
)

input_manifest_file <- file.path(
  qc_dir,
  "somnibus_input_file_manifest.tsv"
)

duplicate_qc_file <- file.path(
  qc_dir,
  "somnibus_duplicate_region_split_rows.tsv"
)

split_files_existing <- list.files(
  split_output_dir,
  pattern = "^16_somnibus_regions_split_[0-9]{3}\\.csv$",
  full.names = TRUE
)

if (!overwrite) {
  existing_outputs <- c(
    all_adjusted_file,
    all_regions_file,
    all_regions_detailed_file,
    split_qc_file,
    status_qc_file,
    input_manifest_file,
    duplicate_qc_file,
    split_files_existing
  )

  existing_outputs <- existing_outputs[file.exists(existing_outputs)]

  if (length(existing_outputs) > 0L) {
    stop(
      "Collector output already exists. Set OVERWRITE=1 to replace it. ",
      "First existing file: ",
      existing_outputs[1L]
    )
  }
} else {
  unlink(
    c(
      all_adjusted_file,
      all_regions_file,
      all_regions_detailed_file,
      split_qc_file,
      status_qc_file,
      input_manifest_file,
      duplicate_qc_file,
      split_files_existing
    ),
    force = TRUE
  )
}

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir)
}

input_files <- list.files(
  input_dir,
  pattern = "^somnibus_training_splits_job_[0-9]{4}\\.tsv$",
  full.names = TRUE
)

if (length(input_files) == 0L) {
  stop("No SOMNiBUS training-split result files were found in: ", input_dir)
}

input_job_id <- suppressWarnings(
  as.integer(
    sub(
      "^somnibus_training_splits_job_([0-9]{4})\\.tsv$",
      "\\1",
      basename(input_files)
    )
  )
)

input_manifest <- data.table(
  source_file = basename(input_files),
  source_path = input_files,
  prep_job_id_from_filename = input_job_id,
  file_size_bytes = file.info(input_files)$size,
  modified_time = as.character(file.info(input_files)$mtime),
  expected_bundle_count = expected_bundles,
  observed_bundle_count = length(input_files)
)

setorder(input_manifest, prep_job_id_from_filename)

missing_bundle_ids <- setdiff(
  seq_len(expected_bundles),
  input_manifest$prep_job_id_from_filename
)

fwrite(
  input_manifest,
  input_manifest_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Input result files found:", length(input_files), "\n")
cat("Expected bundle files:", expected_bundles, "\n")

if (length(missing_bundle_ids) > 0L) {
  warning(
    "Missing result files for ",
    length(missing_bundle_ids),
    " bundle IDs. First missing IDs: ",
    paste(head(missing_bundle_ids, 20L), collapse = ", ")
  )
}

required_columns <- c(
  "prep_job_id",
  "omnibus_region_id",
  "data_chunk_id",
  "parent_region_id",
  "chr",
  "region_start",
  "region_end",
  "splitID",
  "status",
  "aa_pvalue",
  "aa_EDF",
  "max_steps_reached",
  "n_train_FIDs_observed",
  "n_train_samples",
  "n_train_cases",
  "n_train_controls",
  "n_train_cpgs"
)

read_one_result <- function(file_path) {
  header <- names(
    fread(
      file_path,
      nrows = 0L,
      showProgress = FALSE
    )
  )

  missing_columns <- setdiff(required_columns, header)

  if (length(missing_columns) > 0L) {
    stop(
      "Input file is missing required columns: ",
      basename(file_path),
      " [",
      paste(missing_columns, collapse = ", "),
      "]"
    )
  }

  dt <- fread(
    file_path,
    select = required_columns,
    showProgress = FALSE
  )

  dt[, source_file := basename(file_path)]
  dt
}

result_list <- vector("list", length(input_files))

for (i in seq_along(input_files)) {
  if (i %% 50L == 0L || i == 1L || i == length(input_files)) {
    cat("Reading input file", i, "of", length(input_files), "\n")
  }

  result_list[[i]] <- read_one_result(input_files[i])
}

results_all <- rbindlist(
  result_list,
  use.names = TRUE,
  fill = TRUE
)

rm(result_list)
invisible(gc(FALSE))

if (nrow(results_all) == 0L) {
  stop("The collected SOMNiBUS result table is empty.")
}

results_all[, omnibus_region_id := as.character(omnibus_region_id)]
results_all[, splitID := suppressWarnings(as.integer(splitID))]
results_all[, aa_pvalue := suppressWarnings(as.numeric(aa_pvalue))]
results_all[, aa_EDF := suppressWarnings(as.numeric(aa_EDF))]
results_all[, chr := suppressWarnings(as.integer(chr))]
results_all[, region_start := suppressWarnings(as.integer(region_start))]
results_all[, region_end := suppressWarnings(as.integer(region_end))]

results_all <- results_all[
  !is.na(splitID) &
    splitID >= split_start &
    splitID <= split_end
]

if (nrow(results_all) == 0L) {
  stop("No result rows remain in the requested split range.")
}

status_qc <- results_all[, .N, by = .(splitID, status)]
setorder(status_qc, splitID, status)

fwrite(
  status_qc,
  status_qc_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

results_all[, region_split_key := paste(
  omnibus_region_id,
  splitID,
  sep = "::"
)]

duplicate_rows <- results_all[
  duplicated(region_split_key) |
    duplicated(region_split_key, fromLast = TRUE)
]

if (nrow(duplicate_rows) > 0L) {
  duplicate_qc <- duplicate_rows[
    ,
    .(
      n_rows = .N,
      statuses = paste(sort(unique(status)), collapse = " | "),
      source_files = paste(sort(unique(source_file)), collapse = " | "),
      n_fit_ok = sum(status == "FIT_OK", na.rm = TRUE),
      n_finite_pvalues = sum(is.finite(aa_pvalue), na.rm = TRUE)
    ),
    by = .(omnibus_region_id, splitID)
  ]

  fwrite(
    duplicate_qc,
    duplicate_qc_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  results_all[, duplicate_priority := (
    100L * as.integer(status == "FIT_OK") +
      10L * as.integer(is.finite(aa_pvalue)) +
      as.integer(
        is.na(max_steps_reached) |
          max_steps_reached == FALSE
      )
  )]

  setorder(
    results_all,
    region_split_key,
    -duplicate_priority,
    source_file
  )

  results_all <- results_all[!duplicated(region_split_key)]
  results_all[, duplicate_priority := NULL]
} else {
  fwrite(
    data.table(
      omnibus_region_id = character(),
      splitID = integer(),
      n_rows = integer(),
      statuses = character(),
      source_files = character(),
      n_fit_ok = integer(),
      n_finite_pvalues = integer()
    ),
    duplicate_qc_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

# Match the legacy calibration from 14_eval.R, separately within each split.
legacy_adjust_pvalues <- function(raw_pvalues) {
  p_raw <- suppressWarnings(as.numeric(raw_pvalues))

  valid <- is.finite(p_raw) & p_raw >= 0 & p_raw <= 1

  p_work <- rep(NA_real_, length(p_raw))
  p_transformed <- rep(NA_real_, length(p_raw))
  z_score <- rep(NA_real_, length(p_raw))
  z_adjusted <- rep(NA_real_, length(p_raw))
  p_adjusted <- rep(NA_real_, length(p_raw))

  if (sum(valid) < 2L) {
    return(list(
      p_work = p_work,
      p_transformed = p_transformed,
      z_score = z_score,
      z_adjusted = z_adjusted,
      p_adjusted = p_adjusted,
      mean_z = NA_real_,
      sd_z = NA_real_,
      n_valid = sum(valid),
      n_zero = sum(valid & p_raw == 0),
      min_positive = NA_real_,
      calibration_status = "TOO_FEW_VALID_PVALUES"
    ))
  }

  p_work[valid] <- p_raw[valid]

  positive <- p_work[is.finite(p_work) & p_work > 0]
  min_positive <- if (length(positive) > 0L) {
    min(positive)
  } else {
    .Machine$double.xmin
  }

  # Exact behavior of 14_eval.R: replace zero with the minimum nonzero value.
  p_work[valid & p_work == 0] <- min_positive

  p_work[valid] <- pmin(
    pmax(p_work[valid], .Machine$double.xmin),
    1
  )

  p_transformed[valid] <- exp(
    -sqrt(-log(p_work[valid]))
  )

  p_transformed[valid] <- pmin(
    pmax(p_transformed[valid], .Machine$double.xmin),
    1 - .Machine$double.eps
  )

  z_score[valid] <- qnorm(
    1 - p_transformed[valid] / 2
  )

  finite_z <- is.finite(z_score)

  if (sum(finite_z) < 2L) {
    return(list(
      p_work = p_work,
      p_transformed = p_transformed,
      z_score = z_score,
      z_adjusted = z_adjusted,
      p_adjusted = p_adjusted,
      mean_z = NA_real_,
      sd_z = NA_real_,
      n_valid = sum(valid),
      n_zero = sum(valid & p_raw == 0),
      min_positive = min_positive,
      calibration_status = "TOO_FEW_FINITE_Z_SCORES"
    ))
  }

  mean_z <- mean(z_score[finite_z])
  sd_z <- stats::sd(z_score[finite_z])

  if (!is.finite(sd_z) || sd_z <= sqrt(.Machine$double.eps)) {
    return(list(
      p_work = p_work,
      p_transformed = p_transformed,
      z_score = z_score,
      z_adjusted = z_adjusted,
      p_adjusted = p_adjusted,
      mean_z = mean_z,
      sd_z = sd_z,
      n_valid = sum(valid),
      n_zero = sum(valid & p_raw == 0),
      min_positive = min_positive,
      calibration_status = "ZERO_OR_INVALID_Z_SD"
    ))
  }

  z_adjusted[finite_z] <- (
    z_score[finite_z] - mean_z
  ) / sd_z

  p_adjusted[finite_z] <- 2 * pnorm(
    abs(z_adjusted[finite_z]),
    lower.tail = FALSE
  )

  list(
    p_work = p_work,
    p_transformed = p_transformed,
    z_score = z_score,
    z_adjusted = z_adjusted,
    p_adjusted = p_adjusted,
    mean_z = mean_z,
    sd_z = sd_z,
    n_valid = sum(valid),
    n_zero = sum(valid & p_raw == 0),
    min_positive = min_positive,
    calibration_status = "OK"
  )
}

adjusted_split_list <- vector("list", length(split_ids_expected))
significant_split_list <- vector("list", length(split_ids_expected))
split_qc_list <- vector("list", length(split_ids_expected))

names(adjusted_split_list) <- as.character(split_ids_expected)
names(significant_split_list) <- as.character(split_ids_expected)
names(split_qc_list) <- as.character(split_ids_expected)

for (split_id in split_ids_expected) {
  cat("Processing split", split_id, "\n")

  split_all <- results_all[splitID == split_id]

  split_fit <- split_all[
    status == "FIT_OK" &
      is.finite(aa_pvalue) &
      aa_pvalue >= 0 &
      aa_pvalue <= 1 &
      !is.na(chr) &
      !is.na(region_start) &
      !is.na(region_end)
  ]

  calibration <- legacy_adjust_pvalues(split_fit$aa_pvalue)

  if (nrow(split_fit) > 0L) {
    split_fit[, `:=`(
      aa_pvalue_raw = aa_pvalue,
      aa_pvalue_work = calibration$p_work,
      pvalue_transformed = calibration$p_transformed,
      z_score = calibration$z_score,
      z_adjusted = calibration$z_adjusted,
      pval_adjusted = calibration$p_adjusted,
      pval = calibration$p_adjusted,
      calibration_mean_z = calibration$mean_z,
      calibration_sd_z = calibration$sd_z,
      calibration_status = calibration$calibration_status
    )]
  } else {
    split_fit[, `:=`(
      aa_pvalue_raw = numeric(),
      aa_pvalue_work = numeric(),
      pvalue_transformed = numeric(),
      z_score = numeric(),
      z_adjusted = numeric(),
      pval_adjusted = numeric(),
      pval = numeric(),
      calibration_mean_z = numeric(),
      calibration_sd_z = numeric(),
      calibration_status = character()
    )]
  }

  significant <- split_fit[
    is.finite(pval_adjusted) &
      pval_adjusted < region_p_threshold
  ]

  setorder(significant, pval_adjusted)

  split_region_file <- file.path(
    split_output_dir,
    sprintf("16_somnibus_regions_split_%03d.csv", split_id)
  )

  # Exact column format of the previous 16_somnibus_regions.csv.
  split_region_exact <- significant[, .(
    chr,
    region_start,
    region_end,
    pval_adjusted,
    pval
  )]

  fwrite(
    split_region_exact,
    split_region_file,
    sep = ",",
    quote = FALSE,
    na = "NA"
  )

  adjusted_split_list[[as.character(split_id)]] <- split_fit[, .(
    splitID,
    prep_job_id,
    omnibus_region_id,
    data_chunk_id,
    parent_region_id,
    chr,
    region_start,
    region_end,
    status,
    aa_EDF,
    aa_pvalue_raw,
    aa_pvalue_work,
    pvalue_transformed,
    z_score,
    z_adjusted,
    pval_adjusted,
    pval,
    calibration_mean_z,
    calibration_sd_z,
    calibration_status,
    max_steps_reached,
    n_train_FIDs_observed,
    n_train_samples,
    n_train_cases,
    n_train_controls,
    n_train_cpgs,
    source_file
  )]

  significant_split_list[[as.character(split_id)]] <- significant[, .(
    splitID,
    omnibus_region_id,
    data_chunk_id,
    parent_region_id,
    chr,
    region_start,
    region_end,
    aa_pvalue_raw,
    pval_adjusted,
    pval
  )]

  split_qc_list[[as.character(split_id)]] <- data.table(
    splitID = split_id,
    n_rows_all_statuses = nrow(split_all),
    n_unique_regions_all_statuses = uniqueN(split_all$omnibus_region_id),
    n_fit_ok = split_all[status == "FIT_OK", .N],
    n_valid_fit_pvalues = nrow(split_fit),
    n_invalid_fit_pvalues = split_all[
      status == "FIT_OK" &
        (
          !is.finite(aa_pvalue) |
            aa_pvalue < 0 |
            aa_pvalue > 1
        ),
      .N
    ],
    n_zero_raw_pvalues = calibration$n_zero,
    min_positive_raw_pvalue = calibration$min_positive,
    calibration_mean_z = calibration$mean_z,
    calibration_sd_z = calibration$sd_z,
    calibration_status = calibration$calibration_status,
    pval_adjusted_threshold = region_p_threshold,
    n_regions_selected = nrow(significant)
  )
}

adjusted_all_splits <- rbindlist(
  adjusted_split_list,
  use.names = TRUE,
  fill = TRUE
)

regions_all_splits <- rbindlist(
  significant_split_list,
  use.names = TRUE,
  fill = TRUE
)

split_qc <- rbindlist(
  split_qc_list,
  use.names = TRUE,
  fill = TRUE
)

setorder(
  adjusted_all_splits,
  splitID,
  pval_adjusted,
  omnibus_region_id,
  na.last = TRUE
)

setorder(
  regions_all_splits,
  splitID,
  pval_adjusted,
  omnibus_region_id,
  na.last = TRUE
)

setorder(split_qc, splitID)

fwrite(
  adjusted_all_splits,
  all_adjusted_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

regions_all_splits_minimal <- regions_all_splits[, .(
  splitID,
  chr,
  region_start,
  region_end,
  pval_adjusted,
  pval
)]

fwrite(
  regions_all_splits_minimal,
  all_regions_file,
  sep = ",",
  quote = FALSE,
  na = "NA"
)

fwrite(
  regions_all_splits,
  all_regions_detailed_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  split_qc,
  split_qc_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\nCollection finished.\n")
cat("All adjusted results:", all_adjusted_file, "\n")
cat("All significant split-specific regions:", all_regions_file, "\n")
cat("Detailed significant-region table:", all_regions_detailed_file, "\n")
cat("Per-split region directory:", split_output_dir, "\n")
cat("Split QC:", split_qc_file, "\n")
cat("Status QC:", status_qc_file, "\n")
cat("Input manifest:", input_manifest_file, "\n")
cat("Duplicate QC:", duplicate_qc_file, "\n")
cat("Collected rows across all statuses:", nrow(results_all), "\n")
cat("Adjusted FIT_OK rows:", nrow(adjusted_all_splits), "\n")
cat("Selected region-split rows:", nrow(regions_all_splits), "\n")
cat("Missing bundle result files:", length(missing_bundle_ids), "\n")

print(
  split_qc[, .(
    splits = .N,
    median_fit_ok = median(n_fit_ok),
    median_selected = median(n_regions_selected),
    min_selected = min(n_regions_selected),
    max_selected = max(n_regions_selected),
    calibration_failures = sum(calibration_status != "OK")
  )]
)

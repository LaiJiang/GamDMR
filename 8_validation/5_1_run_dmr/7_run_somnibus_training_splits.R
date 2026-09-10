#!/usr/bin/env Rscript

# Fit SOMNiBUS separately within each family-level training split.
# One SLURM array task reads one prepared SOMNiBUS bundle.
# Results are appended to a task-specific TSV and can be resumed safely.

# --------------------------- R library path ----------------------------------
user_lib <- path.expand("~/scratch/R/library")
.libPaths(unique(c(user_lib, .libPaths())))

suppressPackageStartupMessages({
  library(data.table)
  library(SOMNiBUS)
})

cat("R library paths:\n")
print(.libPaths())
cat("SOMNiBUS version:", as.character(packageVersion("SOMNiBUS")), "\n")

# --------------------------- configuration -----------------------------------
PATH_wk <- path.expand("~/scratch/UQAC/meth/")

bundle_dir <- file.path(
  PATH_wk,
  "results/15_revision/5_cv/2_somnibus/prep/bundles"
)

train_split_file <- file.path(
  PATH_wk,
  "scr/11_mgcv/dat/100_family_training_splits.csv"
)

output_dir <- file.path(
  PATH_wk,
  "results/15_revision/5_cv/2_somnibus/training_splits"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Production bundle partitioning.
N_bundles <- as.integer(Sys.getenv("N_BUNDLES", "900"))

# SOMNiBUS settings matching the previous analysis.
min_cpgs <- as.integer(Sys.getenv("SOMNIBUS_MIN_CPGS", "51"))
p0_value <- as.numeric(Sys.getenv("SOMNIBUS_P0", "0.003"))
p1_value <- as.numeric(Sys.getenv("SOMNIBUS_P1", "0.9"))
epsilon_value <- as.numeric(Sys.getenv("SOMNIBUS_EPSILON", "1e-6"))
epsilon_lambda_value <- as.numeric(
  Sys.getenv("SOMNIBUS_EPSILON_LAMBDA", "1e-3")
)
max_step <- as.integer(Sys.getenv("SOMNIBUS_MAX_STEP", "200"))

# Resume and smoke-test controls.
overwrite <- identical(Sys.getenv("OVERWRITE", "0"), "1")
retry_errors <- identical(Sys.getenv("RETRY_ERRORS", "0"), "1")
flush_every <- as.integer(Sys.getenv("FLUSH_EVERY", "10"))
max_regions <- as.integer(Sys.getenv("MAX_REGIONS", "0"))
split_start <- as.integer(Sys.getenv("SPLIT_START", "1"))
split_end <- as.integer(Sys.getenv("SPLIT_END", "2147483647"))
max_runtime_minutes <- as.numeric(Sys.getenv("MAX_RUNTIME_MINUTES", "0"))

if (is.na(N_bundles) || N_bundles < 1L) stop("N_BUNDLES must be >= 1.")
if (is.na(min_cpgs) || min_cpgs < 1L) stop("SOMNIBUS_MIN_CPGS must be >= 1.")
if (!is.finite(p0_value) || p0_value < 0 || p0_value >= 1) {
  stop("SOMNIBUS_P0 must be in [0, 1).")
}
if (!is.finite(p1_value) || p1_value <= 0 || p1_value > 1) {
  stop("SOMNIBUS_P1 must be in (0, 1].")
}
if (p0_value >= p1_value) stop("SOMNIBUS_P0 must be smaller than SOMNIBUS_P1.")
if (is.na(max_step) || max_step < 1L) stop("SOMNIBUS_MAX_STEP must be >= 1.")
if (is.na(flush_every) || flush_every < 1L) stop("FLUSH_EVERY must be >= 1.")
if (is.na(max_regions) || max_regions < 0L) stop("MAX_REGIONS must be >= 0.")
if (is.na(split_start) || split_start < 1L) stop("SPLIT_START must be >= 1.")
if (is.na(split_end) || split_end < split_start) {
  stop("SPLIT_END must be >= SPLIT_START.")
}
if (!is.finite(max_runtime_minutes) || max_runtime_minutes < 0) {
  stop("MAX_RUNTIME_MINUTES must be >= 0.")
}

# --------------------------- task ID -----------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Pass SLURM_ARRAY_TASK_ID as the first argument.")

job_id <- suppressWarnings(as.integer(args[1]))
if (is.na(job_id) || job_id < 1L || job_id > N_bundles) {
  stop(
    "Invalid SLURM_ARRAY_TASK_ID. Expected an integer from 1 to ",
    N_bundles, "."
  )
}

cat("Bundle task ID:", job_id, "\n")
job_start_time <- Sys.time()

bundle_file <- file.path(
  bundle_dir,
  sprintf("somnibus_bundle_job_%04d.rds", job_id)
)

result_file <- file.path(
  output_dir,
  sprintf("somnibus_training_splits_job_%04d.tsv", job_id)
)

summary_file <- file.path(
  output_dir,
  sprintf("somnibus_training_splits_summary_job_%04d.tsv", job_id)
)

session_file <- file.path(
  output_dir,
  sprintf("somnibus_training_splits_session_job_%04d.txt", job_id)
)

if (!file.exists(bundle_file)) {
  stop("Bundle file does not exist: ", bundle_file)
}
if (!file.exists(train_split_file)) {
  stop("Training-split file does not exist: ", train_split_file)
}

if (overwrite) {
  unlink(c(result_file, summary_file, session_file), force = TRUE)
}

# --------------------------- helper functions --------------------------------
clean_text <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_character_)
  x <- paste(unique(as.character(x)), collapse = " | ")
  x <- gsub("[\r\n\t]+", " ", x)
  trimws(x)
}

as_scalar_numeric <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  suppressWarnings(as.numeric(x[1L]))
}

matrix_value <- function(x, row_name, column_name) {
  if (
    is.null(x) || is.null(rownames(x)) || is.null(colnames(x)) ||
      !row_name %in% rownames(x) || !column_name %in% colnames(x)
  ) {
    return(NA_real_)
  }
  as_scalar_numeric(x[row_name, column_name])
}

make_key <- function(region_id, split_id) {
  paste(region_id, split_id, sep = "::")
}

make_result_row <- function(
  prep_job_id = job_id,
  bundle_index = NA_integer_,
  omnibus_region_id = NA_character_,
  data_chunk_id = NA_integer_,
  parent_region_id = NA_integer_,
  chr = NA_integer_,
  region_start = NA_integer_,
  region_end = NA_integer_,
  splitID = NA_integer_,
  status = NA_character_,
  message = NA_character_,
  warnings = NA_character_,
  n_train_FIDs_requested = NA_integer_,
  n_train_FIDs_observed = NA_integer_,
  n_train_samples = NA_integer_,
  n_train_cases = NA_integer_,
  n_train_controls = NA_integer_,
  n_train_cpgs = NA_integer_,
  n_train_rows = NA_integer_,
  n_k = NA_integer_,
  aa_EDF = NA_real_,
  aa_statistic = NA_real_,
  aa_statistic_type = NA_character_,
  aa_pvalue = NA_real_,
  intercept_pvalue = NA_real_,
  phi_fletcher = NA_real_,
  phi_reml = NA_real_,
  phi_gam = NA_real_,
  lambda_intercept = NA_real_,
  lambda_AA = NA_real_,
  n_iterations = NA_integer_,
  max_steps_reached = NA,
  mean_control_prob = NA_real_,
  mean_case_prob = NA_real_,
  max_abs_prob_diff = NA_real_,
  mean_abs_prob_diff = NA_real_,
  mean_signed_prob_diff = NA_real_,
  max_abs_logit_effect = NA_real_,
  mean_abs_logit_effect = NA_real_,
  fit_seconds = NA_real_
) {
  data.table(
    prep_job_id = as.integer(prep_job_id),
    bundle_index = as.integer(bundle_index),
    omnibus_region_id = as.character(omnibus_region_id),
    data_chunk_id = as.integer(data_chunk_id),
    parent_region_id = as.integer(parent_region_id),
    chr = as.integer(chr),
    region_start = as.integer(region_start),
    region_end = as.integer(region_end),
    splitID = as.integer(splitID),
    status = as.character(status),
    message = as.character(message),
    warnings = as.character(warnings),
    n_train_FIDs_requested = as.integer(n_train_FIDs_requested),
    n_train_FIDs_observed = as.integer(n_train_FIDs_observed),
    n_train_samples = as.integer(n_train_samples),
    n_train_cases = as.integer(n_train_cases),
    n_train_controls = as.integer(n_train_controls),
    n_train_cpgs = as.integer(n_train_cpgs),
    n_train_rows = as.integer(n_train_rows),
    n_k = as.integer(n_k),
    aa_EDF = as.numeric(aa_EDF),
    aa_statistic = as.numeric(aa_statistic),
    aa_statistic_type = as.character(aa_statistic_type),
    aa_pvalue = as.numeric(aa_pvalue),
    intercept_pvalue = as.numeric(intercept_pvalue),
    phi_fletcher = as.numeric(phi_fletcher),
    phi_reml = as.numeric(phi_reml),
    phi_gam = as.numeric(phi_gam),
    lambda_intercept = as.numeric(lambda_intercept),
    lambda_AA = as.numeric(lambda_AA),
    n_iterations = as.integer(n_iterations),
    max_steps_reached = as.logical(max_steps_reached),
    mean_control_prob = as.numeric(mean_control_prob),
    mean_case_prob = as.numeric(mean_case_prob),
    max_abs_prob_diff = as.numeric(max_abs_prob_diff),
    mean_abs_prob_diff = as.numeric(mean_abs_prob_diff),
    mean_signed_prob_diff = as.numeric(mean_signed_prob_diff),
    max_abs_logit_effect = as.numeric(max_abs_logit_effect),
    mean_abs_logit_effect = as.numeric(mean_abs_logit_effect),
    fit_seconds = as.numeric(fit_seconds)
  )
}

result_template <- make_result_row()[0]

# Create a header-only file so empty tasks still have a valid output.
if (!file.exists(result_file)) {
  fwrite(result_template, result_file, sep = "\t", quote = FALSE, na = "NA")
}

completed_keys <- new.env(hash = TRUE, parent = emptyenv())

if (file.info(result_file)$size > 0L) {
  previous_results <- tryCatch(
    fread(result_file, showProgress = FALSE),
    error = function(e) e
  )

  if (inherits(previous_results, "error")) {
    stop("Cannot read existing result file: ", conditionMessage(previous_results))
  }

  required_previous <- c("omnibus_region_id", "splitID", "status")
  if (!all(required_previous %in% names(previous_results))) {
    stop("Existing result file has an incompatible format: ", result_file)
  }

  if (nrow(previous_results) > 0L) {
    if (retry_errors) {
      retry_statuses <- c("FIT_ERROR", "FIT_OUTPUT_ERROR")
      if (any(previous_results$status %chin% retry_statuses)) {
        previous_results <- previous_results[
          !status %chin% retry_statuses
        ]
        fwrite(
          previous_results,
          result_file,
          sep = "\t",
          quote = FALSE,
          na = "NA"
        )
      }
    }

    previous_key_values <- make_key(
      previous_results$omnibus_region_id,
      previous_results$splitID
    )

    for (k in previous_key_values) completed_keys[[k]] <- TRUE

    cat(
      "Existing completed region-split rows:",
      length(previous_key_values),
      "\n"
    )
  }

  rm(previous_results)
}

pending_rows <- list()

flush_pending <- function(force = FALSE) {
  if (length(pending_rows) == 0L) return(invisible(NULL))
  if (!force && length(pending_rows) < flush_every) return(invisible(NULL))

  out_dt <- rbindlist(pending_rows, use.names = TRUE, fill = TRUE)
  fwrite(
    out_dt,
    result_file,
    append = TRUE,
    col.names = FALSE,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  pending_rows <<- list()
  invisible(NULL)
}

add_result <- function(row) {
  pending_rows[[length(pending_rows) + 1L]] <<- row
  completed_keys[[make_key(row$omnibus_region_id, row$splitID)]] <- TRUE
  flush_pending(force = FALSE)
  invisible(NULL)
}

runtime_limit_reached <- function() {
  if (max_runtime_minutes <= 0) return(FALSE)
  elapsed_minutes <- as.numeric(
    difftime(Sys.time(), job_start_time, units = "mins")
  )
  elapsed_minutes >= max_runtime_minutes
}

# --------------------------- load inputs -------------------------------------
bundle <- readRDS(bundle_file)

required_bundle <- c("metadata", "sample_map", "manifest_ready", "regions")
missing_bundle <- setdiff(required_bundle, names(bundle))
if (length(missing_bundle) > 0L) {
  stop("Bundle is missing: ", paste(missing_bundle, collapse = ", "))
}

if (!is.list(bundle$regions)) stop("bundle$regions must be a list.")
if (length(bundle$regions) > 0L) {
  if (is.null(names(bundle$regions)) || any(!nzchar(names(bundle$regions)))) {
    stop("Every bundled region must have an omnibus_region_id name.")
  }
  if (anyDuplicated(names(bundle$regions))) {
    stop("Duplicated omnibus_region_id values in bundle$regions.")
  }
}

sample_map <- as.data.table(copy(bundle$sample_map))
required_sample_map <- c("ID", "FID", "AA_only")
missing_sample_map <- setdiff(required_sample_map, names(sample_map))
if (length(missing_sample_map) > 0L) {
  stop("bundle$sample_map is missing: ", paste(missing_sample_map, collapse = ", "))
}

sample_map <- sample_map[, ..required_sample_map]
sample_map[, ID := as.character(ID)]
sample_map[, FID := as.character(FID)]
sample_map[, AA_only := suppressWarnings(as.integer(as.character(AA_only)))]

if (anyNA(sample_map$ID) || any(!nzchar(sample_map$ID))) {
  stop("Missing or empty IDs in bundle$sample_map.")
}
if (anyNA(sample_map$FID) || any(!nzchar(sample_map$FID))) {
  stop("Missing or empty FIDs in bundle$sample_map.")
}
if (anyNA(sample_map$AA_only) || any(!sample_map$AA_only %in% c(0L, 1L))) {
  stop("AA_only in bundle$sample_map must contain only 0 and 1.")
}
if (anyDuplicated(sample_map$ID)) stop("Duplicated IDs in bundle$sample_map.")

manifest_dt <- as.data.table(copy(bundle$manifest_ready))
if (nrow(manifest_dt) > 0L) {
  required_manifest <- c(
    "bundle_index", "omnibus_region_id", "data_chunk_id",
    "parent_region_id", "chr", "region_start", "region_end"
  )
  missing_manifest <- setdiff(required_manifest, names(manifest_dt))
  if (length(missing_manifest) > 0L) {
    stop(
      "bundle$manifest_ready is missing: ",
      paste(missing_manifest, collapse = ", ")
    )
  }
  manifest_dt[, omnibus_region_id := as.character(omnibus_region_id)]
  if (anyDuplicated(manifest_dt$omnibus_region_id)) {
    stop("Duplicated omnibus_region_id values in bundle$manifest_ready.")
  }
  setkey(manifest_dt, omnibus_region_id)
}

train_FID_df <- fread(train_split_file, showProgress = FALSE)
if (!"splitID" %in% names(train_FID_df)) {
  stop("Training-split file must contain splitID.")
}

train_FID_cols <- grep(
  "^train_FID_[0-9]+$",
  names(train_FID_df),
  value = TRUE
)
if (length(train_FID_cols) == 0L) {
  stop("No train_FID_* columns were found in the training-split file.")
}

train_FID_col_number <- as.integer(sub("^train_FID_", "", train_FID_cols))
train_FID_cols <- train_FID_cols[order(train_FID_col_number)]
train_FID_df[, splitID := suppressWarnings(as.integer(splitID))]

if (anyNA(train_FID_df$splitID)) stop("splitID contains missing/non-integer values.")
if (anyDuplicated(train_FID_df$splitID)) stop("splitID values must be unique.")

train_FID_df <- train_FID_df[
  splitID >= split_start & splitID <= split_end
]
setorder(train_FID_df, splitID)

if (nrow(train_FID_df) == 0L) {
  stop("No splitID values remain after applying SPLIT_START and SPLIT_END.")
}

# Precompute requested FIDs and sample IDs once for each split.
split_info <- vector("list", nrow(train_FID_df))
names(split_info) <- as.character(train_FID_df$splitID)

for (i in seq_len(nrow(train_FID_df))) {
  split_id <- train_FID_df$splitID[i]
  raw_fids <- unlist(
    train_FID_df[i, ..train_FID_cols],
    use.names = FALSE
  )
  train_fids <- unique(as.character(raw_fids))
  train_fids <- train_fids[
    !is.na(train_fids) & nzchar(train_fids) & train_fids != "NA"
  ]

  train_ids <- sample_map[FID %chin% train_fids, unique(ID)]

  split_info[[as.character(split_id)]] <- list(
    splitID = split_id,
    train_FIDs = train_fids,
    train_IDs = train_ids
  )
}

region_names <- names(bundle$regions)
if (max_regions > 0L && length(region_names) > max_regions) {
  region_names <- region_names[seq_len(max_regions)]
}

cat("Regions selected for this run:", length(region_names), "\n")
cat("Splits selected for this run:", nrow(train_FID_df), "\n")
cat("Result file:", result_file, "\n")

stop_requested <- FALSE

# --------------------------- fit models --------------------------------------
for (region_name in region_names) {
  if (runtime_limit_reached()) {
    cat("Runtime limit reached before region", region_name, "\n")
    stop_requested <- TRUE
    break
  }

  region_all <- as.data.table(copy(bundle$regions[[region_name]]))
  required_region_cols <- c(
    "Meth_Counts", "Total_Counts", "Position", "ID", "AA_only"
  )
  missing_region_cols <- setdiff(required_region_cols, names(region_all))

  meta <- if (nrow(manifest_dt) > 0L) {
    manifest_dt[J(region_name), nomatch = 0L]
  } else {
    data.table()
  }

  if (nrow(meta) == 0L) {
    meta <- data.table(
      bundle_index = match(region_name, names(bundle$regions)),
      omnibus_region_id = region_name,
      data_chunk_id = NA_integer_,
      parent_region_id = NA_integer_,
      chr = NA_integer_,
      region_start = NA_integer_,
      region_end = NA_integer_
    )
  }

  base_args <- list(
    bundle_index = meta$bundle_index[1L],
    omnibus_region_id = region_name,
    data_chunk_id = meta$data_chunk_id[1L],
    parent_region_id = meta$parent_region_id[1L],
    chr = meta$chr[1L],
    region_start = meta$region_start[1L],
    region_end = meta$region_end[1L]
  )

  if (length(missing_region_cols) > 0L) {
    for (split_id in train_FID_df$splitID) {
      key <- make_key(region_name, split_id)
      if (!is.null(completed_keys[[key]])) next

      add_result(do.call(
        make_result_row,
        c(
          base_args,
          list(
            splitID = split_id,
            status = "INVALID_REGION_DATA",
            message = paste(
              "Missing regional columns:",
              paste(missing_region_cols, collapse = ", ")
            )
          )
        )
      ))
    }
    next
  }

  region_all <- region_all[, ..required_region_cols]
  region_all[, ID := as.character(ID)]
  region_all[, Position := suppressWarnings(as.integer(Position))]
  region_all[, Meth_Counts := suppressWarnings(as.numeric(Meth_Counts))]
  region_all[, Total_Counts := suppressWarnings(as.numeric(Total_Counts))]
  region_all[, AA_only := suppressWarnings(as.integer(as.character(AA_only)))]

  region_all <- region_all[
    !is.na(ID) & nzchar(ID) &
      is.finite(Position) &
      is.finite(Meth_Counts) &
      is.finite(Total_Counts) &
      Total_Counts > 0 &
      Meth_Counts >= 0 &
      Meth_Counts <= Total_Counts &
      !is.na(AA_only) &
      AA_only %in% c(0L, 1L)
  ]

  invalid_region_reason <- NA_character_

  if (nrow(region_all) == 0L) {
    invalid_region_reason <- "No valid observations remain after regional filtering."
  } else if (anyDuplicated(region_all[, .(ID, Position)])) {
    invalid_region_reason <- "Duplicated ID-Position rows are present."
  } else {
    region_sample_aa <- unique(region_all[, .(ID, AA_only)])
    if (anyDuplicated(region_sample_aa$ID)) {
      invalid_region_reason <- "At least one ID has inconsistent AA_only values."
    } else {
      aa_check <- merge(
        region_sample_aa,
        sample_map[, .(ID, AA_only_map = AA_only)],
        by = "ID",
        all.x = TRUE,
        sort = FALSE
      )
      if (anyNA(aa_check$AA_only_map)) {
        invalid_region_reason <- "At least one regional ID is absent from sample_map."
      } else if (any(aa_check$AA_only != aa_check$AA_only_map)) {
        invalid_region_reason <- "AA_only differs between regional data and sample_map."
      }
    }
  }

  if (!is.na(invalid_region_reason)) {
    for (split_id in train_FID_df$splitID) {
      key <- make_key(region_name, split_id)
      if (!is.null(completed_keys[[key]])) next

      add_result(do.call(
        make_result_row,
        c(
          base_args,
          list(
            splitID = split_id,
            status = "INVALID_REGION_DATA",
            message = invalid_region_reason
          )
        )
      ))
    }
    next
  }

  setorder(region_all, Position, ID)

  cat("\nRegion:", region_name, "\n")

  for (split_id in train_FID_df$splitID) {
    if (runtime_limit_reached()) {
      cat("Runtime limit reached within region", region_name, "\n")
      stop_requested <- TRUE
      break
    }

    key <- make_key(region_name, split_id)
    if (!is.null(completed_keys[[key]])) next

    split_obj <- split_info[[as.character(split_id)]]
    train_fids <- split_obj$train_FIDs
    train_ids <- split_obj$train_IDs

    if (length(train_fids) == 0L) {
      add_result(do.call(
        make_result_row,
        c(
          base_args,
          list(
            splitID = split_id,
            status = "NO_TRAIN_FIDS",
            message = "No training FIDs were supplied for this split.",
            n_train_FIDs_requested = 0L
          )
        )
      ))
      next
    }

    if (length(train_ids) == 0L) {
      add_result(do.call(
        make_result_row,
        c(
          base_args,
          list(
            splitID = split_id,
            status = "NO_TRAIN_IDS",
            message = "No sample_map IDs matched the requested training FIDs.",
            n_train_FIDs_requested = length(train_fids)
          )
        )
      ))
      next
    }

    region_train <- region_all[ID %chin% train_ids]

    observed_ids <- unique(region_train$ID)
    observed_map <- sample_map[ID %chin% observed_ids]
    n_train_fids_observed <- uniqueN(observed_map$FID)
    n_train_samples <- uniqueN(region_train$ID)
    n_train_cpgs <- uniqueN(region_train$Position)
    n_train_rows <- nrow(region_train)

    sample_outcomes <- unique(region_train[, .(ID, AA_only)])
    n_train_cases <- sample_outcomes[AA_only == 1L, uniqueN(ID)]
    n_train_controls <- sample_outcomes[AA_only == 0L, uniqueN(ID)]

    common_counts <- list(
      n_train_FIDs_requested = length(train_fids),
      n_train_FIDs_observed = n_train_fids_observed,
      n_train_samples = n_train_samples,
      n_train_cases = n_train_cases,
      n_train_controls = n_train_controls,
      n_train_cpgs = n_train_cpgs,
      n_train_rows = n_train_rows
    )

    if (n_train_rows == 0L) {
      add_result(do.call(
        make_result_row,
        c(
          base_args,
          list(
            splitID = split_id,
            status = "NO_TRAIN_OBSERVATIONS",
            message = "No regional observations matched the training IDs."
          ),
          common_counts
        )
      ))
      next
    }

    if (n_train_cpgs < min_cpgs) {
      add_result(do.call(
        make_result_row,
        c(
          base_args,
          list(
            splitID = split_id,
            status = "BELOW_MIN_CPGS",
            message = paste0(
              "Training CpGs = ", n_train_cpgs,
              "; required = ", min_cpgs, "."
            )
          ),
          common_counts
        )
      ))
      next
    }

    if (n_train_samples < 2L) {
      add_result(do.call(
        make_result_row,
        c(
          base_args,
          list(
            splitID = split_id,
            status = "TOO_FEW_TRAIN_SAMPLES",
            message = "Fewer than two training samples remain."
          ),
          common_counts
        )
      ))
      next
    }

    if (n_train_fids_observed < 2L) {
      add_result(do.call(
        make_result_row,
        c(
          base_args,
          list(
            splitID = split_id,
            status = "TOO_FEW_TRAIN_FIDS",
            message = "Fewer than two training FIDs remain."
          ),
          common_counts
        )
      ))
      next
    }

    if (n_train_cases == 0L || n_train_controls == 0L) {
      add_result(do.call(
        make_result_row,
        c(
          base_args,
          list(
            splitID = split_id,
            status = "ONE_OUTCOME_CLASS",
            message = "AA_only has fewer than two classes in training samples."
          ),
          common_counts
        )
      ))
      next
    }

    n_k_dim <- max(3L, as.integer(n_train_cpgs / 20L))

    # SOMNiBUS requires these first four columns in this exact order.
    region_train_model <- as.data.frame(
      region_train[
        ,
        .(Meth_Counts, Total_Counts, Position, ID, AA_only)
      ]
    )

    fit_warnings <- character()
    fit_start <- proc.time()[["elapsed"]]

    fit <- tryCatch(
      withCallingHandlers(
        SOMNiBUS::binomRegMethModel(
          data = region_train_model,
          n.k = rep(n_k_dim, 2L),
          p0 = p0_value,
          p1 = p1_value,
          Quasi = FALSE,
          epsilon = epsilon_value,
          epsilon.lambda = epsilon_lambda_value,
          maxStep = max_step,
          binom.link = "logit",
          method = "REML",
          covs = "AA_only",
          RanEff = FALSE,
          reml.scale = FALSE,
          scale = -2,
          verbose = FALSE
        ),
        warning = function(w) {
          fit_warnings <<- c(fit_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )

    fit_seconds <- proc.time()[["elapsed"]] - fit_start
    warnings_text <- clean_text(fit_warnings)

    if (inherits(fit, "error")) {
      add_result(do.call(
        make_result_row,
        c(
          base_args,
          list(
            splitID = split_id,
            status = "FIT_ERROR",
            message = clean_text(conditionMessage(fit)),
            warnings = warnings_text,
            n_k = n_k_dim,
            fit_seconds = fit_seconds
          ),
          common_counts
        )
      ))

      cat(
        "  split", split_id, "FIT_ERROR:",
        clean_text(conditionMessage(fit)), "\n"
      )
      rm(region_train, region_train_model, fit)
      invisible(gc(verbose = FALSE))
      next
    }

    required_fit_objects <- c(
      "reg.out", "Beta.out", "uni.pos", "phi_fletcher",
      "phi_reml", "phi_gam", "lambda", "ite.points"
    )
    missing_fit_objects <- setdiff(required_fit_objects, names(fit))

    if (length(missing_fit_objects) > 0L) {
      add_result(do.call(
        make_result_row,
        c(
          base_args,
          list(
            splitID = split_id,
            status = "FIT_OUTPUT_ERROR",
            message = paste(
              "SOMNiBUS output is missing:",
              paste(missing_fit_objects, collapse = ", ")
            ),
            warnings = warnings_text,
            n_k = n_k_dim,
            fit_seconds = fit_seconds
          ),
          common_counts
        )
      ))
      rm(region_train, region_train_model, fit)
      invisible(gc(verbose = FALSE))
      next
    }

    if (
      is.null(rownames(fit$reg.out)) ||
        !"AA_only" %in% rownames(fit$reg.out) ||
        is.null(colnames(fit$Beta.out)) ||
        !all(c("Intercept", "AA_only") %in% colnames(fit$Beta.out))
    ) {
      add_result(do.call(
        make_result_row,
        c(
          base_args,
          list(
            splitID = split_id,
            status = "FIT_OUTPUT_ERROR",
            message = paste(
              "AA_only was not found in reg.out or Beta.out.",
              "reg.out rows:", paste(rownames(fit$reg.out), collapse = ","),
              "Beta.out columns:", paste(colnames(fit$Beta.out), collapse = ",")
            ),
            warnings = warnings_text,
            n_k = n_k_dim,
            fit_seconds = fit_seconds
          ),
          common_counts
        )
      ))
      rm(region_train, region_train_model, fit)
      invisible(gc(verbose = FALSE))
      next
    }

    statistic_type <- intersect(
      c("Chi.sq", "F"),
      colnames(fit$reg.out)
    )
    statistic_type <- if (length(statistic_type) > 0L) {
      statistic_type[1L]
    } else {
      NA_character_
    }

    aa_edf <- matrix_value(fit$reg.out, "AA_only", "EDF")
    aa_statistic <- if (!is.na(statistic_type)) {
      matrix_value(fit$reg.out, "AA_only", statistic_type)
    } else {
      NA_real_
    }
    aa_pvalue <- matrix_value(fit$reg.out, "AA_only", "p-value")
    intercept_pvalue <- matrix_value(fit$reg.out, "Intercept", "p-value")

    beta_matrix <- as.matrix(fit$Beta.out)
    intercept_effect <- as.numeric(beta_matrix[, "Intercept"])
    aa_logit_effect <- as.numeric(beta_matrix[, "AA_only"])

    control_prob <- plogis(intercept_effect)
    case_prob <- plogis(intercept_effect + aa_logit_effect)
    prob_diff <- case_prob - control_prob

    ite_points <- fit$ite.points
    n_iterations <- if (is.null(dim(ite_points))) {
      1L
    } else {
      as.integer(nrow(ite_points))
    }

    max_steps_reached <- isTRUE(n_iterations >= max_step)

    result_row <- do.call(
      make_result_row,
      c(
        base_args,
        list(
          splitID = split_id,
          status = "FIT_OK",
          message = NA_character_,
          warnings = warnings_text,
          n_k = n_k_dim,
          aa_EDF = aa_edf,
          aa_statistic = aa_statistic,
          aa_statistic_type = statistic_type,
          aa_pvalue = aa_pvalue,
          intercept_pvalue = intercept_pvalue,
          phi_fletcher = as_scalar_numeric(fit$phi_fletcher),
          phi_reml = as_scalar_numeric(fit$phi_reml),
          phi_gam = as_scalar_numeric(fit$phi_gam),
          lambda_intercept = as_scalar_numeric(fit$lambda[1L]),
          lambda_AA = as_scalar_numeric(fit$lambda[2L]),
          n_iterations = n_iterations,
          max_steps_reached = max_steps_reached,
          mean_control_prob = mean(control_prob, na.rm = TRUE),
          mean_case_prob = mean(case_prob, na.rm = TRUE),
          max_abs_prob_diff = max(abs(prob_diff), na.rm = TRUE),
          mean_abs_prob_diff = mean(abs(prob_diff), na.rm = TRUE),
          mean_signed_prob_diff = mean(prob_diff, na.rm = TRUE),
          max_abs_logit_effect = max(abs(aa_logit_effect), na.rm = TRUE),
          mean_abs_logit_effect = mean(abs(aa_logit_effect), na.rm = TRUE),
          fit_seconds = fit_seconds
        ),
        common_counts
      )
    )

    add_result(result_row)

    cat(
      "  completed split", split_id,
      "p =", format(aa_pvalue, digits = 4),
      "CpGs =", n_train_cpgs,
      "samples =", n_train_samples,
      "FIDs =", n_train_fids_observed,
      "seconds =", format(fit_seconds, digits = 4),
      "\n"
    )

    rm(
      region_train, region_train_model, fit, beta_matrix,
      intercept_effect, aa_logit_effect, control_prob, case_prob,
      prob_diff, result_row
    )
    invisible(gc(verbose = FALSE))
  }

  rm(region_all)
  invisible(gc(verbose = FALSE))

  if (stop_requested) break
}

flush_pending(force = TRUE)

# --------------------------- task summary ------------------------------------
all_results <- fread(result_file, showProgress = FALSE)

if (nrow(all_results) > 0L) {
  status_summary <- all_results[, .N, by = status][order(status)]
} else {
  status_summary <- data.table(status = character(), N = integer())
}

job_summary <- data.table(
  prep_job_id = job_id,
  bundle_file = bundle_file,
  result_file = result_file,
  n_bundle_regions = length(bundle$regions),
  n_regions_selected_this_run = length(region_names),
  n_splits_selected_this_run = nrow(train_FID_df),
  n_result_rows_total = nrow(all_results),
  n_fit_ok_total = all_results[status == "FIT_OK", .N],
  runtime_limit_reached = stop_requested,
  elapsed_minutes = as.numeric(
    difftime(Sys.time(), job_start_time, units = "mins")
  ),
  completed_at = format(Sys.time(), tz = "America/Toronto", usetz = TRUE)
)

fwrite(job_summary, summary_file, sep = "\t", quote = FALSE, na = "NA")

session_lines <- capture.output({
  cat("Job summary:\n")
  print(job_summary)
  cat("\nStatus summary:\n")
  print(status_summary)
  cat("\nSession information:\n")
  print(sessionInfo())
})
writeLines(session_lines, session_file)

cat("\nFinished bundle task", job_id, "\n")
cat("Results:", result_file, "\n")
cat("Summary:", summary_file, "\n")
cat("Session info:", session_file, "\n")
cat("Total result rows:", nrow(all_results), "\n")
cat("FIT_OK rows:", all_results[status == "FIT_OK", .N], "\n")
if (nrow(status_summary) > 0L) print(status_summary)

if (stop_requested) {
  cat(
    "The configured runtime limit was reached. Resubmit the same task to resume.\n"
  )
}

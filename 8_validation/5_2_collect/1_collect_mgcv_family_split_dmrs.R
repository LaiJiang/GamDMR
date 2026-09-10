#!/usr/bin/env Rscript

# Collect MGCV family-training results from the forward and reverse batches.
# For each splitID, calculate BH FDR across regions and generate a strict DMR
# file using the same criteria as 24_analyze_BMI.R.

suppressPackageStartupMessages({
  library(data.table)
})

setDTthreads(1L)

# ------------------------------- configuration -------------------------------

PATH_wk <- path.expand(
  Sys.getenv("PATH_WK", "~/scratch/UQAC/meth/")
)

PATH_results_B1 <- Sys.getenv(
  "MGCV_BATCH1_DIR",
  file.path(PATH_wk, "results/15_revision/5_cv/1_mgcv")
)

PATH_results_B2 <- Sys.getenv(
  "MGCV_BATCH2_DIR",
  file.path(PATH_wk, "results/15_revision/5_cv/1_mgcv/2_batch")
)

region_file_path <- Sys.getenv(
  "MGCV_REGION_FILE",
  file.path(PATH_wk, "scr/11_mgcv/dat/region_GAM_M12.csv")
)

PATH_output <- Sys.getenv(
  "MGCV_COLLECT_OUTPUT",
  file.path(
    PATH_wk,
    "results/15_revision/5_cv/1_mgcv/collected_by_split"
  )
)

strict_output_dir <- file.path(PATH_output, "strict_by_split")
qc_output_dir <- file.path(PATH_output, "qc")

dir.create(strict_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_output_dir, recursive = TRUE, showWarnings = FALSE)

# These reproduce the strict criteria in 24_analyze_BMI.R.
raw_p_threshold <- 1e-5
fdr_threshold <- 0.01
edf_threshold <- 1
r2_threshold <- 0.30
min_cpgs_threshold <- 10L
mean_diff_threshold <- 0.01

# ---------------------------- result column names -----------------------------

# The model-fitting scripts write headerless files in this exact order.
result_colnames <- c(
  "splitID",
  "(Intercept)",
  "AA_only",
  "AgeCalc",
  "Sex",
  "Non.smoker",
  "EOSINOpc",
  "LYMPHOpc",
  "MONOpc",
  "NEUTROpc",
  "sv1",
  "sv2",
  "sv3",
  "sv4",
  "sv5",
  "BMI",
  "s(start)",
  "s(start):AA_only",
  "s(FID)",
  "edf_s(start)",
  "edf_s(start):AA_only",
  "edf_s(FID)",
  "R2",
  "AIC",
  "Deviance_explained",
  "REML",
  "N_cpgs",
  "N_samples",
  "N_train_FIDs",
  "data_chunk_id",
  "region_id",
  "max_diff",
  "mean_diff"
)

stopifnot(length(result_colnames) == 33L)

# Retain only the columns needed for per-split strict filtering and QC.
keep_result_columns <- c(
  "splitID",
  "s(start):AA_only",
  "edf_s(start):AA_only",
  "R2",
  "N_cpgs",
  "N_samples",
  "N_train_FIDs",
  "data_chunk_id",
  "region_id",
  "max_diff",
  "mean_diff"
)

# ------------------------------- input files ---------------------------------

batch1_files <- list.files(
  PATH_results_B1,
  pattern = "^2_run_mgcv_training_splits_results_[0-9]+\\.txt$",
  full.names = TRUE,
  recursive = FALSE
)

batch2_files <- list.files(
  PATH_results_B2,
  pattern = "^3_run_reverse_results_[0-9]+\\.txt$",
  full.names = TRUE,
  recursive = FALSE
)

cat("Forward-batch files:", length(batch1_files), "\n")
cat("Reverse-batch files:", length(batch2_files), "\n")

if (length(batch1_files) + length(batch2_files) == 0L) {
  stop("No MGCV family-training result files were found.")
}

file_manifest <- rbindlist(
  list(
    data.table(
      source_batch = 1L,
      source_label = "forward",
      source_file = batch1_files
    ),
    data.table(
      source_batch = 2L,
      source_label = "reverse",
      source_file = batch2_files
    )
  ),
  use.names = TRUE
)

file_manifest[, file_size_bytes := file.info(source_file)$size]
fwrite(
  file_manifest,
  file.path(qc_output_dir, "input_file_manifest.tsv"),
  sep = "\t"
)

# ----------------------------- regional metadata -----------------------------

if (!file.exists(region_file_path)) {
  stop("Region index does not exist: ", region_file_path)
}

region_meta <- fread(region_file_path)

required_region_columns <- c(
  "region_id",
  "data_chunk_id",
  "chr",
  "region_start",
  "region_end"
)

missing_region_columns <- setdiff(
  required_region_columns,
  names(region_meta)
)

if (length(missing_region_columns) > 0L) {
  stop(
    "Region index is missing columns: ",
    paste(missing_region_columns, collapse = ", ")
  )
}

region_meta <- unique(
  region_meta[, ..required_region_columns],
  by = "region_id"
)

region_meta[
  ,
  `:=`(
    region_id = as.integer(region_id),
    data_chunk_id = as.integer(data_chunk_id),
    chr = as.character(chr),
    region_start = as.integer(region_start),
    region_end = as.integer(region_end)
  )
]

if (anyDuplicated(region_meta$region_id)) {
  stop("region_id is not unique in the region index.")
}

expected_region_count <- nrow(region_meta)
cat("Regions in region index:", expected_region_count, "\n")

# ----------------------------- read result files ------------------------------

read_one_result <- function(path, source_batch, source_label) {
  error_message <- NA_character_

  raw <- tryCatch(
    fread(
      path,
      header = FALSE,
      sep = "\t",
      fill = TRUE,
      showProgress = FALSE
    ),
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(raw)) {
    return(
      list(
        data = NULL,
        problem = data.table(
          source_batch = source_batch,
          source_label = source_label,
          source_file = path,
          problem = error_message
        )
      )
    )
  }

  if (nrow(raw) == 0L) {
    return(
      list(
        data = NULL,
        problem = data.table(
          source_batch = source_batch,
          source_label = source_label,
          source_file = path,
          problem = "empty result file"
        )
      )
    )
  }

  if (ncol(raw) > length(result_colnames)) {
    return(
      list(
        data = NULL,
        problem = data.table(
          source_batch = source_batch,
          source_label = source_label,
          source_file = path,
          problem = paste0(
            "found ", ncol(raw),
            " columns; expected at most ", length(result_colnames)
          )
        )
      )
    )
  }

  if (ncol(raw) < length(result_colnames)) {
    for (j in seq.int(ncol(raw) + 1L, length(result_colnames))) {
      raw[[j]] <- NA_real_
    }
  }

  setnames(raw, result_colnames)
  raw <- raw[, ..keep_result_columns]

  # All model result fields are numeric in the source files.
  raw[
    ,
    (keep_result_columns) := lapply(.SD, function(z) {
      suppressWarnings(as.numeric(z))
    }),
    .SDcols = keep_result_columns
  ]

  raw[
    ,
    `:=`(
      source_batch = as.integer(source_batch),
      source_label = as.character(source_label),
      source_file = basename(path)
    )
  ]

  list(data = raw, problem = NULL)
}

read_results <- vector("list", nrow(file_manifest))
read_problems <- vector("list", nrow(file_manifest))

for (i in seq_len(nrow(file_manifest))) {
  if (i %% 100L == 0L || i == nrow(file_manifest)) {
    cat("Reading file", i, "of", nrow(file_manifest), "\n")
  }

  z <- read_one_result(
    file_manifest$source_file[i],
    file_manifest$source_batch[i],
    file_manifest$source_label[i]
  )

  read_results[[i]] <- z$data
  read_problems[[i]] <- z$problem
}

read_problems_dt <- rbindlist(
  read_problems,
  use.names = TRUE,
  fill = TRUE
)

if (nrow(read_problems_dt) > 0L) {
  fwrite(
    read_problems_dt,
    file.path(qc_output_dir, "input_read_problems.tsv"),
    sep = "\t"
  )
  warning(nrow(read_problems_dt), " input files could not be used.")
}

results_all <- rbindlist(
  read_results,
  use.names = TRUE,
  fill = TRUE
)

rm(read_results)
invisible(gc())

if (nrow(results_all) == 0L) {
  stop("No usable MGCV result rows were read.")
}

cat("Rows read before cleaning:", nrow(results_all), "\n")

# --------------------------- malformed-row handling ---------------------------

results_all[
  ,
  `:=`(
    splitID = as.integer(splitID),
    region_id = as.integer(region_id),
    data_chunk_id = as.integer(data_chunk_id),
    N_cpgs = as.integer(N_cpgs),
    N_samples = as.integer(N_samples),
    N_train_FIDs = as.integer(N_train_FIDs)
  )
]

malformed_rows <- results_all[
  is.na(splitID) |
    is.na(region_id) |
    splitID < 1L |
    splitID > 100L
]

if (nrow(malformed_rows) > 0L) {
  fwrite(
    malformed_rows,
    file.path(qc_output_dir, "malformed_result_rows.tsv"),
    sep = "\t"
  )

  results_all <- results_all[
    !(
      is.na(splitID) |
        is.na(region_id) |
        splitID < 1L |
        splitID > 100L
    )
  ]
}

# A row with more nonmissing model fields is preferred when duplicates exist.
completeness_columns <- c(
  "s(start):AA_only",
  "edf_s(start):AA_only",
  "R2",
  "N_cpgs",
  "N_samples",
  "N_train_FIDs",
  "data_chunk_id",
  "max_diff",
  "mean_diff"
)

results_all[
  ,
  completeness := rowSums(!is.na(.SD)),
  .SDcols = completeness_columns
]

# ----------------------------- duplicate handling -----------------------------

key_counts <- results_all[
  ,
  .N,
  by = .(splitID, region_id)
]

duplicate_keys <- key_counts[N > 1L]

n_duplicate_keys <- nrow(duplicate_keys)
n_extra_duplicate_rows <- sum(duplicate_keys$N - 1L)

cat("Duplicate split-region keys:", n_duplicate_keys, "\n")
cat("Extra duplicate rows:", n_extra_duplicate_rows, "\n")

if (n_duplicate_keys > 0L) {
  duplicate_rows <- results_all[
    duplicate_keys,
    on = .(splitID, region_id),
    nomatch = 0L
  ]

  spread_or_zero <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) return(0)
    diff(range(x))
  }

  duplicate_qc <- duplicate_rows[
    ,
    .(
      n_rows = .N,
      n_batches = uniqueN(source_batch),
      pvalue_spread = spread_or_zero(`s(start):AA_only`),
      edf_spread = spread_or_zero(`edf_s(start):AA_only`),
      R2_spread = spread_or_zero(R2),
      mean_diff_spread = spread_or_zero(mean_diff),
      max_diff_spread = spread_or_zero(max_diff),
      max_completeness = max(completeness, na.rm = TRUE)
    ),
    by = .(splitID, region_id)
  ]

  fwrite(
    duplicate_qc,
    file.path(qc_output_dir, "duplicate_key_qc.tsv"),
    sep = "\t"
  )

  duplicate_conflicts <- duplicate_qc[
    pvalue_spread > 1e-10 |
      edf_spread > 1e-10 |
      R2_spread > 1e-10 |
      mean_diff_spread > 1e-10 |
      max_diff_spread > 1e-10
  ]

  if (nrow(duplicate_conflicts) > 0L) {
    fwrite(
      duplicate_conflicts,
      file.path(qc_output_dir, "duplicate_conflicts.tsv"),
      sep = "\t"
    )

    warning(
      nrow(duplicate_conflicts),
      " duplicate split-region keys contain conflicting values. ",
      "The most complete row is retained, with the forward batch preferred ",
      "when completeness is tied."
    )
  }
}

# Prefer the most complete row. For ties, prefer the forward batch.
setorder(
  results_all,
  splitID,
  region_id,
  -completeness,
  source_batch,
  source_file
)

results_unique <- unique(
  results_all,
  by = c("splitID", "region_id")
)

rm(results_all, key_counts)
invisible(gc())

cat("Unique split-region models:", nrow(results_unique), "\n")

# Add chromosome and parent-region coordinates.
results_unique <- merge(
  results_unique,
  region_meta,
  by = "region_id",
  all.x = TRUE,
  suffixes = c("", "_index"),
  sort = FALSE
)

# ------------------------- per-split FDR and filtering ------------------------

split_ids <- 1:100
strict_list <- vector("list", length(split_ids))
qc_list <- vector("list", length(split_ids))

for (j in seq_along(split_ids)) {
  sid <- split_ids[j]
  split_dt <- copy(results_unique[splitID == sid])

  p_raw <- split_dt[["s(start):AA_only"]]

  valid_positive <- p_raw[
    is.finite(p_raw) &
      p_raw > 0 &
      p_raw <= 1
  ]

  min_positive <- if (length(valid_positive) > 0L) {
    min(valid_positive)
  } else {
    1e-300
  }

  # Match the legacy handling in 24_analyze_BMI.R while retaining pvals_raw.
  p_clean <- p_raw
  p_clean[is.finite(p_clean) & p_clean == 0] <- max(
    min_positive / 2,
    .Machine$double.xmin
  )
  p_clean[is.finite(p_clean) & p_clean < 0] <- max(
    min_positive / 10,
    .Machine$double.xmin
  )
  p_clean[!is.finite(p_clean) | p_clean > 1] <- NA_real_

  split_dt[, pvals_raw := p_raw]
  split_dt[, pvals := p_clean]
  split_dt[, edf := `edf_s(start):AA_only`]
  split_dt[, pval_FDR := p.adjust(pvals, method = "BH")]
  split_dt[, FDR := pval_FDR]

  base_dt <- split_dt[
    !is.na(pvals) &
      pvals < raw_p_threshold &
      !is.na(FDR) &
      FDR < fdr_threshold &
      !is.na(edf) &
      edf >= edf_threshold &
      !is.na(N_cpgs) &
      N_cpgs >= min_cpgs_threshold
  ]

  strict_dt <- base_dt[
    !is.na(R2) &
      R2 >= r2_threshold &
      !is.na(mean_diff) &
      mean_diff >= mean_diff_threshold
  ]

  strict_columns <- c(
    "splitID",
    "region_id",
    "data_chunk_id",
    "chr",
    "region_start",
    "region_end",
    "edf",
    "pvals",
    "pvals_raw",
    "R2",
    "N_cpgs",
    "N_samples",
    "N_train_FIDs",
    "mean_diff",
    "max_diff",
    "pval_FDR",
    "FDR",
    "source_label",
    "source_file"
  )

  strict_columns <- strict_columns[strict_columns %in% names(strict_dt)]
  strict_dt <- strict_dt[, ..strict_columns]
  setorder(strict_dt, pval_FDR, pvals, -mean_diff, region_id)

  strict_file <- file.path(
    strict_output_dir,
    sprintf("24_dmrs_STRICT_split_%03d.tsv", sid)
  )

  fwrite(
    strict_dt,
    strict_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  strict_list[[j]] <- strict_dt

  qc_list[[j]] <- data.table(
    splitID = sid,
    n_models = nrow(split_dt),
    n_unique_regions = uniqueN(split_dt$region_id),
    n_regions_in_index = expected_region_count,
    n_missing_vs_index = expected_region_count - uniqueN(split_dt$region_id),
    n_pvalue_NA_or_gt1 = sum(is.na(p_clean)),
    n_pvalue_zero_raw = sum(p_raw == 0, na.rm = TRUE),
    n_pvalue_negative_raw = sum(p_raw < 0, na.rm = TRUE),
    min_positive_pvalue = min_positive,
    n_raw_p_lt_1e5 = sum(p_clean < raw_p_threshold, na.rm = TRUE),
    n_base = nrow(base_dt),
    n_strict = nrow(strict_dt)
  )

  cat(
    "Split", sid,
    ": models =", nrow(split_dt),
    "base =", nrow(base_dt),
    "strict =", nrow(strict_dt),
    "\n"
  )
}

strict_all <- rbindlist(
  strict_list,
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  strict_all,
  file.path(PATH_output, "24_dmrs_STRICT_all_splits.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

qc_by_split <- rbindlist(qc_list)

fwrite(
  qc_by_split,
  file.path(PATH_output, "mgcv_collection_qc_by_split.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# Selection stability across the 100 family splits.
if (nrow(strict_all) > 0L) {
  strict_stability <- strict_all[
    ,
    .(
      n_splits_STRICT = uniqueN(splitID),
      selection_rate_STRICT = uniqueN(splitID) / 100,
      median_FDR = median(FDR, na.rm = TRUE),
      min_FDR = min(FDR, na.rm = TRUE),
      median_pvals = median(pvals, na.rm = TRUE),
      median_edf = median(edf, na.rm = TRUE),
      median_R2 = median(R2, na.rm = TRUE),
      median_mean_diff = median(mean_diff, na.rm = TRUE),
      median_max_diff = median(max_diff, na.rm = TRUE)
    ),
    by = .(
      region_id,
      data_chunk_id,
      chr,
      region_start,
      region_end
    )
  ]

  setorder(
    strict_stability,
    -n_splits_STRICT,
    median_FDR,
    -median_mean_diff,
    region_id
  )
} else {
  strict_stability <- data.table(
    region_id = integer(),
    data_chunk_id = integer(),
    chr = character(),
    region_start = integer(),
    region_end = integer(),
    n_splits_STRICT = integer(),
    selection_rate_STRICT = numeric(),
    median_FDR = numeric(),
    min_FDR = numeric(),
    median_pvals = numeric(),
    median_edf = numeric(),
    median_R2 = numeric(),
    median_mean_diff = numeric(),
    median_max_diff = numeric()
  )
}

fwrite(
  strict_stability,
  file.path(PATH_output, "24_dmrs_STRICT_stability_across_splits.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

run_summary <- data.table(
  metric = c(
    "forward_files",
    "reverse_files",
    "files_with_read_problems",
    "rows_after_malformed_removal_before_deduplication",
    "duplicate_split_region_keys",
    "extra_duplicate_rows",
    "unique_split_region_models",
    "regions_in_region_index",
    "total_STRICT_region_split_rows",
    "unique_STRICT_regions"
  ),
  value = c(
    length(batch1_files),
    length(batch2_files),
    nrow(read_problems_dt),
    nrow(results_unique) + n_extra_duplicate_rows,
    n_duplicate_keys,
    n_extra_duplicate_rows,
    nrow(results_unique),
    expected_region_count,
    nrow(strict_all),
    uniqueN(strict_all$region_id)
  )
)

fwrite(
  run_summary,
  file.path(PATH_output, "mgcv_collection_run_summary.tsv"),
  sep = "\t"
)

cat("\nCollection finished.\n")
cat("Per-split strict files:", strict_output_dir, "\n")
cat(
  "Combined strict file:",
  file.path(PATH_output, "24_dmrs_STRICT_all_splits.tsv"),
  "\n"
)
cat(
  "QC summary:",
  file.path(PATH_output, "mgcv_collection_qc_by_split.tsv"),
  "\n"
)

#!/usr/bin/env Rscript

# BSmooth on region_BSmooth_M12.csv using the same family-level training
# splits as the MGCV, SOMNiBUS, and DMRcate analyses.

suppressPackageStartupMessages({
  library(bsseq)
  library(BiocParallel)
  library(data.table)
  library(stringr)
  library(GenomicRanges)
})

data.table::setDTthreads(1L)
BiocParallel::register(BiocParallel::SerialParam())
options(mc.cores = 1L)

env_int <- function(name, default) {
  z <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (is.na(z)) stop(name, " must be an integer.")
  z
}

env_num <- function(name, default) {
  z <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
  if (!is.finite(z)) stop(name, " must be finite numeric.")
  z
}

env_bool <- function(name, default = FALSE) {
  Sys.getenv(name, if (default) "1" else "0") %in%
    c("1", "TRUE", "true", "T", "yes", "YES")
}

mem_rss_mb <- function() {
  fp <- "/proc/self/status"
  if (!file.exists(fp)) return(NA_real_)
  x <- grep("^VmRSS:", readLines(fp, warn = FALSE), value = TRUE)
  if (!length(x)) return(NA_real_)
  as.numeric(sub(".*:\\s+([0-9]+).*", "\\1", x[1L])) / 1024
}

build_region_matrices <- function(dt) {
  dt <- data.table::copy(dt)
  meth_cols <- grep("_meth$", names(dt), value = TRUE)
  if (!length(meth_cols)) stop("No *_meth columns found.")

  ids <- stringr::str_remove(meth_cols, "_meth$")
  cov_cols <- paste0(ids, "_tot")
  missing_cov <- setdiff(cov_cols, names(dt))
  if (length(missing_cov)) {
    stop("Missing paired coverage columns: ",
         paste(head(missing_cov, 10L), collapse = ", "))
  }

  dt[, (meth_cols) := lapply(.SD, as.numeric), .SDcols = meth_cols]
  dt[, (cov_cols) := lapply(.SD, as.numeric), .SDcols = cov_cols]
  data.table::setorder(dt, start)

  meth <- as.matrix(dt[, ..meth_cols])
  cov <- as.matrix(dt[, ..cov_cols])
  storage.mode(meth) <- "numeric"
  storage.mode(cov) <- "numeric"

  # Missing methylation is treated as missing coverage, not as methylation = 0.
  bad <- !is.finite(meth) | !is.finite(cov) | cov < 0
  meth[bad] <- 0
  cov[bad] <- 0

  meth <- pmin(pmax(meth, 0), 1)
  cov <- round(pmax(cov, 0))
  M <- round(meth * cov)
  M <- pmin(pmax(M, 0), cov)

  colnames(M) <- ids
  colnames(cov) <- ids

  list(
    M = M,
    Cov = cov,
    sample_names = ids,
    pos = as.integer(dt$start)
  )
}

extract_tstat <- function(x) {
  st <- bsseq::getStats(x)

  if (is.matrix(st) || is.data.frame(st)) {
    st <- as.matrix(st)
    preferred <- c("tstat.corrected", "tstat", "coef", "stat")
    hit <- intersect(preferred, colnames(st))
    if (!length(hit)) {
      stop("Unknown getStats columns: ",
           paste(colnames(st), collapse = ", "))
    }
    return(list(values = as.numeric(st[, hit[1L]]), column = hit[1L]))
  }

  list(values = as.numeric(st), column = "vector")
}

# -------------------------------- paths --------------------------------------

PATH_wk <- path.expand("~/scratch/UQAC/meth/")
PATH_scr11 <- file.path(PATH_wk, "scr/11_mgcv")

region_file_path <- file.path(
  PATH_wk, "results/15_revision/5_cv/region_BSmooth_M12.csv"
)
train_split_path <- file.path(
  PATH_scr11, "dat/100_family_training_splits.csv"
)
pheno_file_path <- file.path(
  PATH_scr11, "dat/bsmooth_pheno_only.RData"
)
chunk_dir <- file.path(PATH_wk, "data/meth_split")
PATH_output <- file.path(
  PATH_wk, "results/15_revision/5_cv/4_bsmooth/training_splits"
)
dir.create(PATH_output, recursive = TRUE, showWarnings = FALSE)

# ------------------------------ configuration -------------------------------

min_cpgs <- env_int("BSMOOTH_MIN_CPGS", 10L)
smooth_ns <- env_int("BSMOOTH_NS", 25L)
smooth_h <- env_num("BSMOOTH_H", 400)
smooth_max_gap <- env_num("BSMOOTH_MAX_GAP", 5e7)
tstat_cutoff <- env_num("BSMOOTH_TSTAT_CUTOFF", 4.417)
min_group_samples <- env_int("BSMOOTH_MIN_GROUP_SAMPLES", 2L)
estimate_var <- Sys.getenv("BSMOOTH_ESTIMATE_VAR", "group2")
local_correct <- env_bool("BSMOOTH_LOCAL_CORRECT", TRUE)

N_jobs_requested <- env_int("N_JOBS", 990L)
split_start <- env_int("SPLIT_START", 1L)
split_end <- env_int("SPLIT_END", 100L)
max_regions <- env_int("MAX_REGIONS", 0L)
max_runtime_minutes <- env_num("MAX_RUNTIME_MINUTES", 165)
overwrite <- env_bool("OVERWRITE", FALSE)
retry_errors <- env_bool("RETRY_ERRORS", FALSE)

if (min_cpgs < 2L) stop("BSMOOTH_MIN_CPGS must be at least 2.")
if (smooth_ns < 1L) stop("BSMOOTH_NS must be positive.")
if (smooth_h <= 0 || smooth_max_gap <= 0) {
  stop("BSMOOTH_H and BSMOOTH_MAX_GAP must be positive.")
}
if (tstat_cutoff <= 0) stop("BSMOOTH_TSTAT_CUTOFF must be positive.")
if (min_group_samples < 2L) {
  stop("BSMOOTH_MIN_GROUP_SAMPLES must be at least 2.")
}
if (N_jobs_requested < 1L) stop("N_JOBS must be positive.")
if (split_start < 1L || split_end < split_start) {
  stop("Invalid SPLIT_START/SPLIT_END.")
}

# -------------------------------- arguments ---------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Pass SLURM_ARRAY_TASK_ID.")
i_job <- suppressWarnings(as.integer(args[1L]))
if (is.na(i_job) || i_job < 1L) {
  stop("SLURM_ARRAY_TASK_ID must be a positive integer.")
}
cat("BSmooth family-training job ID:", i_job, "\n")

# ----------------------------- selected regions ------------------------------

if (!file.exists(region_file_path)) {
  stop("Region file does not exist: ", region_file_path)
}
region_file <- data.table::fread(region_file_path)

region_cols <- c(
  "data_chunk_id", "chr", "region_start",
  "region_end", "n_cpgs", "region_id"
)
missing_cols <- setdiff(region_cols, names(region_file))
if (length(missing_cols)) {
  stop("Region file is missing: ", paste(missing_cols, collapse = ", "))
}

region_file <- region_file[complete.cases(region_file[, ..region_cols])]
region_file[, (region_cols) := lapply(.SD, as.integer), .SDcols = region_cols]

if (!nrow(region_file)) stop("region_BSmooth_M12.csv has no complete rows.")
if (anyDuplicated(region_file$region_id)) stop("region_id must be unique.")
if (any(region_file$region_start > region_file$region_end)) {
  stop("At least one region has start > end.")
}

data.table::setorder(
  region_file, data_chunk_id, region_start, region_end, region_id
)

N_jobs_effective <- min(N_jobs_requested, nrow(region_file))
job_id_for_row <- ceiling(
  seq_len(nrow(region_file)) * N_jobs_effective / nrow(region_file)
)
region_row_list <- split(seq_len(nrow(region_file)), job_id_for_row)

if (i_job > N_jobs_effective) {
  cat("No selected regions assigned to this task.\n")
  quit(save = "no", status = 0L)
}

rows_this_job <- region_row_list[[as.character(i_job)]]
if (is.null(rows_this_job) || !length(rows_this_job)) {
  cat("No selected regions assigned to this task.\n")
  quit(save = "no", status = 0L)
}

# Keep original region_id values; never replace them with row numbers.
regions_this_job <- region_file[rows_this_job]
if (max_regions > 0L) regions_this_job <- head(regions_this_job, max_regions)

cat(
  "Selected regions in file:", nrow(region_file), "\n",
  "Requested jobs:", N_jobs_requested, "\n",
  "Effective jobs:", N_jobs_effective, "\n",
  "Regions selected for this run:", nrow(regions_this_job), "\n"
)

# ----------------------------- family split file -----------------------------

if (!file.exists(train_split_path)) {
  stop("Training split file does not exist: ", train_split_path)
}
train_df <- data.table::fread(train_split_path)
if (!"splitID" %in% names(train_df)) stop("splitID column is missing.")

train_cols <- grep("^train_FID_[0-9]+$", names(train_df), value = TRUE)
if (!length(train_cols)) stop("No train_FID_* columns found.")
train_cols <- train_cols[
  order(as.integer(sub("^train_FID_", "", train_cols)))
]

train_df[, splitID := as.integer(splitID)]
if (anyNA(train_df$splitID)) stop("splitID contains NA/nonnumeric values.")
if (anyDuplicated(train_df$splitID)) stop("splitID must be unique.")

split_ids <- train_df[
  splitID >= split_start & splitID <= split_end,
  splitID
]
if (!length(split_ids)) stop("No splitIDs selected.")

train_FID_list <- setNames(
  lapply(split_ids, function(sid) {
    z <- unlist(train_df[splitID == sid, ..train_cols], use.names = FALSE)
    z <- trimws(as.character(z))
    unique(z[!is.na(z) & nzchar(z) & z != "NA"])
  }),
  as.character(split_ids)
)

cat(
  "Splits selected for this run:", length(split_ids), "\n",
  "Training-FID columns per split:", length(train_cols), "\n"
)

# -------------------------------- phenotype ---------------------------------

if (!file.exists(pheno_file_path)) {
  stop("Phenotype file does not exist: ", pheno_file_path)
}
e <- new.env(parent = emptyenv())
load(pheno_file_path, envir = e, verbose = TRUE)
if (!exists("pheno_file", envir = e, inherits = FALSE)) {
  stop("Phenotype RData does not contain pheno_file.")
}
pheno_file <- data.table::as.data.table(
  data.table::copy(get("pheno_file", envir = e, inherits = FALSE))
)
rm(e)

pheno_cols <- c("ID", "FID", "AA_only")
missing_pheno <- setdiff(pheno_cols, names(pheno_file))
if (length(missing_pheno)) {
  stop("pheno_file is missing: ", paste(missing_pheno, collapse = ", "))
}

pheno_file <- pheno_file[, ..pheno_cols]
pheno_file[, ID := as.character(ID)]
pheno_file[, FID := as.character(FID)]
pheno_file[, AA_only := suppressWarnings(as.integer(as.character(AA_only)))]

if (any(!is.na(pheno_file$AA_only) &
        !pheno_file$AA_only %in% c(0L, 1L))) {
  stop("AA_only must contain only 0, 1, or NA.")
}
if (anyDuplicated(pheno_file$ID)) stop("Duplicated IDs in pheno_file.")


# -------------------------------- outputs -----------------------------------

summary_out <- file.path(
  PATH_output,
  sprintf("bsmooth_training_summary_job_%04d.tsv", i_job)
)
dmr_out <- file.path(
  PATH_output,
  sprintf("bsmooth_training_dmrs_job_%04d.tsv", i_job)
)

if (overwrite) unlink(c(summary_out, dmr_out), force = TRUE)

error_statuses <- c(
  "MATRIX_ERROR", "BSSEQ_ERROR", "BSMOOTH_ERROR",
  "TSTAT_ERROR", "TSTAT_INVALID", "DMRFINDER_ERROR"
)

completed_keys <- character()

if (file.exists(summary_out)) {
  previous <- data.table::fread(summary_out)

  if (!all(c("region_id", "splitID", "status") %in% names(previous))) {
    stop("Existing summary has incompatible columns: ", summary_out)
  }

  if (retry_errors) {
    retry_rows <- previous$status %in% error_statuses

    if (any(retry_rows)) {
      previous <- previous[!retry_rows]

      if (!nrow(previous)) {
        unlink(summary_out, force = TRUE)
      } else {
        data.table::fwrite(
          previous, summary_out,
          sep = "\t", quote = FALSE, na = "NA"
        )
      }
    }
  }

  completed_keys <- unique(
    paste(previous$region_id, previous$splitID, sep = "::")
  )
  cat("Existing completed region-split rows:",
      length(completed_keys), "\n")
}

existing_dmr_keys <- character()

if (file.exists(dmr_out)) {
  old_dmrs <- tryCatch(data.table::fread(dmr_out), error = function(e) NULL)

  if (!is.null(old_dmrs) &&
      all(c("region_id", "splitID") %in% names(old_dmrs))) {
    existing_dmr_keys <- unique(
      paste(old_dmrs$region_id, old_dmrs$splitID, sep = "::")
    )
  }
}

append_summary <- function(
  info, split_id, status, message = NA_character_,
  n_requested_fids = NA_integer_,
  n_observed_fids = NA_integer_,
  n_samples = NA_integer_,
  n_case = NA_integer_,
  n_ctrl = NA_integer_,
  n_cpgs = NA_integer_,
  n_finite_stats = NA_integer_,
  stat_column = NA_character_,
  abs_t_q95 = NA_real_,
  max_abs_t = NA_real_,
  n_dmrs = NA_integer_,
  elapsed_seconds = NA_real_
) {
  row <- data.table::data.table(
    prep_job_id = i_job,
    region_id = as.integer(info$region_id),
    data_chunk_id = as.integer(info$data_chunk_id),
    chr = as.character(info$chr),
    region_start = as.integer(info$region_start),
    region_end = as.integer(info$region_end),
    splitID = as.integer(split_id),
    status = as.character(status),
    message = as.character(message),
    n_train_FIDs_requested = as.integer(n_requested_fids),
    n_train_FIDs_observed = as.integer(n_observed_fids),
    n_train_samples = as.integer(n_samples),
    n_train_cases = as.integer(n_case),
    n_train_controls = as.integer(n_ctrl),
    n_cpgs = as.integer(n_cpgs),
    n_finite_stats = as.integer(n_finite_stats),
    tstat_column = as.character(stat_column),
    abs_t_q95 = as.numeric(abs_t_q95),
    max_abs_t = as.numeric(max_abs_t),
    n_dmrs = as.integer(n_dmrs),
    fit_seconds = as.numeric(elapsed_seconds)
  )

  data.table::fwrite(
    row,
    file = summary_out,
    sep = "\t",
    append = file.exists(summary_out),
    col.names = !file.exists(summary_out),
    quote = FALSE,
    na = "NA"
  )
}

cat("Summary file:", summary_out, "\n")
cat("DMR file:", dmr_out, "\n")
cat("Initial RSS (MB):", round(mem_rss_mb(), 1), "\n")

# -------------------------------- main loop ---------------------------------

run_start <- Sys.time()
stop_requested <- FALSE
chunks <- split(regions_this_job, regions_this_job$data_chunk_id)

for (chunk_name in names(chunks)) {
  if (stop_requested) break

  chunk_regions <- data.table::as.data.table(chunks[[chunk_name]])
  chunk_id <- as.integer(chunk_name)
  chunk_fp <- file.path(chunk_dir, sprintf("chunk_%04d.csv", chunk_id))

  if (!file.exists(chunk_fp)) {
    for (row_i in seq_len(nrow(chunk_regions))) {
      info <- chunk_regions[row_i]

      for (split_id in split_ids) {
        key <- paste(info$region_id, split_id, sep = "::")
        if (key %in% completed_keys) next

        append_summary(
          info, split_id,
          status = "MISSING_CHUNK_FILE",
          message = chunk_fp,
          n_requested_fids = length(
            train_FID_list[[as.character(split_id)]]
          )
        )
        completed_keys <- c(completed_keys, key)
      }
    }
    next
  }

  header <- data.table::fread(chunk_fp, nrows = 0L, showProgress = FALSE)
  all_cols <- names(header)

  raw_ids <- unique(
    stringr::str_remove(
      setdiff(all_cols, c("chr", "start", "end")),
      "_meth$|_tot$"
    )
  )

  paired_ids <- raw_ids[
    paste0(raw_ids, "_meth") %in% all_cols &
      paste0(raw_ids, "_tot") %in% all_cols
  ]
  phenotype_ids <- intersect(paired_ids, pheno_file$ID)

  if (!length(phenotype_ids)) {
    for (row_i in seq_len(nrow(chunk_regions))) {
      info <- chunk_regions[row_i]

      for (split_id in split_ids) {
        key <- paste(info$region_id, split_id, sep = "::")
        if (key %in% completed_keys) next

        append_summary(
          info, split_id,
          status = "NO_PHENOTYPE_SAMPLES",
          n_requested_fids = length(
            train_FID_list[[as.character(split_id)]]
          )
        )
        completed_keys <- c(completed_keys, key)
      }
    }
    next
  }

  selected_cols <- c(
    "chr", "start", "end",
    paste0(phenotype_ids, "_meth"),
    paste0(phenotype_ids, "_tot")
  )

  cat(sprintf(
    "Reading chunk %d with %d selected columns (%d paired samples)\n",
    chunk_id, length(selected_cols), length(phenotype_ids)
  ))

  chunk_data <- data.table::fread(
    chunk_fp, select = selected_cols, showProgress = FALSE
  )
  chunk_data[, chr := as.character(chr)]
  chunk_data[, start := as.integer(start)]

  for (row_i in seq_len(nrow(chunk_regions))) {
    if (stop_requested) break

    info <- chunk_regions[row_i]
    info[, chr := as.character(chr)]

    cat("\nRegion", info$region_id, "chunk", info$data_chunk_id, "\n")

    region_all <- chunk_data[
      chr == info$chr &
        start >= info$region_start &
        start <= info$region_end
    ]
    n_cpgs_region <- nrow(region_all)

    if (n_cpgs_region < min_cpgs) {
      for (split_id in split_ids) {
        key <- paste(info$region_id, split_id, sep = "::")
        if (key %in% completed_keys) next

        append_summary(
          info, split_id,
          status = "BELOW_MIN_CPGS",
          message = paste0(
            "Observed CpGs = ", n_cpgs_region,
            "; required = ", min_cpgs
          ),
          n_requested_fids = length(
            train_FID_list[[as.character(split_id)]]
          ),
          n_cpgs = n_cpgs_region
        )
        completed_keys <- c(completed_keys, key)
      }
      next
    }

    matrices <- tryCatch(
      build_region_matrices(region_all),
      error = function(e) e
    )

    if (inherits(matrices, "error")) {
      for (split_id in split_ids) {
        key <- paste(info$region_id, split_id, sep = "::")
        if (key %in% completed_keys) next

        append_summary(
          info, split_id,
          status = "MATRIX_ERROR",
          message = conditionMessage(matrices),
          n_requested_fids = length(
            train_FID_list[[as.character(split_id)]]
          ),
          n_cpgs = n_cpgs_region
        )
        completed_keys <- c(completed_keys, key)
      }
      next
    }

    for (split_id in split_ids) {
      elapsed_minutes <- as.numeric(
        difftime(Sys.time(), run_start, units = "mins")
      )

      if (max_runtime_minutes > 0 &&
          elapsed_minutes >= max_runtime_minutes) {
        cat("Reached MAX_RUNTIME_MINUTES =",
            max_runtime_minutes, "; stopping cleanly.\n")
        stop_requested <- TRUE
        break
      }

      key <- paste(info$region_id, split_id, sep = "::")
      if (key %in% completed_keys) next

      fit_start <- proc.time()[["elapsed"]]
      train_fids <- train_FID_list[[as.character(split_id)]]
      n_requested_fids <- length(train_fids)

      if (!n_requested_fids) {
        append_summary(
          info, split_id,
          status = "NO_TRAIN_FIDS",
          n_requested_fids = 0L,
          n_cpgs = n_cpgs_region
        )
        completed_keys <- c(completed_keys, key)
        next
      }

      ph_train <- pheno_file[
        FID %chin% train_fids & !is.na(AA_only)
      ]
      ph_train <- ph_train[
        ID %chin% matrices$sample_names
      ]

      if (!nrow(ph_train)) {
        append_summary(
          info, split_id,
          status = "NO_TRAINING_SAMPLES",
          n_requested_fids = n_requested_fids,
          n_cpgs = n_cpgs_region
        )
        completed_keys <- c(completed_keys, key)
        next
      }

      sample_index <- match(ph_train$ID, matrices$sample_names)
      keep_match <- !is.na(sample_index)
      ph_sub <- ph_train[keep_match]
      sample_index <- sample_index[keep_match]

      M2 <- matrices$M[, sample_index, drop = FALSE]
      Cov2 <- matrices$Cov[, sample_index, drop = FALSE]

      # Remove samples with no coverage anywhere in this region.
      has_coverage <- colSums(Cov2, na.rm = TRUE) > 0
      M2 <- M2[, has_coverage, drop = FALSE]
      Cov2 <- Cov2[, has_coverage, drop = FALSE]
      ph_sub <- ph_sub[has_coverage]

      sample_ids <- ph_sub$ID
      colnames(M2) <- sample_ids
      colnames(Cov2) <- sample_ids

      n_samples <- ncol(M2)
      n_observed_fids <- data.table::uniqueN(ph_sub$FID)
      grp <- as.integer(ph_sub$AA_only)
      n_case <- sum(grp == 1L)
      n_ctrl <- sum(grp == 0L)

      if (n_samples < 2L) {
        append_summary(
          info, split_id,
          status = "TOO_FEW_TRAINING_SAMPLES",
          n_requested_fids = n_requested_fids,
          n_observed_fids = n_observed_fids,
          n_samples = n_samples,
          n_case = n_case,
          n_ctrl = n_ctrl,
          n_cpgs = n_cpgs_region,
          elapsed_seconds = proc.time()[["elapsed"]] - fit_start
        )
        completed_keys <- c(completed_keys, key)
        next
      }

      if (n_case < min_group_samples || n_ctrl < min_group_samples) {
        append_summary(
          info, split_id,
          status = "INSUFFICIENT_GROUP_SAMPLES",
          message = paste0(
            "Cases = ", n_case,
            "; controls = ", n_ctrl,
            "; required per group = ", min_group_samples
          ),
          n_requested_fids = n_requested_fids,
          n_observed_fids = n_observed_fids,
          n_samples = n_samples,
          n_case = n_case,
          n_ctrl = n_ctrl,
          n_cpgs = n_cpgs_region,
          elapsed_seconds = proc.time()[["elapsed"]] - fit_start
        )
        completed_keys <- c(completed_keys, key)
        next
      }


      BS <- tryCatch(
        bsseq::BSseq(
          M = M2,
          Cov = Cov2,
          chr = rep(as.character(info$chr), nrow(M2)),
          pos = matrices$pos,
          sampleNames = sample_ids
        ),
        error = function(e) e
      )

      if (inherits(BS, "error")) {
        append_summary(
          info, split_id,
          status = "BSSEQ_ERROR",
          message = conditionMessage(BS),
          n_requested_fids = n_requested_fids,
          n_observed_fids = n_observed_fids,
          n_samples = n_samples,
          n_case = n_case,
          n_ctrl = n_ctrl,
          n_cpgs = n_cpgs_region,
          elapsed_seconds = proc.time()[["elapsed"]] - fit_start
        )
        completed_keys <- c(completed_keys, key)
        next
      }

      BS_sm <- tryCatch(
        bsseq::BSmooth(
          BS,
          ns = smooth_ns,
          h = smooth_h,
          maxGap = smooth_max_gap,
          verbose = FALSE,
          BPPARAM = BiocParallel::SerialParam()
        ),
        error = function(e) e
      )

      if (inherits(BS_sm, "error")) {
        append_summary(
          info, split_id,
          status = "BSMOOTH_ERROR",
          message = conditionMessage(BS_sm),
          n_requested_fids = n_requested_fids,
          n_observed_fids = n_observed_fids,
          n_samples = n_samples,
          n_case = n_case,
          n_ctrl = n_ctrl,
          n_cpgs = n_cpgs_region,
          elapsed_seconds = proc.time()[["elapsed"]] - fit_start
        )
        completed_keys <- c(completed_keys, key)
        rm(BS)
        invisible(gc(FALSE))
        next
      }

      group1 <- which(grp == 1L)
      group2 <- which(grp == 0L)

      tstat <- tryCatch(
        bsseq::BSmooth.tstat(
          BS_sm,
          group1 = group1,
          group2 = group2,
          estimate.var = estimate_var,
          local.correct = local_correct,
          verbose = FALSE
        ),
        error = function(e) e
      )

      if (inherits(tstat, "error")) {
        append_summary(
          info, split_id,
          status = "TSTAT_ERROR",
          message = conditionMessage(tstat),
          n_requested_fids = n_requested_fids,
          n_observed_fids = n_observed_fids,
          n_samples = n_samples,
          n_case = n_case,
          n_ctrl = n_ctrl,
          n_cpgs = n_cpgs_region,
          elapsed_seconds = proc.time()[["elapsed"]] - fit_start
        )
        completed_keys <- c(completed_keys, key)
        rm(BS, BS_sm)
        invisible(gc(FALSE))
        next
      }

      stat_result <- tryCatch(
        extract_tstat(tstat),
        error = function(e) e
      )

      if (inherits(stat_result, "error")) {
        append_summary(
          info, split_id,
          status = "TSTAT_INVALID",
          message = conditionMessage(stat_result),
          n_requested_fids = n_requested_fids,
          n_observed_fids = n_observed_fids,
          n_samples = n_samples,
          n_case = n_case,
          n_ctrl = n_ctrl,
          n_cpgs = n_cpgs_region,
          elapsed_seconds = proc.time()[["elapsed"]] - fit_start
        )
        completed_keys <- c(completed_keys, key)
        rm(BS, BS_sm, tstat)
        invisible(gc(FALSE))
        next
      }

      stat_vec <- stat_result$values
      stat_column <- stat_result$column
      n_granges <- length(GenomicRanges::granges(tstat))

      if (!length(stat_vec) ||
          length(stat_vec) != n_granges ||
          !any(is.finite(stat_vec))) {
        append_summary(
          info, split_id,
          status = "TSTAT_INVALID",
          message = paste0(
            "stat length = ", length(stat_vec),
            "; granges length = ", n_granges,
            "; finite = ", sum(is.finite(stat_vec))
          ),
          n_requested_fids = n_requested_fids,
          n_observed_fids = n_observed_fids,
          n_samples = n_samples,
          n_case = n_case,
          n_ctrl = n_ctrl,
          n_cpgs = n_cpgs_region,
          n_finite_stats = sum(is.finite(stat_vec)),
          stat_column = stat_column,
          elapsed_seconds = proc.time()[["elapsed"]] - fit_start
        )
        completed_keys <- c(completed_keys, key)
        rm(BS, BS_sm, tstat, stat_result, stat_vec)
        invisible(gc(FALSE))
        next
      }

      finite_abs_t <- abs(stat_vec[is.finite(stat_vec)])
      abs_t_q95 <- as.numeric(
        stats::quantile(
          finite_abs_t, 0.95,
          na.rm = TRUE, names = FALSE
        )
      )
      max_abs_t <- max(finite_abs_t, na.rm = TRUE)

      dmrs <- tryCatch(
        bsseq::dmrFinder(
          tstat,
          cutoff = c(-abs(tstat_cutoff), abs(tstat_cutoff))
        ),
        error = function(e) e
      )

      if (inherits(dmrs, "error")) {
        append_summary(
          info, split_id,
          status = "DMRFINDER_ERROR",
          message = conditionMessage(dmrs),
          n_requested_fids = n_requested_fids,
          n_observed_fids = n_observed_fids,
          n_samples = n_samples,
          n_case = n_case,
          n_ctrl = n_ctrl,
          n_cpgs = n_cpgs_region,
          n_finite_stats = length(finite_abs_t),
          stat_column = stat_column,
          abs_t_q95 = abs_t_q95,
          max_abs_t = max_abs_t,
          elapsed_seconds = proc.time()[["elapsed"]] - fit_start
        )
        completed_keys <- c(completed_keys, key)
        rm(BS, BS_sm, tstat, stat_result, stat_vec, finite_abs_t)
        invisible(gc(FALSE))
        next
      }

      n_dmrs <- if (is.null(dmrs) || NROW(dmrs) == 0L) 0L else NROW(dmrs)
      fit_seconds <- proc.time()[["elapsed"]] - fit_start

      if (n_dmrs > 0L && !key %in% existing_dmr_keys) {
        dmr_table <- data.table::as.data.table(dmrs)

        dmr_table[, `:=`(
          prep_job_id = i_job,
          region_id = as.integer(info$region_id),
          data_chunk_id = as.integer(info$data_chunk_id),
          parent_chr = as.character(info$chr),
          parent_start = as.integer(info$region_start),
          parent_end = as.integer(info$region_end),
          splitID = as.integer(split_id),
          n_train_FIDs_requested = as.integer(n_requested_fids),
          n_train_FIDs_observed = as.integer(n_observed_fids),
          n_train_samples = as.integer(n_samples),
          n_train_cases = as.integer(n_case),
          n_train_controls = as.integer(n_ctrl),
          n_cpgs_parent = as.integer(n_cpgs_region),
          tstat_column = stat_column,
          tstat_cutoff = as.numeric(tstat_cutoff)
        )]

        if ("areaStat" %in% names(dmr_table)) {
          dmr_table[, score := abs(areaStat)]
          dmr_table[, direction := sign(areaStat)]
        } else if ("maxStat" %in% names(dmr_table)) {
          dmr_table[, score := abs(maxStat)]
          dmr_table[, direction := sign(maxStat)]
        } else {
          dmr_table[, score := NA_real_]
          dmr_table[, direction := NA_integer_]
        }

        if (all(c("start", "end") %in% names(dmr_table))) {
          dmr_table[, midpoint := floor((start + end) / 2)]
          if (!"width" %in% names(dmr_table)) {
            dmr_table[, width := end - start + 1L]
          }
        }

        dmr_table[, dmr_index := seq_len(.N)]

        data.table::fwrite(
          dmr_table,
          file = dmr_out,
          sep = "\t",
          append = file.exists(dmr_out),
          col.names = !file.exists(dmr_out),
          quote = FALSE,
          na = "NA"
        )
        existing_dmr_keys <- c(existing_dmr_keys, key)
      }

      append_summary(
        info, split_id,
        status = "FIT_OK",
        n_requested_fids = n_requested_fids,
        n_observed_fids = n_observed_fids,
        n_samples = n_samples,
        n_case = n_case,
        n_ctrl = n_ctrl,
        n_cpgs = n_cpgs_region,
        n_finite_stats = length(finite_abs_t),
        stat_column = stat_column,
        abs_t_q95 = abs_t_q95,
        max_abs_t = max_abs_t,
        n_dmrs = n_dmrs,
        elapsed_seconds = fit_seconds
      )

      completed_keys <- c(completed_keys, key)

      cat(
        "  completed split", split_id,
        "DMRs =", n_dmrs,
        "CpGs =", n_cpgs_region,
        "samples =", n_samples,
        "FIDs =", n_observed_fids,
        "|t|max =", signif(max_abs_t, 4),
        "seconds =", round(fit_seconds, 2),
        "\n"
      )

      rm(
        M2, Cov2, BS, BS_sm, tstat,
        stat_result, stat_vec, finite_abs_t, dmrs
      )
      invisible(gc(FALSE))
    }

    rm(region_all, matrices)
    invisible(gc(FALSE))
  }

  rm(
    chunk_data, header, all_cols, raw_ids,
    paired_ids, phenotype_ids, selected_cols, chunk_regions
  )
  invisible(gc(FALSE))
}

cat("Finished BSmooth family-level training analysis.\n")
cat("Summary:", summary_out, "\n")
cat("DMRs:", dmr_out, "\n")
cat("Final RSS (MB):", round(mem_rss_mb(), 1), "\n")

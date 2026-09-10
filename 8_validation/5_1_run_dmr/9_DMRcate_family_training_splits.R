#!/usr/bin/env Rscript

# Run DMRcate on selected regions using family-level training splits.
# Every model uses training families only. Original region_id values are kept.

safe_extract_gr <- function(
  dmrc_obj,
  cutoff_fdr = 0.01,
  min_cpgs = 10L,
  min_abs_md = 0.05,
  min_width = 250L
) {
  gr_try <- try(
    DMRcate::extractRanges(
      dmrc_obj,
      genome = "hg19",
      cutoff = cutoff_fdr,
      min.cpgs = min_cpgs
    ),
    silent = TRUE
  )

  if (!inherits(gr_try, "try-error")) {
    if (!is.null(S4Vectors::mcols(gr_try)$meandiff)) {
      gr_try <- gr_try[
        abs(S4Vectors::mcols(gr_try)$meandiff) >= min_abs_md
      ]
    }
    gr_try <- gr_try[IRanges::width(gr_try) >= min_width]
    return(gr_try)
  }

  if (methods::is(dmrc_obj, "DMResults")) {
    coords <- methods::slot(dmrc_obj, "coord")
    no.cpgs <- methods::slot(dmrc_obj, "no.cpgs")
    minfdr <- methods::slot(dmrc_obj, "min_smoothed_fdr")

    get_slot_or_na <- function(object, slot_name) {
      if (slot_name %in% methods::slotNames(object)) {
        methods::slot(object, slot_name)
      } else {
        NA_real_
      }
    }

    stouffer <- get_slot_or_na(dmrc_obj, "Stouffer")
    hmfdr <- get_slot_or_na(dmrc_obj, "HMFDR")
    fisher <- get_slot_or_na(dmrc_obj, "Fisher")
    maxdiff <- get_slot_or_na(dmrc_obj, "maxdiff")
    meandiff <- get_slot_or_na(dmrc_obj, "meandiff")

    parse_coord <- function(z) {
      m <- regexec("^([^:]+):(\\d+)-(\\d+)$", z)
      r <- regmatches(z, m)[[1L]]
      if (length(r) != 4L) {
        return(c(NA_character_, NA_character_, NA_character_))
      }
      c(r[2L], r[3L], r[4L])
    }

    parsed <- t(
      vapply(
        as.character(coords),
        parse_coord,
        FUN.VALUE = character(3L)
      )
    )

    gr <- GenomicRanges::GRanges(
      seqnames = parsed[, 1L],
      ranges = IRanges::IRanges(
        start = suppressWarnings(as.integer(parsed[, 2L])),
        end = suppressWarnings(as.integer(parsed[, 3L]))
      )
    )

    S4Vectors::mcols(gr)$no.cpgs <- as.integer(no.cpgs)
    S4Vectors::mcols(gr)$minfdr <- as.numeric(minfdr)
    S4Vectors::mcols(gr)$Stouffer <- as.numeric(stouffer)
    S4Vectors::mcols(gr)$HMFDR <- as.numeric(hmfdr)
    S4Vectors::mcols(gr)$Fisher <- as.numeric(fisher)
    S4Vectors::mcols(gr)$maxdiff <- as.numeric(maxdiff)
    S4Vectors::mcols(gr)$meandiff <- as.numeric(meandiff)

    keep <- (
      !is.na(S4Vectors::mcols(gr)$minfdr) &
        S4Vectors::mcols(gr)$minfdr <= cutoff_fdr &
        S4Vectors::mcols(gr)$no.cpgs >= min_cpgs &
        IRanges::width(gr) >= min_width
    )

    if (!all(is.na(S4Vectors::mcols(gr)$meandiff))) {
      keep <- keep &
        abs(S4Vectors::mcols(gr)$meandiff) >= min_abs_md
    }

    return(gr[keep])
  }

  GenomicRanges::GRanges()
}

suppressPackageStartupMessages({
  library(DMRcate)
  library(BiocParallel)
  library(bsseq)
  library(edgeR)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(GenomicRanges)
})

data.table::setDTthreads(1L)
BiocParallel::register(BiocParallel::SerialParam())
options(mc.cores = 1L)

# ----------------------------- configuration ---------------------------------
PATH_wk <- path.expand("~/scratch/UQAC/meth/")
PATH_scr11 <- file.path(PATH_wk, "scr/11_mgcv")

region_file_path <- file.path(
  PATH_wk,
  "results/15_revision/5_cv/region_DMRcate_M12.csv"
)

train_split_path <- file.path(
  PATH_scr11,
  "dat/100_family_training_splits.csv"
)

pheno_file_path <- file.path(
  PATH_scr11,
  "dat/bsmooth_pheno_only.RData"
)

chunk_dir <- file.path(PATH_wk, "data/meth_split")

PATH_out <- file.path(
  PATH_wk,
  "results/15_revision/5_cv/3_dmrcate/training_splits"
)

dir.create(PATH_out, recursive = TRUE, showWarnings = FALSE)

min_cpgs <- as.integer(Sys.getenv("DMRCATE_MIN_CPGS", "10"))
min_cov <- as.numeric(Sys.getenv("DMRCATE_MIN_COV", "5"))
coverage_fraction <- as.numeric(
  Sys.getenv("DMRCATE_COVERAGE_FRACTION", "0.70")
)
seed_fdr <- as.numeric(Sys.getenv("DMRCATE_SEED_FDR", "0.05"))
min_seed_cpgs <- as.integer(
  Sys.getenv("DMRCATE_MIN_SEED_CPGS", "5")
)
extract_fdr <- as.numeric(
  Sys.getenv("DMRCATE_EXTRACT_FDR", "0.01")
)
min_abs_md <- as.numeric(
  Sys.getenv("DMRCATE_MIN_ABS_MD", "0.05")
)
min_width <- as.integer(
  Sys.getenv("DMRCATE_MIN_WIDTH", "250")
)
lambda_value <- as.numeric(Sys.getenv("DMRCATE_LAMBDA", "1000"))
C_value <- as.numeric(Sys.getenv("DMRCATE_C", "2"))

N_jobs_requested <- as.integer(Sys.getenv("N_JOBS", "900"))
split_start <- as.integer(Sys.getenv("SPLIT_START", "1"))
split_end <- as.integer(Sys.getenv("SPLIT_END", "100"))
max_regions <- as.integer(Sys.getenv("MAX_REGIONS", "0"))
max_runtime_minutes <- as.numeric(
  Sys.getenv("MAX_RUNTIME_MINUTES", "1380")
)
overwrite <- identical(Sys.getenv("OVERWRITE", "0"), "1")
retry_errors <- identical(Sys.getenv("RETRY_ERRORS", "0"), "1")

if (
  is.na(N_jobs_requested) ||
    N_jobs_requested < 1L
) {
  stop("N_JOBS must be a positive integer.")
}

if (
  is.na(split_start) ||
    is.na(split_end) ||
    split_start < 1L ||
    split_end < split_start
) {
  stop("Invalid SPLIT_START/SPLIT_END values.")
}

if (
  !is.finite(coverage_fraction) ||
    coverage_fraction <= 0 ||
    coverage_fraction > 1
) {
  stop("DMRCATE_COVERAGE_FRACTION must be in (0, 1].")
}

# -------------------------------- arguments ----------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1L) {
  stop("Pass SLURM_ARRAY_TASK_ID.")
}

i_job <- suppressWarnings(as.integer(args[1L]))

if (is.na(i_job) || i_job < 1L) {
  stop("SLURM_ARRAY_TASK_ID must be a positive integer.")
}

cat("DMRcate family-training job ID:", i_job, "\n")

# ------------------------------- region file ---------------------------------
if (!file.exists(region_file_path)) {
  stop("Region file does not exist: ", region_file_path)
}

region_file <- data.table::fread(region_file_path)

required_region_columns <- c(
  "data_chunk_id",
  "chr",
  "region_start",
  "region_end",
  "n_cpgs",
  "region_id"
)

missing_region_columns <- setdiff(
  required_region_columns,
  names(region_file)
)

if (length(missing_region_columns) > 0L) {
  stop(
    "Region file is missing columns: ",
    paste(missing_region_columns, collapse = ", ")
  )
}

region_file <- region_file[
  complete.cases(region_file[, ..required_region_columns])
]

region_file[
  ,
  (required_region_columns) := lapply(.SD, as.integer),
  .SDcols = required_region_columns
]

if (nrow(region_file) == 0L) {
  stop("region_DMRcate_M12.csv contains no complete regions.")
}

if (anyDuplicated(region_file$region_id)) {
  stop("region_id must be unique.")
}

if (any(region_file$region_start > region_file$region_end)) {
  stop("At least one region has region_start > region_end.")
}

data.table::setorder(
  region_file,
  data_chunk_id,
  region_start,
  region_end,
  region_id
)

N_jobs_effective <- min(N_jobs_requested, nrow(region_file))

row_job_id <- ceiling(
  seq_len(nrow(region_file)) *
    N_jobs_effective /
    nrow(region_file)
)

region_row_list <- split(
  seq_len(nrow(region_file)),
  row_job_id
)

if (i_job > N_jobs_effective) {
  cat("No selected regions assigned to this task; exiting.\n")
  quit(save = "no", status = 0L)
}

i_region_rows <- region_row_list[[as.character(i_job)]]

if (is.null(i_region_rows) || length(i_region_rows) == 0L) {
  cat("No selected regions assigned to this task; exiting.\n")
  quit(save = "no", status = 0L)
}

regions_this_job <- region_file[i_region_rows]

if (!is.na(max_regions) && max_regions > 0L) {
  regions_this_job <- head(regions_this_job, max_regions)
}

cat(
  "Selected regions in file:", nrow(region_file), "\n",
  "Requested jobs:", N_jobs_requested, "\n",
  "Effective jobs:", N_jobs_effective, "\n",
  "Regions selected for this run:", nrow(regions_this_job), "\n"
)

# ---------------------------- family split file -------------------------------
if (!file.exists(train_split_path)) {
  stop("Training-split file does not exist: ", train_split_path)
}

train_FID_df <- data.table::fread(train_split_path)

if (!"splitID" %in% names(train_FID_df)) {
  stop("Training-split file must contain splitID.")
}

train_FID_cols <- grep(
  "^train_FID_[0-9]+$",
  names(train_FID_df),
  value = TRUE
)

if (length(train_FID_cols) == 0L) {
  stop("No train_FID_* columns found in training-split file.")
}

train_FID_cols <- train_FID_cols[
  order(
    as.integer(
      sub("^train_FID_", "", train_FID_cols)
    )
  )
]

train_FID_df[, splitID := as.integer(splitID)]

if (anyNA(train_FID_df$splitID)) {
  stop("At least one splitID is missing or nonnumeric.")
}

if (anyDuplicated(train_FID_df$splitID)) {
  stop("splitID values must be unique.")
}

split_ids <- train_FID_df[
  splitID >= split_start & splitID <= split_end,
  splitID
]

if (length(split_ids) == 0L) {
  stop("No splitID values fall within SPLIT_START:SPLIT_END.")
}

train_FID_list <- setNames(
  lapply(split_ids, function(sid) {
    z <- unlist(
      train_FID_df[
        splitID == sid,
        ..train_FID_cols
      ],
      use.names = FALSE
    )

    unique(
      as.character(
        z[
          !is.na(z) &
            nzchar(as.character(z))
        ]
      )
    )
  }),
  as.character(split_ids)
)

cat(
  "Splits selected for this run:", length(split_ids), "\n",
  "Training-FID columns per split:", length(train_FID_cols), "\n"
)

# ------------------------------- phenotype -----------------------------------
if (!file.exists(pheno_file_path)) {
  stop("Phenotype file does not exist: ", pheno_file_path)
}

pheno_env <- new.env(parent = emptyenv())
load(pheno_file_path, envir = pheno_env, verbose = TRUE)

if (!exists("pheno_file", envir = pheno_env, inherits = FALSE)) {
  stop("Phenotype RData does not contain pheno_file.")
}

pheno_file <- data.table::as.data.table(
  data.table::copy(
    get("pheno_file", envir = pheno_env, inherits = FALSE)
  )
)

rm(pheno_env)

design_variables <- c(
  "ID",
  "FID",
  "AA_only",
  "Sex",
  "AgeCalc",
  "Non.smoker",
  "EOSINOpc",
  "LYMPHOpc",
  "MONOpc",
  "NEUTROpc",
  "BMI",
  "sv1",
  "sv2",
  "sv3",
  "sv4",
  "sv5"
)

missing_pheno_columns <- setdiff(
  design_variables,
  names(pheno_file)
)

if (length(missing_pheno_columns) > 0L) {
  stop(
    "pheno_file is missing columns: ",
    paste(missing_pheno_columns, collapse = ", ")
  )
}

pheno_file <- pheno_file[, ..design_variables]
pheno_file[, ID := as.character(ID)]
pheno_file[, FID := as.character(FID)]
pheno_file[
  ,
  AA_only := suppressWarnings(
    as.integer(as.character(AA_only))
  )
]

if (any(!is.na(pheno_file$AA_only) &
        !pheno_file$AA_only %in% c(0L, 1L))) {
  stop("AA_only must contain only 0, 1, or NA.")
}

if (anyDuplicated(pheno_file$ID)) {
  stop("Duplicated IDs in pheno_file.")
}

# -------------------------------- helpers -------------------------------------
build_bsseq_from_region <- function(dt_region) {
  dt_region <- data.table::copy(dt_region)

  meth_cols <- grep("_meth$", names(dt_region), value = TRUE)

  if (length(meth_cols) == 0L) {
    stop("No *_meth columns found.")
  }

  sample_ids <- stringr::str_remove(meth_cols, "_meth$")
  cov_cols <- paste0(sample_ids, "_tot")
  missing_cov <- setdiff(cov_cols, names(dt_region))

  if (length(missing_cov) > 0L) {
    stop(
      "Missing matching coverage columns: ",
      paste(head(missing_cov, 10L), collapse = ", ")
    )
  }

  dt_region[
    ,
    (meth_cols) := lapply(.SD, as.numeric),
    .SDcols = meth_cols
  ]

  dt_region[
    ,
    (cov_cols) := lapply(.SD, as.numeric),
    .SDcols = cov_cols
  ]

  data.table::setorder(dt_region, start)

  meth_prop <- as.matrix(dt_region[, ..meth_cols])
  coverage <- as.matrix(dt_region[, ..cov_cols])

  storage.mode(meth_prop) <- "numeric"
  storage.mode(coverage) <- "numeric"

  meth_prop[is.na(meth_prop)] <- 0
  coverage[is.na(coverage)] <- 0
  meth_prop <- pmin(pmax(meth_prop, 0), 1)
  methylated_counts <- round(meth_prop * coverage)

  colnames(methylated_counts) <- sample_ids
  colnames(coverage) <- sample_ids

  bsseq::BSseq(
    M = methylated_counts,
    Cov = coverage,
    chr = dt_region$chr,
    pos = dt_region$start,
    sampleNames = sample_ids
  )
}

make_summary_row <- function(
  info,
  split_id,
  status,
  message = NA_character_,
  n_requested_fids = NA_integer_,
  n_observed_fids = NA_integer_,
  n_samples = NA_integer_,
  n_case = NA_integer_,
  n_ctrl = NA_integer_,
  n_cpg_before = NA_integer_,
  n_cpg_after = NA_integer_,
  n_seed_cpgs = NA_integer_,
  n_dmrs = NA_integer_,
  elapsed_seconds = NA_real_
) {
  data.table::data.table(
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
    n_cpgs_before_coverage = as.integer(n_cpg_before),
    n_cpgs_after_coverage = as.integer(n_cpg_after),
    n_seed_cpgs = as.integer(n_seed_cpgs),
    n_dmrs = as.integer(n_dmrs),
    fit_seconds = as.numeric(elapsed_seconds)
  )
}

append_summary <- function(row) {
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

# -------------------------------- outputs -------------------------------------
summary_out <- file.path(
  PATH_out,
  sprintf("dmrcate_training_summary_job_%04d.tsv", i_job)
)

dmr_out <- file.path(
  PATH_out,
  sprintf("dmrcate_training_dmrs_job_%04d.tsv", i_job)
)

if (overwrite) {
  unlink(c(summary_out, dmr_out), force = TRUE)
}

completed_keys <- character()

if (file.exists(summary_out)) {
  previous_summary <- data.table::fread(summary_out)

  required_previous <- c("region_id", "splitID", "status")

  if (all(required_previous %in% names(previous_summary))) {
    if (retry_errors) {
      previous_summary <- previous_summary[
        !status %in% c(
          "ANNOTATE_ERROR",
          "DMRCATE_ERROR",
          "EXTRACT_ERROR",
          "BSSEQ_ERROR",
          "DESIGN_ERROR"
        )
      ]
    }

    completed_keys <- unique(
      paste(
        previous_summary$region_id,
        previous_summary$splitID,
        sep = "::"
      )
    )
  }

  cat(
    "Existing completed region-split rows:",
    length(completed_keys),
    "\n"
  )
}

existing_dmr_keys <- character()

if (file.exists(dmr_out)) {
  previous_dmrs <- tryCatch(
    data.table::fread(dmr_out),
    error = function(e) NULL
  )

  if (
    !is.null(previous_dmrs) &&
      all(c("region_id", "splitID") %in% names(previous_dmrs))
  ) {
    existing_dmr_keys <- unique(
      paste(
        previous_dmrs$region_id,
        previous_dmrs$splitID,
        sep = "::"
      )
    )
  }
}

cat("Summary file:", summary_out, "\n")
cat("DMR file:", dmr_out, "\n")

# ------------------------------- main loop ------------------------------------
run_start_time <- Sys.time()
stop_requested <- FALSE

chunks <- split(
  regions_this_job,
  regions_this_job$data_chunk_id
)

for (chunk_name in names(chunks)) {
  if (stop_requested) break

  chunk_regions <- data.table::as.data.table(chunks[[chunk_name]])
  chunk_id_num <- as.integer(chunk_name)

  chunk_fp <- file.path(
    chunk_dir,
    sprintf("chunk_%04d.csv", chunk_id_num)
  )

  if (!file.exists(chunk_fp)) {
    for (row_i in seq_len(nrow(chunk_regions))) {
      info <- chunk_regions[row_i]

      for (split_id in split_ids) {
        key <- paste(info$region_id, split_id, sep = "::")
        if (key %in% completed_keys) next

        append_summary(
          make_summary_row(
            info,
            split_id,
            status = "MISSING_CHUNK_FILE",
            message = chunk_fp,
            n_requested_fids = length(
              train_FID_list[[as.character(split_id)]]
            )
          )
        )
        completed_keys <- c(completed_keys, key)
      }
    }
    next
  }

  header <- data.table::fread(
    chunk_fp,
    nrows = 0L,
    showProgress = FALSE
  )

  all_cols <- names(header)

  file_sample_ids <- unique(
    gsub(
      "_meth$|_tot$",
      "",
      setdiff(all_cols, c("chr", "start", "end"))
    )
  )

  phenotype_ids <- intersect(file_sample_ids, pheno_file$ID)

  if (length(phenotype_ids) == 0L) {
    cat("No overlapping phenotype samples in chunk", chunk_id_num, "\n")
    next
  }

  selected_cols <- intersect(
    c(
      "chr",
      "start",
      "end",
      paste0(phenotype_ids, "_meth"),
      paste0(phenotype_ids, "_tot")
    ),
    all_cols
  )

  cat(
    sprintf(
      "Reading chunk %d with %d selected columns (%d samples)\n",
      chunk_id_num,
      length(selected_cols),
      length(phenotype_ids)
    )
  )

  chunk_data <- data.table::fread(
    chunk_fp,
    select = selected_cols,
    showProgress = FALSE
  )

  chunk_data[, chr := as.character(chr)]
  chunk_data[, start := as.integer(start)]

  for (row_i in seq_len(nrow(chunk_regions))) {
    if (stop_requested) break

    info <- chunk_regions[row_i]
    info[, chr := as.character(chr)]

    cat(
      "\nRegion", info$region_id,
      "chunk", info$data_chunk_id,
      "\n"
    )

    region_all <- chunk_data[
      chr == info$chr &
        start >= info$region_start &
        start <= info$region_end
    ]

    n_cpg_before <- nrow(region_all)

    if (n_cpg_before < min_cpgs) {
      for (split_id in split_ids) {
        key <- paste(info$region_id, split_id, sep = "::")
        if (key %in% completed_keys) next

        append_summary(
          make_summary_row(
            info,
            split_id,
            status = "BELOW_MIN_CPGS_BEFORE",
            message = paste0(
              "CpGs before coverage filtering = ",
              n_cpg_before
            ),
            n_requested_fids = length(
              train_FID_list[[as.character(split_id)]]
            ),
            n_cpg_before = n_cpg_before
          )
        )
        completed_keys <- c(completed_keys, key)
      }
      next
    }

    for (split_id in split_ids) {
      elapsed_minutes <- as.numeric(
        difftime(
          Sys.time(),
          run_start_time,
          units = "mins"
        )
      )

      if (
        is.finite(max_runtime_minutes) &&
          max_runtime_minutes > 0 &&
          elapsed_minutes >= max_runtime_minutes
      ) {
        cat(
          "Reached MAX_RUNTIME_MINUTES =",
          max_runtime_minutes,
          "; stopping cleanly.\n"
        )
        stop_requested <- TRUE
        break
      }

      key <- paste(info$region_id, split_id, sep = "::")

      if (key %in% completed_keys) {
        next
      }

      fit_start <- proc.time()[["elapsed"]]
      train_fids <- train_FID_list[[as.character(split_id)]]
      n_requested_fids <- length(train_fids)

      if (n_requested_fids == 0L) {
        append_summary(
          make_summary_row(
            info,
            split_id,
            status = "NO_TRAIN_FIDS",
            n_requested_fids = 0L,
            n_cpg_before = n_cpg_before
          )
        )
        completed_keys <- c(completed_keys, key)
        next
      }

      ph_train <- pheno_file[
        FID %in% train_fids
      ]

      ph_train <- ph_train[
        complete.cases(ph_train[, ..design_variables])
      ]

      ph_train <- ph_train[
        ID %in% phenotype_ids
      ]

      if (nrow(ph_train) == 0L) {
        append_summary(
          make_summary_row(
            info,
            split_id,
            status = "NO_COMPLETE_TRAINING_PHENO",
            n_requested_fids = n_requested_fids,
            n_cpg_before = n_cpg_before
          )
        )
        completed_keys <- c(completed_keys, key)
        next
      }

      # Keep one phenotype row per sample and preserve phenotype order.
      if (anyDuplicated(ph_train$ID)) {
        stop("Duplicated training IDs after phenotype filtering.")
      }

      candidate_ids <- ph_train$ID
      meth_cols <- paste0(candidate_ids, "_meth")
      cov_cols <- paste0(candidate_ids, "_tot")
      paired <- meth_cols %in% names(region_all) &
        cov_cols %in% names(region_all)

      candidate_ids <- candidate_ids[paired]
      meth_cols <- paste0(candidate_ids, "_meth")
      cov_cols <- paste0(candidate_ids, "_tot")
      ph_train <- ph_train[match(candidate_ids, ID)]

      n_train_samples <- length(candidate_ids)
      n_train_fids_observed <- data.table::uniqueN(ph_train$FID)
      n_case <- sum(ph_train$AA_only == 1L)
      n_ctrl <- sum(ph_train$AA_only == 0L)

      if (n_train_samples < 2L) {
        append_summary(
          make_summary_row(
            info,
            split_id,
            status = "TOO_FEW_TRAINING_SAMPLES",
            n_requested_fids = n_requested_fids,
            n_observed_fids = n_train_fids_observed,
            n_samples = n_train_samples,
            n_case = n_case,
            n_ctrl = n_ctrl,
            n_cpg_before = n_cpg_before
          )
        )
        completed_keys <- c(completed_keys, key)
        next
      }

      if (n_case < 1L || n_ctrl < 1L) {
        append_summary(
          make_summary_row(
            info,
            split_id,
            status = "ONE_CLASS_ONLY",
            n_requested_fids = n_requested_fids,
            n_observed_fids = n_train_fids_observed,
            n_samples = n_train_samples,
            n_case = n_case,
            n_ctrl = n_ctrl,
            n_cpg_before = n_cpg_before
          )
        )
        completed_keys <- c(completed_keys, key)
        next
      }

      region_train <- region_all[
        ,
        c(
          intersect(c("chr", "start", "end"), names(region_all)),
          meth_cols,
          cov_cols
        ),
        with = FALSE
      ]

      region_train[
        ,
        (cov_cols) := lapply(.SD, as.numeric),
        .SDcols = cov_cols
      ]

      min_samples_cov <- max(
        2L,
        floor(coverage_fraction * length(cov_cols))
      )

      coverage_matrix <- as.matrix(
        region_train[, ..cov_cols]
      )
      storage.mode(coverage_matrix) <- "numeric"

      keep_row <- rowSums(
        coverage_matrix >= min_cov,
        na.rm = TRUE
      ) >= min_samples_cov

      region_train <- region_train[keep_row]
      n_cpg_after <- nrow(region_train)
      rm(coverage_matrix, keep_row)

      if (n_cpg_after < min_cpgs) {
        append_summary(
          make_summary_row(
            info,
            split_id,
            status = "BELOW_MIN_CPGS_AFTER",
            message = paste0(
              "CpGs after training-only coverage filtering = ",
              n_cpg_after
            ),
            n_requested_fids = n_requested_fids,
            n_observed_fids = n_train_fids_observed,
            n_samples = n_train_samples,
            n_case = n_case,
            n_ctrl = n_ctrl,
            n_cpg_before = n_cpg_before,
            n_cpg_after = n_cpg_after
          )
        )
        completed_keys <- c(completed_keys, key)
        next
      }

      # Build the BSseq object only after training-family restriction and
      # training-only coverage filtering.
      bs <- tryCatch(
        build_bsseq_from_region(region_train),
        error = function(e) e
      )

      if (inherits(bs, "error")) {
        append_summary(
          make_summary_row(
            info,
            split_id,
            status = "BSSEQ_ERROR",
            message = conditionMessage(bs),
            n_requested_fids = n_requested_fids,
            n_observed_fids = n_train_fids_observed,
            n_samples = n_train_samples,
            n_case = n_case,
            n_ctrl = n_ctrl,
            n_cpg_before = n_cpg_before,
            n_cpg_after = n_cpg_after,
            elapsed_seconds = proc.time()[["elapsed"]] - fit_start
          )
        )
        completed_keys <- c(completed_keys, key)
        next
      }

      sample_ids_bs <- colnames(
        SummarizedExperiment::assay(bs, "M")
      )

      ph_sub <- ph_train[match(sample_ids_bs, ID)]

      if (anyNA(ph_sub$ID)) {
        stop(
          "Internal phenotype/BSseq sample mismatch for region ",
          info$region_id,
          " split ",
          split_id
        )
      }

      design <- tryCatch(
        model.matrix(
          ~ AA_only +
            Sex +
            AgeCalc +
            Non.smoker +
            EOSINOpc +
            LYMPHOpc +
            MONOpc +
            NEUTROpc +
            BMI +
            sv1 +
            sv2 +
            sv3 +
            sv4 +
            sv5,
          data = ph_sub
        ),
        error = function(e) e
      )

      if (inherits(design, "error")) {
        append_summary(
          make_summary_row(
            info,
            split_id,
            status = "DESIGN_ERROR",
            message = conditionMessage(design),
            n_requested_fids = n_requested_fids,
            n_observed_fids = n_train_fids_observed,
            n_samples = n_train_samples,
            n_case = n_case,
            n_ctrl = n_ctrl,
            n_cpg_before = n_cpg_before,
            n_cpg_after = n_cpg_after,
            elapsed_seconds = proc.time()[["elapsed"]] - fit_start
          )
        )
        completed_keys <- c(completed_keys, key)
        rm(bs)
        next
      }

      if (nrow(design) != ncol(bs)) {
        stop(
          "Design/BSseq sample mismatch for region ",
          info$region_id,
          " split ",
          split_id
        )
      }

      colnames(design)[2L] <- "AA_only"
      rownames(design) <- colnames(bs)

      methdesign <- edgeR::modelMatrixMeth(design)
      colnames(methdesign) <- make.names(
        colnames(methdesign),
        unique = TRUE
      )

      aa_col <- make.names("AA_only")

      if (!aa_col %in% colnames(methdesign)) {
        append_summary(
          make_summary_row(
            info,
            split_id,
            status = "DESIGN_ERROR",
            message = paste0(
              "AA_only not found in methdesign; columns: ",
              paste(colnames(methdesign), collapse = ", ")
            ),
            n_requested_fids = n_requested_fids,
            n_observed_fids = n_train_fids_observed,
            n_samples = n_train_samples,
            n_case = n_case,
            n_ctrl = n_ctrl,
            n_cpg_before = n_cpg_before,
            n_cpg_after = n_cpg_after,
            elapsed_seconds = proc.time()[["elapsed"]] - fit_start
          )
        )
        completed_keys <- c(completed_keys, key)
        rm(bs, design, methdesign)
        next
      }

      contrast_matrix <- limma::makeContrasts(
        AA_vs_ctrl = AA_only,
        levels = methdesign
      )

      seq_annot <- tryCatch(
        DMRcate::sequencing.annotate(
          obj = bs,
          methdesign = methdesign,
          contrasts = TRUE,
          cont.matrix = contrast_matrix,
          coef = "AA_vs_ctrl",
          all.cov = FALSE,
          fdr = seed_fdr
        ),
        error = function(e) e
      )

      if (inherits(seq_annot, "error")) {
        append_summary(
          make_summary_row(
            info,
            split_id,
            status = "ANNOTATE_ERROR",
            message = conditionMessage(seq_annot),
            n_requested_fids = n_requested_fids,
            n_observed_fids = n_train_fids_observed,
            n_samples = n_train_samples,
            n_case = n_case,
            n_ctrl = n_ctrl,
            n_cpg_before = n_cpg_before,
            n_cpg_after = n_cpg_after,
            elapsed_seconds = proc.time()[["elapsed"]] - fit_start
          )
        )
        completed_keys <- c(completed_keys, key)
        rm(bs, design, methdesign, contrast_matrix)
        next
      }

      ranges_obj <- try(
        methods::slot(seq_annot, "ranges"),
        silent = TRUE
      )

      if (
        inherits(ranges_obj, "try-error") ||
          is.null(ranges_obj) ||
          length(ranges_obj) == 0L
      ) {
        append_summary(
          make_summary_row(
            info,
            split_id,
            status = "NO_ANNOTATION",
            n_requested_fids = n_requested_fids,
            n_observed_fids = n_train_fids_observed,
            n_samples = n_train_samples,
            n_case = n_case,
            n_ctrl = n_ctrl,
            n_cpg_before = n_cpg_before,
            n_cpg_after = n_cpg_after,
            n_seed_cpgs = 0L,
            n_dmrs = 0L,
            elapsed_seconds = proc.time()[["elapsed"]] - fit_start
          )
        )
        completed_keys <- c(completed_keys, key)
        rm(
          bs,
          design,
          methdesign,
          contrast_matrix,
          seq_annot,
          ranges_obj
        )
        next
      }

      metadata_df <- as.data.frame(
        S4Vectors::mcols(ranges_obj)
      )

      if ("is.sig" %in% names(metadata_df)) {
        n_seed_cpgs <- sum(metadata_df$is.sig, na.rm = TRUE)
      } else {
        p_values <- if ("ind.p" %in% names(metadata_df)) {
          metadata_df$ind.p
        } else if ("rawpval" %in% names(metadata_df)) {
          metadata_df$rawpval
        } else {
          numeric()
        }

        if (length(p_values) > 0L && any(is.finite(p_values))) {
          n_seed_cpgs <- sum(
            p.adjust(p_values, method = "BH") <= seed_fdr,
            na.rm = TRUE
          )
        } else {
          n_seed_cpgs <- 0L
        }
      }

      if (n_seed_cpgs < min_seed_cpgs) {
        append_summary(
          make_summary_row(
            info,
            split_id,
            status = "BELOW_MIN_SEEDS",
            message = paste0(
              "Significant seed CpGs = ",
              n_seed_cpgs,
              "; required = ",
              min_seed_cpgs
            ),
            n_requested_fids = n_requested_fids,
            n_observed_fids = n_train_fids_observed,
            n_samples = n_train_samples,
            n_case = n_case,
            n_ctrl = n_ctrl,
            n_cpg_before = n_cpg_before,
            n_cpg_after = n_cpg_after,
            n_seed_cpgs = n_seed_cpgs,
            n_dmrs = 0L,
            elapsed_seconds = proc.time()[["elapsed"]] - fit_start
          )
        )
        completed_keys <- c(completed_keys, key)
        rm(
          bs,
          design,
          methdesign,
          contrast_matrix,
          seq_annot,
          ranges_obj,
          metadata_df
        )
        next
      }

      dmrc <- try(
        DMRcate::dmrcate(
          seq_annot,
          lambda = lambda_value,
          C = C_value,
          min.cpgs = min_cpgs
        ),
        silent = TRUE
      )

      if (inherits(dmrc, "try-error")) {
        append_summary(
          make_summary_row(
            info,
            split_id,
            status = "DMRCATE_ERROR",
            message = as.character(dmrc),
            n_requested_fids = n_requested_fids,
            n_observed_fids = n_train_fids_observed,
            n_samples = n_train_samples,
            n_case = n_case,
            n_ctrl = n_ctrl,
            n_cpg_before = n_cpg_before,
            n_cpg_after = n_cpg_after,
            n_seed_cpgs = n_seed_cpgs,
            elapsed_seconds = proc.time()[["elapsed"]] - fit_start
          )
        )
        completed_keys <- c(completed_keys, key)
        rm(
          bs,
          design,
          methdesign,
          contrast_matrix,
          seq_annot,
          ranges_obj,
          metadata_df,
          dmrc
        )
        next
      }

      dmrs_gr <- tryCatch(
        safe_extract_gr(
          dmrc,
          cutoff_fdr = extract_fdr,
          min_cpgs = min_cpgs,
          min_abs_md = min_abs_md,
          min_width = min_width
        ),
        error = function(e) e
      )

      if (inherits(dmrs_gr, "error")) {
        append_summary(
          make_summary_row(
            info,
            split_id,
            status = "EXTRACT_ERROR",
            message = conditionMessage(dmrs_gr),
            n_requested_fids = n_requested_fids,
            n_observed_fids = n_train_fids_observed,
            n_samples = n_train_samples,
            n_case = n_case,
            n_ctrl = n_ctrl,
            n_cpg_before = n_cpg_before,
            n_cpg_after = n_cpg_after,
            n_seed_cpgs = n_seed_cpgs,
            elapsed_seconds = proc.time()[["elapsed"]] - fit_start
          )
        )
        completed_keys <- c(completed_keys, key)
        rm(
          bs,
          design,
          methdesign,
          contrast_matrix,
          seq_annot,
          ranges_obj,
          metadata_df,
          dmrc
        )
        next
      }

      n_dmrs <- length(dmrs_gr)
      fit_seconds <- proc.time()[["elapsed"]] - fit_start

      if (n_dmrs > 0L && !key %in% existing_dmr_keys) {
        dmr_table <- data.table::as.data.table(dmrs_gr)

        if (!"no.cpgs" %in% names(dmr_table)) {
          dmr_table[, no.cpgs := NA_integer_]
        }

        if (!"minfdr" %in% names(dmr_table)) {
          dmr_table[, minfdr := NA_real_]
        }

        dmr_table[
          ,
          `:=`(
            prep_job_id = i_job,
            region_id = as.integer(info$region_id),
            data_chunk_id = as.integer(info$data_chunk_id),
            parent_chr = as.character(info$chr),
            parent_start = as.integer(info$region_start),
            parent_end = as.integer(info$region_end),
            splitID = as.integer(split_id),
            n_train_FIDs_requested = n_requested_fids,
            n_train_FIDs_observed = n_train_fids_observed,
            n_train_samples = n_train_samples,
            n_train_cases = n_case,
            n_train_controls = n_ctrl,
            n_cpgs_after_coverage = n_cpg_after,
            n_seed_cpgs = n_seed_cpgs,
            score = -log10(
              pmax(minfdr, .Machine$double.xmin)
            )
          )
        ]

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
        make_summary_row(
          info,
          split_id,
          status = "FIT_OK",
          n_requested_fids = n_requested_fids,
          n_observed_fids = n_train_fids_observed,
          n_samples = n_train_samples,
          n_case = n_case,
          n_ctrl = n_ctrl,
          n_cpg_before = n_cpg_before,
          n_cpg_after = n_cpg_after,
          n_seed_cpgs = n_seed_cpgs,
          n_dmrs = n_dmrs,
          elapsed_seconds = fit_seconds
        )
      )

      completed_keys <- c(completed_keys, key)

      cat(
        "  completed split", split_id,
        "DMRs =", n_dmrs,
        "seeds =", n_seed_cpgs,
        "CpGs =", n_cpg_after,
        "samples =", n_train_samples,
        "FIDs =", n_train_fids_observed,
        "seconds =", round(fit_seconds, 2),
        "\n"
      )

      rm(
        bs,
        design,
        methdesign,
        contrast_matrix,
        seq_annot,
        ranges_obj,
        metadata_df,
        dmrc,
        dmrs_gr
      )
      invisible(gc(FALSE))
    }

    rm(region_all)
    invisible(gc(FALSE))
  }

  rm(
    chunk_data,
    header,
    all_cols,
    file_sample_ids,
    phenotype_ids,
    selected_cols,
    chunk_regions
  )
  invisible(gc(FALSE))
}

cat("Finished DMRcate family-level training analysis.\n")
cat("Summary:", summary_out, "\n")
cat("DMRs:", dmr_out, "\n")

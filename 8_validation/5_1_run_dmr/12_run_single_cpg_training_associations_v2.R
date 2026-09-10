#!/usr/bin/env Rscript

# Single-CpG split-specific association analysis for predictive benchmarking.
#
# One SLURM array task processes a deterministic contiguous block of rows from
# final_table. For every assigned CpG/model row, the script:
#   1. locates the CpG methylation vector by chr + start (end is never used);
#   2. uses data_chunk_id, or if that is NA, data_chunk_id_prev then
#      data_chunk_id_next;
#   3. repeats the row-specific association model separately in each of the
#      100 outer TRAINING-family splits;
#   4. writes one result row per final_table row x splitID.
#
# Models:
#   M1: two-stage logistic mixed model:
#       glmmLasso AA_only ~ methylation + covariates + (1|FID), lambda chosen
#       by minimum BIC, followed by unpenalized glmer refit of selected terms.
#
#   M2: standard LMM:
#       ASR methylation ~ AA_only + covariates + (1|FID).
#
#   M3: two-stage Gaussian mixed model:
#       glmmLasso ASR methylation ~ AA_only + covariates + (1|FID), lambda
#       chosen by minimum BIC, followed by unpenalized lmer refit of selected
#       terms.
#
# IMPORTANT:
# final_table used for the PRIMARY predictive benchmark must not itself be
# selected using AA outcomes from the full dataset. A full-data outcome-based
# prescreen would leak outer-test information. A restricted final_table is fine
# for smoke/pilot testing.

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(glmmLasso)
})

data.table::setDTthreads(1L)
options(mc.cores = 1L)

# ------------------------------- configuration -------------------------------

PATH_wk <- path.expand(Sys.getenv("PATH_WK", "~/scratch/UQAC/meth/"))

final_table_file <- path.expand(Sys.getenv(
  "FINAL_TABLE_FILE",
  file.path(
    PATH_wk,
    "results/15_revision/5_cv/6_single_cpg/final_table_for_cv.rds"
  )
))

train_split_file <- path.expand(Sys.getenv(
  "TRAIN_SPLIT_FILE",
  file.path(PATH_wk, "scr/11_mgcv/dat/100_family_training_splits.csv")
))

pheno_file_path <- path.expand(Sys.getenv(
  "PHENO_FILE",
  file.path(PATH_wk, "scr/11_mgcv/dat/18_pheno_BMI.RData")
))

chunk_dir <- path.expand(Sys.getenv(
  "METH_CHUNK_DIR",
  file.path(PATH_wk, "data/meth_split")
))

output_root <- path.expand(Sys.getenv(
  "SINGLE_CPG_ASSOC_OUTPUT",
  file.path(
    PATH_wk,
    "results/15_revision/5_cv/6_single_cpg/1_training_associations"
  )
))

n_jobs <- suppressWarnings(as.integer(Sys.getenv("N_JOBS", "990")))
min_complete_n <- suppressWarnings(as.integer(Sys.getenv("MIN_COMPLETE_N", "31")))
cache_chunks <- suppressWarnings(as.integer(Sys.getenv("CACHE_CHUNKS", "3")))
max_rows <- suppressWarnings(as.integer(Sys.getenv("MAX_ROWS", "0")))
overwrite <- Sys.getenv("OVERWRITE", "0") %in%
  c("1", "TRUE", "true", "T", "yes", "YES")

# Optional smoke-test controls:
#   SPLIT_IDS=1,2,3
#   MAX_SPLITS=2
split_ids_env <- trimws(Sys.getenv("SPLIT_IDS", ""))
max_splits <- suppressWarnings(as.integer(Sys.getenv("MAX_SPLITS", "0")))

lambda_grid <- c(5, 10, 20, 40, 50, 60, 80, 100)

if (is.na(n_jobs) || n_jobs < 1L) stop("N_JOBS must be >= 1.")
if (is.na(min_complete_n) || min_complete_n < 10L) {
  stop("MIN_COMPLETE_N must be >= 10.")
}
if (is.na(cache_chunks) || cache_chunks < 1L) cache_chunks <- 1L
if (is.na(max_rows) || max_rows < 0L) max_rows <- 0L
if (is.na(max_splits) || max_splits < 0L) max_splits <- 0L

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
result_dir <- file.path(output_root, "job_results")
qc_dir <- file.path(output_root, "job_qc")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------- arguments ----------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Pass SLURM_ARRAY_TASK_ID / job_id as the first argument.")
}

job_id <- suppressWarnings(as.integer(args[1L]))
if (is.na(job_id) || job_id < 1L || job_id > n_jobs) {
  stop("job_id must be between 1 and ", n_jobs, ".")
}

run_start <- Sys.time()

cat("Single-CpG training-split associations\n")
cat("job_id:", job_id, "of", n_jobs, "\n")
cat("final_table:", final_table_file, "\n")

result_file <- file.path(
  result_dir,
  sprintf("single_cpg_training_assoc_job_%04d.tsv", job_id)
)

qc_file <- file.path(
  qc_dir,
  sprintf("single_cpg_training_assoc_job_%04d_qc.tsv", job_id)
)

if (!overwrite && file.exists(result_file) && file.exists(qc_file)) {
  qc_old <- tryCatch(fread(qc_file), error = function(e) NULL)
  if (!is.null(qc_old) &&
      "job_status" %in% names(qc_old) &&
      any(qc_old$job_status == "COMPLETED")) {
    cat("Completed outputs already exist; exiting.\n")
    quit(save = "no", status = 0L)
  }
}

if (!overwrite && file.exists(result_file) && !file.exists(qc_file)) {
  stop(
    "A result file exists without completed QC. It may be partial. ",
    "Rerun with OVERWRITE=1: ", result_file
  )
}

if (overwrite) {
  unlink(c(result_file, qc_file), force = TRUE)
}

# -------------------------------- utilities ----------------------------------

normalize_id <- function(x) {
  x <- as.character(x)
  x <- sub("^X", "", x)
  x <- sub("_meth$", "", x)
  gsub("\\.", "-", x)
}

normalize_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  sub("\\.0$", "", x)
}

finite_numeric <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  out[!is.finite(out)] <- NA_real_
  out
}

first_finite <- function(x) {
  x <- finite_numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else x[1L]
}

safe_fid_sd <- function(fit) {
  out <- tryCatch({
    vc <- lme4::VarCorr(fit)
    if (!"FID" %in% names(vc)) return(NA_real_)
    as.numeric(attr(vc$FID, "stddev"))[1L]
  }, error = function(e) NA_real_)
  out
}

safe_glmmLasso_sd <- function(fit) {
  out <- tryCatch(as.numeric(fit$StdDev)[1L], error = function(e) NA_real_)
  if (!is.finite(out)) NA_real_ else out
}

find_col <- function(tab, patterns) {
  if (is.null(tab) || is.null(colnames(tab))) return(NA_character_)
  nms <- colnames(tab)
  for (pat in patterns) {
    hit <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0L) return(hit[1L])
  }
  NA_character_
}

extract_term_from_table <- function(tab, term) {
  out <- list(
    estimate = NA_real_,
    std_error = NA_real_,
    statistic = NA_real_,
    df = NA_real_,
    p_value = NA_real_
  )

  if (is.null(tab) || is.null(rownames(tab)) || !term %in% rownames(tab)) {
    return(out)
  }

  est_col <- find_col(tab, c("^Estimate$", "^Coef", "estimate"))
  se_col <- find_col(tab, c("Std\\. Error", "Std\\.Err", "Std\\.Error", "^SE$"))
  stat_col <- find_col(tab, c("z value", "t value", "Wald", "stat"))
  df_col <- find_col(tab, c("^df$", "d\\.f"))
  p_col <- find_col(tab, c("Pr\\(", "p\\.value", "p-value", "^p$"))

  if (!is.na(est_col)) out$estimate <- first_finite(tab[term, est_col])
  if (!is.na(se_col)) out$std_error <- first_finite(tab[term, se_col])
  if (!is.na(stat_col)) out$statistic <- first_finite(tab[term, stat_col])
  if (!is.na(df_col)) out$df <- first_finite(tab[term, df_col])
  if (!is.na(p_col)) out$p_value <- first_finite(tab[term, p_col])

  out
}

extract_glmmLasso_term <- function(fit, term) {
  co <- tryCatch(stats::coef(fit), error = function(e) NULL)
  estimate <- NA_real_
  if (!is.null(co) && !is.null(names(co)) && term %in% names(co)) {
    estimate <- first_finite(co[term])
  }

  tab <- tryCatch(
    suppressWarnings(stats::coef(summary(fit))),
    error = function(e) NULL
  )
  ext <- extract_term_from_table(tab, term)

  if (!is.finite(ext$estimate) && is.finite(estimate)) {
    ext$estimate <- estimate
  }
  ext
}

make_base_result <- function(
  info,
  split_id,
  resolved_chunk_id = NA_integer_,
  chunk_source = NA_character_,
  cpg_row_in_chunk = NA_integer_
) {
  data.table(
    job_id = job_id,
    final_table_row_id = as.integer(info$final_table_row_id),
    splitID = as.integer(split_id),
    CpG = as.character(info$CpG),
    requested_model = as.character(info$model),
    cpg_chr = as.character(info$cpg_chr),
    cpg_start = as.integer(info$cpg_start),
    data_chunk_id = as.integer(info$data_chunk_id),
    data_chunk_id_prev = as.integer(info$data_chunk_id_prev),
    data_chunk_id_next = as.integer(info$data_chunk_id_next),
    resolved_chunk_id = as.integer(resolved_chunk_id),
    chunk_source = as.character(chunk_source),
    cpg_row_in_chunk = as.integer(cpg_row_in_chunk),
    original_full_data_pval = finite_numeric(info$pval)[1L],
    original_full_data_coef = finite_numeric(info$coef)[1L],
    genes_entrez = as.character(info$genes_entrez),
    genes_symbol = as.character(info$genes_symbol),
    association_term = NA_character_,
    coefficient = NA_real_,
    std_error = NA_real_,
    statistic = NA_real_,
    df = NA_real_,
    p_value = NA_real_,
    FID_sd = NA_real_,
    optimal_lambda = NA_real_,
    optimal_lambda_BIC = NA_real_,
    optimal_lambda_AIC = NA_real_,
    stage1_coefficient = NA_real_,
    stage1_p_value = NA_real_,
    selected_stage1 = NA,
    final_fit_stage = NA_character_,
    n_train_total = NA_integer_,
    n_train_complete = NA_integer_,
    n_train_FIDs = NA_integer_,
    n_cases = NA_integer_,
    n_controls = NA_integer_,
    n_meth_observed = NA_integer_,
    meth_sd = NA_real_,
    fit_status = NA_character_,
    error_message = NA_character_
  )
}

# ----------------------------- load final table -------------------------------

load_final_table <- function(path) {
  if (!file.exists(path)) stop("final_table file does not exist: ", path)

  low <- tolower(path)

  if (grepl("\\.rds$", low)) {
    x <- readRDS(path)
  } else if (grepl("\\.(tsv|txt|csv)(\\.gz)?$", low)) {
    x <- fread(path)
  } else if (grepl("\\.(rdata|rda)$", low)) {
    e <- new.env(parent = emptyenv())
    load(path, envir = e)
    if (!exists("final_table", envir = e, inherits = FALSE)) {
      stop("RData/RDA must contain an object named final_table.")
    }
    x <- get("final_table", envir = e, inherits = FALSE)
  } else {
    stop("Unsupported FINAL_TABLE_FILE extension: ", path)
  }

  as.data.table(copy(x))
}

final_table <- load_final_table(final_table_file)

required_final_cols <- c(
  "CpG", "pval", "coef", "genes_entrez", "genes_symbol", "model",
  "cpg_chr", "cpg_start",
  "data_chunk_id", "data_chunk_id_prev", "data_chunk_id_next"
)

missing_final <- setdiff(required_final_cols, names(final_table))
if (length(missing_final) > 0L) {
  stop(
    "final_table is missing required columns: ",
    paste(missing_final, collapse = ", ")
  )
}

final_table[, final_table_row_id := .I]
final_table[, cpg_chr := normalize_chr(cpg_chr)]
final_table[, cpg_start := suppressWarnings(as.integer(cpg_start))]
final_table[, model := toupper(trimws(as.character(model)))]

for (z in c("data_chunk_id", "data_chunk_id_prev", "data_chunk_id_next")) {
  set(final_table, j = z, value = suppressWarnings(as.integer(final_table[[z]])))
}

if (any(!final_table$model %in% c("M1", "M2", "M3"))) {
  bad <- unique(final_table[!model %in% c("M1", "M2", "M3"), model])
  stop("Unsupported model labels in final_table: ", paste(bad, collapse = ", "))
}

if (any(is.na(final_table$cpg_chr) | !nzchar(final_table$cpg_chr) |
        is.na(final_table$cpg_start))) {
  stop("final_table contains missing cpg_chr or cpg_start.")
}

n_total_rows <- nrow(final_table)

# Deterministic contiguous partition. Every row is assigned exactly once.
row_start <- floor((job_id - 1L) * n_total_rows / n_jobs) + 1L
row_end <- floor(job_id * n_total_rows / n_jobs)

if (row_end < row_start || n_total_rows == 0L) {
  job_rows <- final_table[0]
} else {
  job_rows <- final_table[row_start:row_end]
}

if (max_rows > 0L && nrow(job_rows) > max_rows) {
  job_rows <- head(job_rows, max_rows)
}

cat(
  "Total final_table rows:", n_total_rows,
  "; assigned rows:", nrow(job_rows),
  "; original range:", row_start, "-", row_end, "\n"
)

# ----------------------- load split and phenotype data ------------------------

if (!file.exists(train_split_file)) {
  stop("Training-split file does not exist: ", train_split_file)
}
train_split_dt <- fread(train_split_file)

train_cols <- grep("^train_FID_[0-9]+$", names(train_split_dt), value = TRUE)
if (!"splitID" %in% names(train_split_dt) || length(train_cols) == 0L) {
  stop("Training-split file lacks splitID or train_FID_* columns.")
}

split_ids <- sort(unique(as.integer(train_split_dt$splitID)))
split_ids <- split_ids[is.finite(split_ids)]

if (nzchar(split_ids_env)) {
  requested <- suppressWarnings(
    as.integer(strsplit(split_ids_env, ",", fixed = TRUE)[[1L]])
  )
  requested <- requested[is.finite(requested)]
  split_ids <- intersect(split_ids, requested)
}

if (max_splits > 0L && length(split_ids) > max_splits) {
  split_ids <- head(split_ids, max_splits)
}

if (length(split_ids) == 0L) stop("No splitIDs selected.")

train_fid_list <- lapply(split_ids, function(sid) {
  rr <- train_split_dt[splitID == sid]
  if (nrow(rr) != 1L) stop("Expected one row for splitID ", sid)
  fids <- trimws(as.character(unlist(rr[, ..train_cols], use.names = FALSE)))
  unique(fids[!is.na(fids) & nzchar(fids) & fids != "NA"])
})
names(train_fid_list) <- as.character(split_ids)

if (!file.exists(pheno_file_path)) {
  stop("Phenotype file does not exist: ", pheno_file_path)
}

pe <- new.env(parent = emptyenv())
load(pheno_file_path, envir = pe)
if (!exists("pheno_file", envir = pe, inherits = FALSE)) {
  stop("Phenotype RData does not contain pheno_file.")
}
pheno_dt <- as.data.table(copy(get("pheno_file", envir = pe, inherits = FALSE)))
rm(pe)

covariates <- c(
  "AgeCalc", "Sex", "Non.smoker",
  "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc",
  "sv1", "sv2", "sv3", "sv4", "sv5", "BMI"
)

required_pheno <- c("ID", "FID", "AA_only", covariates)
missing_pheno <- setdiff(required_pheno, names(pheno_dt))
if (length(missing_pheno) > 0L) {
  stop("pheno_file is missing: ", paste(missing_pheno, collapse = ", "))
}

pheno_dt <- pheno_dt[, ..required_pheno]
pheno_dt[, ID := normalize_id(ID)]
pheno_dt[, FID := as.character(FID)]
pheno_dt[, AA_only := suppressWarnings(as.integer(as.character(AA_only)))]

pheno_dt <- pheno_dt[
  !is.na(ID) & nzchar(ID) &
    !is.na(FID) & nzchar(FID) &
    AA_only %in% c(0L, 1L)
]

if (anyDuplicated(pheno_dt$ID)) {
  stop("Duplicated normalized IDs in phenotype data.")
}

pheno_ids <- pheno_dt$ID

# --------------------------- methylation chunk cache -------------------------

chunk_cache <- new.env(parent = emptyenv())
cache_order <- integer()

get_chunk <- function(chunk_id) {
  key <- as.character(as.integer(chunk_id))

  if (exists(key, envir = chunk_cache, inherits = FALSE)) {
    # refresh LRU order
    cache_order <<- c(cache_order[cache_order != as.integer(chunk_id)],
                      as.integer(chunk_id))
    return(get(key, envir = chunk_cache, inherits = FALSE))
  }

  chunk_file <- file.path(chunk_dir, sprintf("chunk_%04d.csv", as.integer(chunk_id)))
  if (!file.exists(chunk_file)) {
    return(structure(
      list(message = paste0("Missing chunk file: ", chunk_file)),
      class = "chunk_error"
    ))
  }

  header <- tryCatch(
    fread(chunk_file, nrows = 0L, showProgress = FALSE),
    error = function(e) e
  )
  if (inherits(header, "error")) {
    return(structure(
      list(message = conditionMessage(header)),
      class = "chunk_error"
    ))
  }

  hnames <- names(header)
  if (!all(c("chr", "start") %in% hnames)) {
    return(structure(
      list(message = "Chunk lacks chr/start columns."),
      class = "chunk_error"
    ))
  }

  meth_cols <- grep("_meth$", hnames, value = TRUE)
  meth_ids <- normalize_id(meth_cols)
  keep_meth <- meth_ids %chin% pheno_ids

  meth_cols <- meth_cols[keep_meth]
  meth_ids <- meth_ids[keep_meth]

  if (length(meth_cols) == 0L) {
    return(structure(
      list(message = "No methylation columns match phenotype IDs."),
      class = "chunk_error"
    ))
  }

  if (anyDuplicated(meth_ids)) {
    return(structure(
      list(message = "Duplicated normalized methylation IDs in chunk."),
      class = "chunk_error"
    ))
  }

  dt <- tryCatch(
    fread(
      chunk_file,
      select = c("chr", "start", meth_cols),
      showProgress = FALSE
    ),
    error = function(e) e
  )

  if (inherits(dt, "error")) {
    return(structure(
      list(message = conditionMessage(dt)),
      class = "chunk_error"
    ))
  }

  dt[, chr := normalize_chr(chr)]
  dt[, start := suppressWarnings(as.integer(start))]

  obj <- list(
    dt = dt,
    meth_cols = meth_cols,
    meth_ids = meth_ids,
    chunk_file = chunk_file
  )

  assign(key, obj, envir = chunk_cache)
  cache_order <<- c(cache_order[cache_order != as.integer(chunk_id)],
                    as.integer(chunk_id))

  while (length(cache_order) > cache_chunks) {
    old <- cache_order[1L]
    cache_order <<- cache_order[-1L]
    rm(list = as.character(old), envir = chunk_cache)
    invisible(gc(FALSE))
  }

  obj
}

resolve_cpg <- function(info) {
  primary <- suppressWarnings(as.integer(info$data_chunk_id))
  prev <- suppressWarnings(as.integer(info$data_chunk_id_prev))
  next_id <- suppressWarnings(as.integer(info$data_chunk_id_next))

  if (is.finite(primary)) {
    candidate_ids <- primary
    candidate_sources <- "data_chunk_id"
  } else {
    candidate_ids <- c(prev, next_id)
    candidate_sources <- c("data_chunk_id_prev", "data_chunk_id_next")
    keep <- is.finite(candidate_ids)
    candidate_ids <- candidate_ids[keep]
    candidate_sources <- candidate_sources[keep]

    # Avoid reading the same chunk twice when prev == next.
    if (length(candidate_ids) > 0L) {
      keep_unique <- !duplicated(candidate_ids)
      candidate_ids <- candidate_ids[keep_unique]
      candidate_sources <- candidate_sources[keep_unique]
    }
  }

  if (length(candidate_ids) == 0L) {
    return(list(
      status = "NO_CHUNK_ID_AVAILABLE",
      message = "data_chunk_id, prev, and next are all NA.",
      meth = NULL,
      resolved_chunk_id = NA_integer_,
      chunk_source = NA_character_,
      cpg_row_in_chunk = NA_integer_
    ))
  }

  failure_messages <- character()

  for (k in seq_along(candidate_ids)) {
    cid <- as.integer(candidate_ids[k])
    source_name <- candidate_sources[k]

    obj <- get_chunk(cid)
    if (inherits(obj, "chunk_error")) {
      failure_messages <- c(
        failure_messages,
        paste0(source_name, "=", cid, ": ", obj$message)
      )
      next
    }

    hit <- which(
      obj$dt$chr == as.character(info$cpg_chr) &
        obj$dt$start == as.integer(info$cpg_start)
    )

    if (length(hit) == 0L) {
      failure_messages <- c(
        failure_messages,
        paste0(
          source_name, "=", cid,
          ": CpG not found at ",
          info$cpg_chr, ":", info$cpg_start
        )
      )
      next
    }

    # If duplicated coordinate rows exist, use the first deterministically and
    # record the condition in the status.
    rr <- hit[1L]

    # IMPORTANT: first copy the dynamic column vector into a simple local
    # variable. This mirrors the working regional scripts:
    #     source_cols <- meth_map$source_column
    #     meth_source <- as.matrix(region_dt[, ..source_cols])
    #
    # Using ..obj$meth_cols directly inside data.table j is not reliable and
    # can yield non-numeric output / all-NA coercion.
    source_cols <- obj$meth_cols
    vals_mat <- as.matrix(obj$dt[rr, ..source_cols])
    storage.mode(vals_mat) <- "numeric"

    vals <- as.numeric(vals_mat[1L, ])
    vals[!is.finite(vals) | vals < 0 | vals > 1] <- NA_real_
    names(vals) <- obj$meth_ids

    return(list(
      status = if (length(hit) > 1L) "FOUND_DUPLICATE_COORD_FIRST_USED" else "FOUND",
      message = if (length(hit) > 1L) {
        paste0("Coordinate occurred ", length(hit), " times; first row used.")
      } else {
        NA_character_
      },
      meth = vals,
      resolved_chunk_id = cid,
      chunk_source = source_name,
      cpg_row_in_chunk = rr
    ))
  }

  list(
    status = "CPG_NOT_FOUND_IN_CANDIDATE_CHUNKS",
    message = paste(failure_messages, collapse = " | "),
    meth = NULL,
    resolved_chunk_id = NA_integer_,
    chunk_source = NA_character_,
    cpg_row_in_chunk = NA_integer_
  )
}

# ------------------------------ model helpers --------------------------------

fixed_cov_formula <- paste(covariates, collapse = " + ")

fit_model1 <- function(dat) {
  target <- "meth_response"

  bic_values <- rep(NA_real_, length(lambda_grid))
  aic_values <- rep(NA_real_, length(lambda_grid))

  for (i in seq_along(lambda_grid)) {
    tmp <- tryCatch(
      suppressWarnings(
        glmmLasso::glmmLasso(
          fix = as.formula(
            paste0("AA_only ~ meth_response + ", fixed_cov_formula)
          ),
          rnd = list(FID = ~1),
          data = dat,
          family = binomial(link = "logit"),
          lambda = lambda_grid[i]
        )
      ),
      error = function(e) NULL
    )

    if (!is.null(tmp)) {
      bic_values[i] <- first_finite(tmp$bic)
      aic_values[i] <- first_finite(tmp$aic)
    }
  }

  ok <- which(is.finite(bic_values))
  if (length(ok) == 0L) {
    stop("All Model 1 lambda fits failed or returned non-finite BIC.")
  }

  best_i <- ok[which.min(bic_values[ok])]
  optimal_lambda <- lambda_grid[best_i]

  mod <- tryCatch(
    suppressWarnings(
      glmmLasso::glmmLasso(
        fix = as.formula(
          paste0("AA_only ~ meth_response + ", fixed_cov_formula)
        ),
        rnd = list(FID = ~1),
        data = dat,
        family = binomial(link = "logit"),
        lambda = optimal_lambda,
        final.re = TRUE
      )
    ),
    error = function(e) e
  )

  if (inherits(mod, "error")) {
    stop("Model 1 optimal glmmLasso failed: ", conditionMessage(mod))
  }

  stage1 <- extract_glmmLasso_term(mod, target)
  co <- tryCatch(stats::coef(mod), error = function(e) numeric())
  selected <- if (!is.null(names(co))) {
    setdiff(names(co)[is.finite(co) & co != 0], "(Intercept)")
  } else {
    character()
  }
  selected_target <- target %in% selected

  final <- stage1
  fid_sd <- safe_glmmLasso_sd(mod)
  final_stage <- "glmmLasso_stage1"

  if (length(selected) > 0L) {
    f2 <- as.formula(
      paste0("AA_only ~ ", paste(selected, collapse = " + "), " + (1|FID)")
    )

    refit <- tryCatch(
      suppressWarnings(
        lme4::glmer(
          formula = f2,
          data = dat,
          family = binomial(link = "logit"),
          control = lme4::glmerControl(
            optimizer = "bobyqa",
            optCtrl = list(maxfun = 100000)
          )
        )
      ),
      error = function(e) e
    )

    if (!inherits(refit, "error")) {
      tab <- summary(refit)$coefficients
      if (target %in% rownames(tab)) {
        final <- extract_term_from_table(tab, target)
        final_stage <- "glmer_stage2"
      } else {
        final_stage <- "glmer_stage2_target_not_selected"
      }
      fid_sd <- safe_fid_sd(refit)
    } else {
      final_stage <- paste0(
        "glmmLasso_stage1_refit_failed: ",
        conditionMessage(refit)
      )
    }
  } else {
    final_stage <- "glmmLasso_stage1_no_terms_selected"
  }

  list(
    association_term = target,
    coefficient = final$estimate,
    std_error = final$std_error,
    statistic = final$statistic,
    df = final$df,
    p_value = final$p_value,
    FID_sd = fid_sd,
    optimal_lambda = optimal_lambda,
    optimal_lambda_BIC = bic_values[best_i],
    optimal_lambda_AIC = aic_values[best_i],
    stage1_coefficient = stage1$estimate,
    stage1_p_value = stage1$p_value,
    selected_stage1 = selected_target,
    final_fit_stage = final_stage
  )
}

fit_model2 <- function(dat) {
  target <- "AA_only"

  fit <- tryCatch(
    suppressWarnings(
      lmerTest::lmer(
        formula = as.formula(
          paste0(
            "meth_response ~ AA_only + ",
            fixed_cov_formula,
            " + (1|FID)"
          )
        ),
        data = dat,
        REML = FALSE
      )
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    stop("Model 2 lmer failed: ", conditionMessage(fit))
  }

  tab <- summary(fit)$coefficients
  ext <- extract_term_from_table(tab, target)

  list(
    association_term = target,
    coefficient = ext$estimate,
    std_error = ext$std_error,
    statistic = ext$statistic,
    df = ext$df,
    p_value = ext$p_value,
    FID_sd = safe_fid_sd(fit),
    optimal_lambda = NA_real_,
    optimal_lambda_BIC = NA_real_,
    optimal_lambda_AIC = NA_real_,
    stage1_coefficient = NA_real_,
    stage1_p_value = NA_real_,
    selected_stage1 = NA,
    final_fit_stage = "lmer_direct"
  )
}

fit_model3 <- function(dat) {
  target <- "AA_only"

  bic_values <- rep(NA_real_, length(lambda_grid))
  aic_values <- rep(NA_real_, length(lambda_grid))

  for (i in seq_along(lambda_grid)) {
    tmp <- tryCatch(
      suppressWarnings(
        glmmLasso::glmmLasso(
          fix = as.formula(
            paste0("meth_response ~ AA_only + ", fixed_cov_formula)
          ),
          rnd = list(FID = ~1),
          data = dat,
          family = gaussian(link = "identity"),
          lambda = lambda_grid[i]
        )
      ),
      error = function(e) NULL
    )

    if (!is.null(tmp)) {
      bic_values[i] <- first_finite(tmp$bic)
      aic_values[i] <- first_finite(tmp$aic)
    }
  }

  ok <- which(is.finite(bic_values))
  if (length(ok) == 0L) {
    stop("All Model 3 lambda fits failed or returned non-finite BIC.")
  }

  best_i <- ok[which.min(bic_values[ok])]
  optimal_lambda <- lambda_grid[best_i]

  mod <- tryCatch(
    suppressWarnings(
      glmmLasso::glmmLasso(
        fix = as.formula(
          paste0("meth_response ~ AA_only + ", fixed_cov_formula)
        ),
        rnd = list(FID = ~1),
        data = dat,
        family = gaussian(link = "identity"),
        lambda = optimal_lambda,
        final.re = TRUE
      )
    ),
    error = function(e) e
  )

  if (inherits(mod, "error")) {
    stop("Model 3 optimal glmmLasso failed: ", conditionMessage(mod))
  }

  stage1 <- extract_glmmLasso_term(mod, target)
  co <- tryCatch(stats::coef(mod), error = function(e) numeric())
  selected <- if (!is.null(names(co))) {
    setdiff(names(co)[is.finite(co) & co != 0], "(Intercept)")
  } else {
    character()
  }
  selected_target <- target %in% selected

  final <- stage1
  fid_sd <- safe_glmmLasso_sd(mod)
  final_stage <- "glmmLasso_stage1"

  if (length(selected) > 0L) {
    f2 <- as.formula(
      paste0(
        "meth_response ~ ",
        paste(selected, collapse = " + "),
        " + (1|FID)"
      )
    )

    refit <- tryCatch(
      suppressWarnings(
        lmerTest::lmer(
          formula = f2,
          data = dat,
          REML = FALSE
        )
      ),
      error = function(e) e
    )

    if (!inherits(refit, "error")) {
      tab <- summary(refit)$coefficients
      if (target %in% rownames(tab)) {
        final <- extract_term_from_table(tab, target)
        final_stage <- "lmer_stage2"
      } else {
        final_stage <- "lmer_stage2_target_not_selected"
      }
      fid_sd <- safe_fid_sd(refit)
    } else {
      final_stage <- paste0(
        "glmmLasso_stage1_refit_failed: ",
        conditionMessage(refit)
      )
    }
  } else {
    final_stage <- "glmmLasso_stage1_no_terms_selected"
  }

  list(
    association_term = target,
    coefficient = final$estimate,
    std_error = final$std_error,
    statistic = final$statistic,
    df = final$df,
    p_value = final$p_value,
    FID_sd = fid_sd,
    optimal_lambda = optimal_lambda,
    optimal_lambda_BIC = bic_values[best_i],
    optimal_lambda_AIC = aic_values[best_i],
    stage1_coefficient = stage1$estimate,
    stage1_p_value = stage1$p_value,
    selected_stage1 = selected_target,
    final_fit_stage = final_stage
  )
}

fit_requested_model <- function(model_label, dat) {
  switch(
    model_label,
    M1 = fit_model1(dat),
    M2 = fit_model2(dat),
    M3 = fit_model3(dat),
    stop("Unsupported model label: ", model_label)
  )
}

# ------------------------------ analysis loop --------------------------------

n_output <- 0L
n_fit_ok <- 0L
n_fit_error <- 0L
n_cpg_not_found <- 0L
n_rows_processed <- 0L

for (j in seq_len(nrow(job_rows))) {
  info <- job_rows[j]

  cat(
    sprintf(
      "[job %04d] row %d/%d; final_table_row_id=%d; %s; %s\n",
      job_id, j, nrow(job_rows),
      info$final_table_row_id,
      info$CpG,
      info$model
    )
  )

  cpg_obj <- resolve_cpg(info)
  n_rows_processed <- n_rows_processed + 1L

  if (is.null(cpg_obj$meth)) {
    n_cpg_not_found <- n_cpg_not_found + 1L

    failed_rows <- rbindlist(lapply(split_ids, function(sid) {
      out <- make_base_result(
        info = info,
        split_id = sid,
        resolved_chunk_id = cpg_obj$resolved_chunk_id,
        chunk_source = cpg_obj$chunk_source,
        cpg_row_in_chunk = cpg_obj$cpg_row_in_chunk
      )
      out[, `:=`(
        fit_status = cpg_obj$status,
        error_message = cpg_obj$message
      )]
      out
    }))

    fwrite(
      failed_rows,
      result_file,
      sep = "\t",
      quote = FALSE,
      na = "NA",
      append = file.exists(result_file),
      col.names = !file.exists(result_file)
    )

    n_output <- n_output + nrow(failed_rows)
    next
  }

  meth_vec <- cpg_obj$meth

  # Explicit ID matching is easier to audit than character subscripting.
  meth_match_idx <- match(pheno_dt$ID, names(meth_vec))
  meth_by_pheno <- meth_vec[meth_match_idx]
  meth_by_pheno[!is.finite(meth_by_pheno)] <- NA_real_

  cat(
    "  methylation matched phenotype IDs:",
    sum(!is.na(meth_match_idx)),
    "; finite beta values:",
    sum(is.finite(meth_by_pheno)),
    "\n"
  )

  if (any(is.finite(meth_by_pheno))) {
    cat(
      "  beta range:",
      paste(signif(range(meth_by_pheno, na.rm = TRUE), 5), collapse = " to "),
      "\n"
    )
  }

  meth_asr <- asin(sqrt(meth_by_pheno))

  split_results <- vector("list", length(split_ids))

  for (ss in seq_along(split_ids)) {
    sid <- split_ids[ss]
    train_fids <- train_fid_list[[as.character(sid)]]

    idx_train <- pheno_dt$FID %chin% train_fids
    dat_dt <- copy(pheno_dt[idx_train])
    dat_dt[, meth_response := meth_asr[idx_train]]

    required_model_cols <- c(
      "AA_only", "FID", "meth_response", covariates
    )

    complete_idx <- complete.cases(dat_dt[, ..required_model_cols])
    dat <- as.data.frame(dat_dt[complete_idx, ..required_model_cols])
    dat$FID <- factor(dat$FID)

    out <- make_base_result(
      info = info,
      split_id = sid,
      resolved_chunk_id = cpg_obj$resolved_chunk_id,
      chunk_source = cpg_obj$chunk_source,
      cpg_row_in_chunk = cpg_obj$cpg_row_in_chunk
    )

    out[, `:=`(
      n_train_total = sum(idx_train),
      n_train_complete = nrow(dat),
      n_train_FIDs = length(unique(dat$FID)),
      n_cases = sum(dat$AA_only == 1L),
      n_controls = sum(dat$AA_only == 0L),
      n_meth_observed = sum(is.finite(dat_dt$meth_response)),
      meth_sd = if (nrow(dat) > 1L) stats::sd(dat$meth_response) else NA_real_
    )]

    # Basic fit guards.
    if (nrow(dat) < min_complete_n) {
      out[, `:=`(
        fit_status = "SKIPPED_TOO_FEW_COMPLETE_SUBJECTS",
        error_message = paste0(
          "Complete training N=", nrow(dat),
          "; required >= ", min_complete_n
        )
      )]
      split_results[[ss]] <- out
      next
    }

    if (length(unique(dat$AA_only)) < 2L) {
      out[, `:=`(
        fit_status = "SKIPPED_ONE_AA_CLASS",
        error_message = "Training data contain only one AA class after complete-case filtering."
      )]
      split_results[[ss]] <- out
      next
    }

    if (!is.finite(stats::sd(dat$meth_response)) ||
        stats::sd(dat$meth_response) == 0) {
      out[, `:=`(
        fit_status = "SKIPPED_ZERO_METHYLATION_VARIANCE",
        error_message = "ASR methylation has zero/non-finite SD in this training split."
      )]
      split_results[[ss]] <- out
      next
    }

    if (length(unique(dat$FID)) < 2L) {
      out[, `:=`(
        fit_status = "SKIPPED_TOO_FEW_FAMILIES",
        error_message = "Fewer than two training families after complete-case filtering."
      )]
      split_results[[ss]] <- out
      next
    }

    fit_res <- tryCatch(
      fit_requested_model(as.character(info$model), dat),
      error = function(e) e
    )

    if (inherits(fit_res, "error")) {
      n_fit_error <- n_fit_error + 1L
      out[, `:=`(
        fit_status = "MODEL_ERROR",
        error_message = conditionMessage(fit_res)
      )]
      split_results[[ss]] <- out
      next
    }

    n_fit_ok <- n_fit_ok + 1L

    out[, `:=`(
      association_term = fit_res$association_term,
      coefficient = fit_res$coefficient,
      std_error = fit_res$std_error,
      statistic = fit_res$statistic,
      df = fit_res$df,
      p_value = fit_res$p_value,
      FID_sd = fit_res$FID_sd,
      optimal_lambda = fit_res$optimal_lambda,
      optimal_lambda_BIC = fit_res$optimal_lambda_BIC,
      optimal_lambda_AIC = fit_res$optimal_lambda_AIC,
      stage1_coefficient = fit_res$stage1_coefficient,
      stage1_p_value = fit_res$stage1_p_value,
      selected_stage1 = fit_res$selected_stage1,
      final_fit_stage = fit_res$final_fit_stage,
      fit_status = "FIT_OK",
      error_message = cpg_obj$message
    )]

    split_results[[ss]] <- out
  }

  cpg_results <- rbindlist(split_results, use.names = TRUE, fill = TRUE)

  file_exists_before <- file.exists(result_file)
  fwrite(
    cpg_results,
    result_file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    append = file_exists_before,
    col.names = !file_exists_before
  )

  n_output <- n_output + nrow(cpg_results)

  rm(cpg_results, split_results)
  if (j %% 10L == 0L) invisible(gc(FALSE))
}

# ----------------------------------- QC --------------------------------------

qc <- data.table(
  job_id = job_id,
  n_jobs = n_jobs,
  total_final_table_rows = n_total_rows,
  assigned_row_start = row_start,
  assigned_row_end = row_end,
  n_job_rows_requested = nrow(job_rows),
  n_rows_processed = n_rows_processed,
  n_split_ids = length(split_ids),
  expected_output_rows = nrow(job_rows) * length(split_ids),
  actual_output_rows = n_output,
  n_successful_model_fits = n_fit_ok,
  n_model_errors = n_fit_error,
  n_cpgs_not_found = n_cpg_not_found,
  result_file = result_file,
  elapsed_minutes = as.numeric(
    difftime(Sys.time(), run_start, units = "mins")
  ),
  completed_at = format(
    Sys.time(),
    tz = "America/Toronto",
    usetz = TRUE
  ),
  job_status = "COMPLETED"
)

fwrite(qc, qc_file, sep = "\t", quote = FALSE, na = "NA")

cat("\nCompleted job", job_id, "\n")
print(qc)

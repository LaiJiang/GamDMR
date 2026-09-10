#!/usr/bin/env Rscript

# =============================================================================
# 13_run_primary_single_cpg_training_associations.R
#
# PRIMARY single-CpG association analysis for predictive benchmarking.
#
# IMPORTANT DESIGN:
#   * Its CpG universe comes from the coordinate-only manifest built directly
#     from the methylation chunk files.
#   * For every CpG, M1/M2/M3 are re-estimated separately within each OUTER
#     TRAINING-family split.
#   * Outer-test families are never used in the association screen.
#
# Output column names/order and default result location intentionally match the
# pilot pipeline so these files can overwrite the pilot outputs.
#
# Default output:
#   .../6_single_cpg/1_training_associations/job_results/
#       single_cpg_training_assoc_job_XXXX.tsv
#   .../6_single_cpg/1_training_associations/job_qc/
#       single_cpg_training_assoc_job_XXXX_qc.tsv
#
# M1/M2/M3 reproduce the pilot implementation:
#   M1: logistic glmmLasso + BIC lambda + unpenalized glmer refit
#   M2: direct LMM of ASR methylation on AA + covariates + (1|FID)
#   M3: Gaussian glmmLasso + BIC lambda + unpenalized lmer refit
#
# Unselected target terms are reported as coefficient=0 and p_value=1,
# matching the pilot v3 behavior.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(glmmLasso)
})

data.table::setDTthreads(1L)
options(mc.cores = 1L)

# ------------------------------- configuration -------------------------------

PATH_WK <- path.expand(Sys.getenv("PATH_WK", "~/scratch/UQAC/meth/"))

manifest_root <- path.expand(Sys.getenv(
  "PRIMARY_CPG_MANIFEST_ROOT",
  file.path(
    PATH_WK,
    "results/15_revision/5_cv/6_single_cpg/0_primary_all_cpg_manifest"
  )
))

train_split_file <- path.expand(Sys.getenv(
  "TRAIN_SPLIT_FILE",
  file.path(PATH_WK, "scr/11_mgcv/dat/100_family_training_splits.csv")
))

pheno_file_path <- path.expand(Sys.getenv(
  "PHENO_FILE",
  file.path(PATH_WK, "scr/11_mgcv/dat/18_pheno_BMI.RData")
))

chunk_dir <- path.expand(Sys.getenv(
  "METH_CHUNK_DIR",
  file.path(PATH_WK, "data/meth_split")
))

# EXACTLY the same default output location as the pilot.
output_root <- path.expand(Sys.getenv(
  "SINGLE_CPG_ASSOC_OUTPUT",
  file.path(
    PATH_WK,
    "results/15_revision/5_cv/6_single_cpg/1_training_associations"
  )
))

n_jobs <- suppressWarnings(as.integer(Sys.getenv("N_JOBS", "990")))
min_complete_n <- suppressWarnings(as.integer(Sys.getenv("MIN_COMPLETE_N", "31")))
max_cpgs <- suppressWarnings(as.integer(Sys.getenv("MAX_CPGS", "0")))
max_splits <- suppressWarnings(as.integer(Sys.getenv("MAX_SPLITS", "0")))
progress_every <- suppressWarnings(as.integer(Sys.getenv("PROGRESS_EVERY", "25")))
write_buffer_cpgs <- suppressWarnings(as.integer(
  Sys.getenv("WRITE_BUFFER_CPGS", "10")
))

overwrite <- Sys.getenv("OVERWRITE", "1") %in%
  c("1", "TRUE", "true", "T", "yes", "YES")

models_env <- toupper(trimws(Sys.getenv("MODELS", "M1,M2,M3")))
models_to_run <- unique(trimws(strsplit(models_env, ",", fixed = TRUE)[[1L]]))
models_to_run <- models_to_run[nzchar(models_to_run)]

valid_models <- c("M1", "M2", "M3")
if (length(models_to_run) == 0L ||
    any(!models_to_run %in% valid_models)) {
  stop(
    "MODELS must be a comma-separated subset of M1,M2,M3. Observed: ",
    models_env
  )
}

split_ids_env <- trimws(Sys.getenv("SPLIT_IDS", ""))

lambda_grid <- c(5, 10, 20, 40, 50, 60, 80, 100)

if (is.na(n_jobs) || n_jobs < 1L) stop("N_JOBS must be >= 1.")
if (is.na(min_complete_n) || min_complete_n < 10L) {
  stop("MIN_COMPLETE_N must be >= 10.")
}
if (is.na(max_cpgs) || max_cpgs < 0L) max_cpgs <- 0L
if (is.na(max_splits) || max_splits < 0L) max_splits <- 0L
if (is.na(progress_every) || progress_every < 1L) progress_every <- 25L
if (is.na(write_buffer_cpgs) || write_buffer_cpgs < 1L) {
  write_buffer_cpgs <- 10L
}

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

job_manifest_file <- file.path(
  manifest_root,
  "job_manifests",
  sprintf("primary_cpg_manifest_job_%04d.tsv", job_id)
)

if (!file.exists(job_manifest_file)) {
  stop("Job manifest does not exist: ", job_manifest_file)
}

result_file <- file.path(
  result_dir,
  sprintf("single_cpg_training_assoc_job_%04d.tsv", job_id)
)

qc_file <- file.path(
  qc_dir,
  sprintf("single_cpg_training_assoc_job_%04d_qc.tsv", job_id)
)

cat("PRIMARY single-CpG training-split associations\n")
cat("job_id:", job_id, "of", n_jobs, "\n")
cat("manifest:", job_manifest_file, "\n")
cat("models:", paste(models_to_run, collapse = ","), "\n")
cat("output:", result_file, "\n")

if (!overwrite && file.exists(result_file) && file.exists(qc_file)) {
  qc_old <- tryCatch(fread(qc_file), error = function(e) NULL)

  if (!is.null(qc_old) &&
      "job_status" %in% names(qc_old) &&
      any(qc_old$job_status == "COMPLETED")) {
    cat("Completed outputs already exist; exiting.\n")
    quit(save = "no", status = 0L)
  }
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
  tryCatch({
    vc <- lme4::VarCorr(fit)
    if (!"FID" %in% names(vc)) return(NA_real_)
    as.numeric(attr(vc$FID, "stddev"))[1L]
  }, error = function(e) NA_real_)
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

  if (is.null(tab) ||
      is.null(rownames(tab)) ||
      !term %in% rownames(tab)) {
    return(out)
  }

  est_col <- find_col(tab, c("^Estimate$", "^Coef", "estimate"))
  se_col <- find_col(
    tab,
    c("Std\\. Error", "Std\\.Err", "Std\\.Error", "^SE$")
  )
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

  if (!is.null(co) &&
      !is.null(names(co)) &&
      term %in% names(co)) {
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

model_index <- function(model_label) {
  match(model_label, c("M1", "M2", "M3"))
}

logical_row_id <- function(primary_cpg_id, model_label) {
  # Unique compatibility ID corresponding to an implicit expansion:
  # CpG1:M1, CpG1:M2, CpG1:M3, CpG2:M1, ...
  as.numeric((primary_cpg_id - 1) * 3 + model_index(model_label))
}

make_base_result <- function(info, split_id, model_label) {
  data.table(
    job_id = job_id,
    final_table_row_id = logical_row_id(
      as.numeric(info$primary_cpg_id),
      model_label
    ),
    splitID = as.integer(split_id),
    CpG = as.character(info$CpG),
    requested_model = as.character(model_label),
    cpg_chr = as.character(info$cpg_chr),
    cpg_start = as.integer(info$cpg_start),

    # Same output columns as the pilot.
    data_chunk_id = as.integer(info$data_chunk_id),
    data_chunk_id_prev = NA_integer_,
    data_chunk_id_next = NA_integer_,
    resolved_chunk_id = as.integer(info$data_chunk_id),
    chunk_source = "primary_raw_chunk",
    cpg_row_in_chunk = as.integer(info$cpg_row_in_chunk),

    # Deliberately NA: the primary universe was NOT selected using full-data
    # associations and no full-data p-value/coef is imported.
    original_full_data_pval = NA_real_,
    original_full_data_coef = NA_real_,
    genes_entrez = NA_character_,
    genes_symbol = NA_character_,

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

append_results <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0L) return(invisible(NULL))

  exists_before <- file.exists(result_file)

  fwrite(
    dt,
    result_file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    append = exists_before,
    col.names = !exists_before
  )

  invisible(NULL)
}

# ----------------------------- model functions -------------------------------

covariates <- c(
  "AgeCalc", "Sex", "Non.smoker",
  "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc",
  "sv1", "sv2", "sv3", "sv4", "sv5", "BMI"
)

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
      paste0(
        "AA_only ~ ",
        paste(selected, collapse = " + "),
        " + (1|FID)"
      )
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
      final <- list(
        estimate = 0,
        std_error = NA_real_,
        statistic = NA_real_,
        df = NA_real_,
        p_value = 1
      )
      fid_sd <- 0
      final_stage <- paste0(
        "glmmLasso_stage1_refit_failed: ",
        conditionMessage(refit)
      )
    }
  } else {
    final_stage <- "glmmLasso_stage1_no_terms_selected"
  }

  if (!selected_target) {
    final <- list(
      estimate = 0,
      std_error = NA_real_,
      statistic = NA_real_,
      df = NA_real_,
      p_value = 1
    )
  }

  stage1_p_out <- stage1$p_value
  if (!is.finite(stage1_p_out) && !selected_target) stage1_p_out <- 1

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
    stage1_p_value = stage1_p_out,
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
      final <- list(
        estimate = 0,
        std_error = NA_real_,
        statistic = NA_real_,
        df = NA_real_,
        p_value = 1
      )
      fid_sd <- 0
      final_stage <- paste0(
        "glmmLasso_stage1_refit_failed: ",
        conditionMessage(refit)
      )
    }
  } else {
    final_stage <- "glmmLasso_stage1_no_terms_selected"
  }

  if (!selected_target) {
    final <- list(
      estimate = 0,
      std_error = NA_real_,
      statistic = NA_real_,
      df = NA_real_,
      p_value = 1
    )
  }

  stage1_p_out <- stage1$p_value
  if (!is.finite(stage1_p_out) && !selected_target) stage1_p_out <- 1

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
    stage1_p_value = stage1_p_out,
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

# ------------------------- load job CpG manifest -----------------------------

job_cpgs <- fread(job_manifest_file)

required_manifest_cols <- c(
  "primary_cpg_id",
  "CpG",
  "cpg_chr",
  "cpg_start",
  "cpg_end",
  "data_chunk_id",
  "cpg_row_in_chunk"
)

missing_manifest <- setdiff(required_manifest_cols, names(job_cpgs))

if (length(missing_manifest) > 0L) {
  stop(
    "Job manifest is missing: ",
    paste(missing_manifest, collapse = ", ")
  )
}

job_cpgs[, cpg_chr := normalize_chr(cpg_chr)]
job_cpgs[, cpg_start := suppressWarnings(as.integer(cpg_start))]
job_cpgs[, cpg_end := suppressWarnings(as.integer(cpg_end))]
job_cpgs[, data_chunk_id := suppressWarnings(as.integer(data_chunk_id))]
job_cpgs[, cpg_row_in_chunk := suppressWarnings(as.integer(cpg_row_in_chunk))]

if (max_cpgs > 0L && nrow(job_cpgs) > max_cpgs) {
  job_cpgs <- head(job_cpgs, max_cpgs)
}

cat("Assigned primary CpGs:", nrow(job_cpgs), "\n")

# ----------------------- load split and phenotype data -----------------------

if (!file.exists(train_split_file)) {
  stop("Training-split file does not exist: ", train_split_file)
}

train_split_dt <- fread(train_split_file)

train_cols <- grep("^train_FID_[0-9]+$", names(train_split_dt), value = TRUE)

if (!"splitID" %in% names(train_split_dt) ||
    length(train_cols) == 0L) {
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

if (!file.exists(pheno_file_path)) {
  stop("Phenotype file does not exist: ", pheno_file_path)
}

pe <- new.env(parent = emptyenv())
load(pheno_file_path, envir = pe)

if (!exists("pheno_file", envir = pe, inherits = FALSE)) {
  stop("Phenotype RData does not contain pheno_file.")
}

pheno_dt <- as.data.table(
  copy(get("pheno_file", envir = pe, inherits = FALSE))
)
rm(pe)

required_pheno <- c("ID", "FID", "AA_only", covariates)
missing_pheno <- setdiff(required_pheno, names(pheno_dt))

if (length(missing_pheno) > 0L) {
  stop(
    "pheno_file is missing: ",
    paste(missing_pheno, collapse = ", ")
  )
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

# Precompute training phenotype subsets once. This avoids repeating FID
# membership operations millions of times.
train_base_list <- vector("list", length(split_ids))
names(train_base_list) <- as.character(split_ids)

for (ss in seq_along(split_ids)) {
  sid <- split_ids[ss]

  rr <- train_split_dt[splitID == sid]

  if (nrow(rr) != 1L) {
    stop("Expected one row for splitID ", sid)
  }

  train_fids <- trimws(
    as.character(unlist(rr[, ..train_cols], use.names = FALSE))
  )

  train_fids <- unique(
    train_fids[
      !is.na(train_fids) &
        nzchar(train_fids) &
        train_fids != "NA"
    ]
  )

  idx <- which(pheno_dt$FID %chin% train_fids)

  base <- copy(pheno_dt[idx])
  base[, pheno_row_index__ := idx]

  # Covariate/outcome completeness is independent of the CpG.
  base_required <- c("AA_only", "FID", covariates)
  base[, base_complete__ := complete.cases(base[, ..base_required])]

  train_base_list[[ss]] <- base
}

cat("Training splits loaded:", length(split_ids), "\n")

# --------------------------- methylation chunk cache -------------------------

chunk_cache <- new.env(parent = emptyenv())
cache_order <- integer()
cache_chunks <- 2L

get_chunk <- function(chunk_id) {
  key <- as.character(as.integer(chunk_id))

  if (exists(key, envir = chunk_cache, inherits = FALSE)) {
    cache_order <<- c(
      cache_order[cache_order != as.integer(chunk_id)],
      as.integer(chunk_id)
    )

    return(get(key, envir = chunk_cache, inherits = FALSE))
  }

  chunk_file <- file.path(
    chunk_dir,
    sprintf("chunk_%04d.csv", as.integer(chunk_id))
  )

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

  select_cols <- c("chr", "start")
  if ("end" %in% hnames) select_cols <- c(select_cols, "end")
  select_cols <- c(select_cols, meth_cols)

  dt <- tryCatch(
    fread(
      chunk_file,
      select = select_cols,
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

  if (!"end" %in% names(dt)) dt[, end := NA_integer_]

  dt[, chr := normalize_chr(chr)]
  dt[, start := suppressWarnings(as.integer(start))]
  dt[, end := suppressWarnings(as.integer(end))]

  obj <- list(
    dt = dt,
    meth_cols = meth_cols,
    meth_ids = meth_ids,
    chunk_file = chunk_file
  )

  assign(key, obj, envir = chunk_cache)

  cache_order <<- c(
    cache_order[cache_order != as.integer(chunk_id)],
    as.integer(chunk_id)
  )

  while (length(cache_order) > cache_chunks) {
    old <- cache_order[1L]
    cache_order <<- cache_order[-1L]
    rm(list = as.character(old), envir = chunk_cache)
    invisible(gc(FALSE))
  }

  obj
}

extract_cpg_methylation <- function(info) {
  chunk_id <- as.integer(info$data_chunk_id)
  obj <- get_chunk(chunk_id)

  if (inherits(obj, "chunk_error")) {
    return(list(
      status = "CHUNK_ERROR",
      message = obj$message,
      meth = NULL
    ))
  }

  rr <- as.integer(info$cpg_row_in_chunk)

  use_rr <- is.finite(rr) &&
    rr >= 1L &&
    rr <= nrow(obj$dt) &&
    identical(
      as.character(obj$dt$chr[rr]),
      as.character(info$cpg_chr)
    ) &&
    !is.na(obj$dt$start[rr]) &&
    obj$dt$start[rr] == as.integer(info$cpg_start)

  if (!use_rr) {
    hit <- which(
      obj$dt$chr == as.character(info$cpg_chr) &
        obj$dt$start == as.integer(info$cpg_start)
    )

    if (length(hit) == 0L) {
      return(list(
        status = "CPG_NOT_FOUND_IN_CHUNK",
        message = paste0(
          "CpG not found at ",
          info$cpg_chr, ":", info$cpg_start,
          " in chunk ", chunk_id
        ),
        meth = NULL
      ))
    }

    rr <- hit[1L]
  }

  source_cols <- obj$meth_cols

  vals_mat <- as.matrix(
    obj$dt[rr, ..source_cols]
  )
  storage.mode(vals_mat) <- "numeric"

  vals <- as.numeric(vals_mat[1L, ])
  vals[!is.finite(vals) | vals < 0 | vals > 1] <- NA_real_
  names(vals) <- obj$meth_ids

  meth_match_idx <- match(pheno_dt$ID, names(vals))
  meth_by_pheno <- vals[meth_match_idx]
  meth_by_pheno[!is.finite(meth_by_pheno)] <- NA_real_

  list(
    status = "FOUND",
    message = NA_character_,
    meth = meth_by_pheno
  )
}

# ------------------------------ analysis loop --------------------------------

n_output <- 0L
n_fit_ok <- 0L
n_fit_error <- 0L
n_cpg_not_found <- 0L
n_cpgs_processed <- 0L

result_buffer <- list()
buffer_cpg_count <- 0L

flush_buffer <- function() {
  if (length(result_buffer) == 0L) return(invisible(NULL))

  dt <- rbindlist(
    result_buffer,
    use.names = TRUE,
    fill = TRUE
  )

  append_results(dt)

  n_output <<- n_output + nrow(dt)
  result_buffer <<- list()
  buffer_cpg_count <<- 0L

  invisible(gc(FALSE))
}

for (jj in seq_len(nrow(job_cpgs))) {
  info <- job_cpgs[jj]
  n_cpgs_processed <- n_cpgs_processed + 1L

  if (jj == 1L ||
      jj %% progress_every == 0L ||
      jj == nrow(job_cpgs)) {
    cat(
      sprintf(
        "[job %04d] CpG %d/%d; primary_cpg_id=%s; %s; chunk=%d row=%d\n",
        job_id,
        jj,
        nrow(job_cpgs),
        as.character(info$primary_cpg_id),
        info$CpG,
        info$data_chunk_id,
        info$cpg_row_in_chunk
      )
    )
  }

  cpg_obj <- extract_cpg_methylation(info)

  cpg_rows <- list()

  if (is.null(cpg_obj$meth)) {
    n_cpg_not_found <- n_cpg_not_found + 1L

    kk <- 0L

    for (sid in split_ids) {
      for (model_label in models_to_run) {
        kk <- kk + 1L

        out <- make_base_result(info, sid, model_label)
        out[, `:=`(
          fit_status = cpg_obj$status,
          error_message = cpg_obj$message
        )]

        cpg_rows[[kk]] <- out
      }
    }

    result_buffer[[length(result_buffer) + 1L]] <- rbindlist(
      cpg_rows,
      use.names = TRUE,
      fill = TRUE
    )

    buffer_cpg_count <- buffer_cpg_count + 1L

    if (buffer_cpg_count >= write_buffer_cpgs) flush_buffer()
    next
  }

  meth_by_pheno <- cpg_obj$meth
  meth_asr <- asin(sqrt(meth_by_pheno))

  kk <- 0L

  for (ss in seq_along(split_ids)) {
    sid <- split_ids[ss]
    base <- train_base_list[[ss]]

    meth_train <- meth_asr[base$pheno_row_index__]

    complete_idx <- base$base_complete__ &
      is.finite(meth_train)

    dat_dt <- copy(base[complete_idx])
    dat_dt[, meth_response := meth_train[complete_idx]]

    dat_dt[, c("pheno_row_index__", "base_complete__") := NULL]

    dat <- as.data.frame(dat_dt)
    dat$FID <- factor(dat$FID)

    n_train_total <- nrow(base)
    n_train_complete <- nrow(dat)
    n_train_FIDs <- length(unique(dat$FID))
    n_cases <- sum(dat$AA_only == 1L)
    n_controls <- sum(dat$AA_only == 0L)
    n_meth_observed <- sum(is.finite(meth_train))
    meth_sd <- if (nrow(dat) > 1L) {
      stats::sd(dat$meth_response)
    } else {
      NA_real_
    }

    for (model_label in models_to_run) {
      kk <- kk + 1L
      out <- make_base_result(info, sid, model_label)

      out[, `:=`(
        n_train_total = as.integer(n_train_total),
        n_train_complete = as.integer(n_train_complete),
        n_train_FIDs = as.integer(n_train_FIDs),
        n_cases = as.integer(n_cases),
        n_controls = as.integer(n_controls),
        n_meth_observed = as.integer(n_meth_observed),
        meth_sd = meth_sd
      )]

      if (nrow(dat) < min_complete_n) {
        out[, `:=`(
          fit_status = "SKIPPED_TOO_FEW_COMPLETE_SUBJECTS",
          error_message = paste0(
            "Complete training N=", nrow(dat),
            "; required >= ", min_complete_n
          )
        )]

        cpg_rows[[kk]] <- out
        next
      }

      if (length(unique(dat$AA_only)) < 2L) {
        out[, `:=`(
          fit_status = "SKIPPED_ONE_AA_CLASS",
          error_message = paste0(
            "Training data contain one AA class after complete-case filtering."
          )
        )]

        cpg_rows[[kk]] <- out
        next
      }

      if (!is.finite(meth_sd) || meth_sd == 0) {
        out[, `:=`(
          fit_status = "SKIPPED_ZERO_METHYLATION_VARIANCE",
          error_message = "ASR methylation has zero/non-finite SD."
        )]

        cpg_rows[[kk]] <- out
        next
      }

      if (n_train_FIDs < 2L) {
        out[, `:=`(
          fit_status = "SKIPPED_TOO_FEW_FAMILIES",
          error_message = "Fewer than two training families."
        )]

        cpg_rows[[kk]] <- out
        next
      }

      fit_res <- tryCatch(
        fit_requested_model(model_label, dat),
        error = function(e) e
      )

      if (inherits(fit_res, "error")) {
        n_fit_error <- n_fit_error + 1L

        out[, `:=`(
          fit_status = "MODEL_ERROR",
          error_message = conditionMessage(fit_res)
        )]

        cpg_rows[[kk]] <- out
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
        error_message = NA_character_
      )]

      cpg_rows[[kk]] <- out
    }
  }

  result_buffer[[length(result_buffer) + 1L]] <- rbindlist(
    cpg_rows,
    use.names = TRUE,
    fill = TRUE
  )

  buffer_cpg_count <- buffer_cpg_count + 1L

  if (buffer_cpg_count >= write_buffer_cpgs) {
    flush_buffer()
  }
}

flush_buffer()

# ----------------------------------- QC --------------------------------------

expected_output_rows <- (
  nrow(job_cpgs) *
    length(split_ids) *
    length(models_to_run)
)

qc <- data.table(
  job_id = job_id,
  n_jobs = n_jobs,

  # Keep pilot QC columns in the same order/name where possible.
  total_final_table_rows = NA_integer_,
  assigned_row_start = if (nrow(job_cpgs) > 0L) {
    min(job_cpgs$primary_cpg_id)
  } else {
    NA_real_
  },
  assigned_row_end = if (nrow(job_cpgs) > 0L) {
    max(job_cpgs$primary_cpg_id)
  } else {
    NA_real_
  },
  n_job_rows_requested = nrow(job_cpgs),
  n_rows_processed = n_cpgs_processed,
  n_split_ids = length(split_ids),
  expected_output_rows = expected_output_rows,
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
  job_status = if (n_output == expected_output_rows) {
    "COMPLETED"
  } else {
    "COMPLETED_ROWCOUNT_MISMATCH"
  }
)

fwrite(
  qc,
  qc_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\nCompleted PRIMARY job", job_id, "\n")
print(qc)

if (n_output != expected_output_rows) {
  warning(
    "Output row count mismatch: expected ",
    expected_output_rows,
    ", observed ",
    n_output
  )
}

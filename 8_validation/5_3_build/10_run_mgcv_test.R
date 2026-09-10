#!/usr/bin/env Rscript

#this file test run different strategies of prediction with mgcv DMRs.!!!!!!!!!!!!!!!!!!!

# Primary predictive evaluation for split-specific MGCV/GAM-DMR selections.
#
# One SLURM array task processes one outer family split:
#   1. Load the DMRs selected using that split's training families.
#   2. Construct DMR PC1/PC2 features using training data only.
#   3. Project the corresponding untouched test families into the training PCA.
#   4. Tune logistic LASSO by family-grouped fivefold inner CV.
#   5. Evaluate the untouched outer test set.
#   6. Fit and evaluate a covariate-only reference model.
#
# Primary methylation model:
#   MGCV DMR PC1/PC2 + age + sex + smoking + BMI
#   Logistic LASSO, alpha = 1, lambda = lambda.1se
#   DMR features penalized; baseline covariates unpenalized.
#
# Critical leakage controls:
#   - DMR selection is split-specific and training-derived.
#   - CpG filtering and imputation use training subjects only.
#   - PCA center, scale, and rotation use training subjects only.
#   - Test subjects are projected with the training PCA.
#   - Inner folds keep all members of a family together.
#   - The classification threshold is chosen from inner-CV predictions only.

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(pROC)
  library(PRROC)
})

data.table::setDTthreads(1L)
options(mc.cores = 1L)


lambda_rule <- Sys.getenv(
  "LAMBDA_RULE",
  "lambda.min"
)

if (!lambda_rule %in% c("lambda.1se", "lambda.min")) {
  stop(
    "LAMBDA_RULE must be either lambda.1se or lambda.min."
  )
}

# ------------------------------- configuration -------------------------------

PATH_wk <- path.expand(
  Sys.getenv("PATH_WK", "~/scratch/UQAC/meth/")
)

selected_dmr_dir <- Sys.getenv(
  "MGCV_SPLIT_DMR_DIR",
  file.path(
    PATH_wk,
    "results/15_revision/5_cv/1_mgcv/collected_by_split/strict_by_split"
  )
)

train_split_file <- Sys.getenv(
  "TRAIN_SPLIT_FILE",
  file.path(
    PATH_wk,
    "scr/11_mgcv/dat/100_family_training_splits.csv"
  )
)

pheno_file_path <- Sys.getenv(
  "PHENO_FILE",
  file.path(
    PATH_wk,
    "scr/11_mgcv/dat/18_pheno_BMI.RData"
  )
)

chunk_dir <- Sys.getenv(
  "METH_CHUNK_DIR",
  file.path(PATH_wk, "data/meth_split")
)

output_root <- Sys.getenv(
  "MGCV_PREDICTION_OUTPUT",
  file.path(
    PATH_wk,
    "results/15_revision/5_cv/5_prediction/1_mgcv_primary"
  )
)

performance_dir <- file.path(output_root, "performance")
prediction_dir <- file.path(output_root, "predictions")
coefficient_dir <- file.path(output_root, "coefficients")
feature_qc_dir <- file.path(output_root, "feature_qc")
split_qc_dir <- file.path(output_root, "split_qc")
shared_fold_dir <- file.path(
  PATH_wk,
  "results/15_revision/5_cv/5_prediction/shared_inner_folds"
)

for (d in c(
  performance_dir,
  prediction_dir,
  coefficient_dir,
  feature_qc_dir,
  split_qc_dir,
  shared_fold_dir
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

inner_folds_requested <- suppressWarnings(
  as.integer(Sys.getenv("INNER_FOLDS", "5"))
)

seed_base <- suppressWarnings(
  as.integer(Sys.getenv("SEED_BASE", "20260803"))
)

min_train_cpg_observed_fraction <- suppressWarnings(
  as.numeric(Sys.getenv("MIN_TRAIN_CPG_OBSERVED_FRACTION", "0.50"))
)

min_variable_cpgs <- suppressWarnings(
  as.integer(Sys.getenv("MIN_VARIABLE_CPGS", "1"))
)

pca_components_requested <- suppressWarnings(
  as.integer(Sys.getenv("N_DMR_PCS", "2"))
)

max_dmrs <- suppressWarnings(
  as.integer(Sys.getenv("MAX_DMRS", "0"))
)

variance_epsilon <- suppressWarnings(
  as.numeric(Sys.getenv("VARIANCE_EPSILON", "1e-8"))
)

overwrite <- Sys.getenv("OVERWRITE", "0") %in%
  c("1", "TRUE", "true", "T", "yes", "YES")

if (is.na(inner_folds_requested) || inner_folds_requested < 3L) {
  stop("INNER_FOLDS must be at least 3.")
}

if (!is.finite(min_train_cpg_observed_fraction) ||
    min_train_cpg_observed_fraction <= 0 ||
    min_train_cpg_observed_fraction > 1) {
  stop("MIN_TRAIN_CPG_OBSERVED_FRACTION must be in (0, 1].")
}

if (is.na(min_variable_cpgs) || min_variable_cpgs < 1L) {
  stop("MIN_VARIABLE_CPGS must be at least 1.")
}

if (is.na(pca_components_requested) || pca_components_requested < 1L) {
  stop("N_DMR_PCS must be at least 1.")
}

if (is.na(max_dmrs) || max_dmrs < 0L) {
  stop("MAX_DMRS must be nonnegative.")
}

if (!is.finite(variance_epsilon) || variance_epsilon < 0) {
  stop("VARIANCE_EPSILON must be nonnegative.")
}

# -------------------------------- arguments ----------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1L) {
  stop("Pass splitID or SLURM_ARRAY_TASK_ID as the first argument.")
}

split_id <- suppressWarnings(as.integer(args[1L]))

if (is.na(split_id) || split_id < 1L || split_id > 100L) {
  stop("splitID must be an integer from 1 to 100.")
}

cat("MGCV predictive analysis: outer split", split_id, "\n")
run_start <- Sys.time()

# --------------------------------- outputs -----------------------------------

performance_file <- file.path(
  performance_dir,
  sprintf("mgcv_primary_performance_split_%03d.tsv", split_id)
)

prediction_file <- file.path(
  prediction_dir,
  sprintf("mgcv_primary_predictions_split_%03d.tsv", split_id)
)

coefficient_file <- file.path(
  coefficient_dir,
  sprintf("mgcv_primary_coefficients_split_%03d.tsv", split_id)
)

feature_qc_file <- file.path(
  feature_qc_dir,
  sprintf("mgcv_primary_feature_qc_split_%03d.tsv", split_id)
)

split_qc_file <- file.path(
  split_qc_dir,
  sprintf("mgcv_primary_split_qc_%03d.tsv", split_id)
)

fold_file <- file.path(
  shared_fold_dir,
  sprintf("inner_family_folds_split_%03d.tsv", split_id)
)

output_files <- c(
  performance_file,
  prediction_file,
  coefficient_file,
  feature_qc_file,
  split_qc_file
)

if (!overwrite && all(file.exists(output_files))) {
  cat("All split outputs already exist; exiting without replacement.\n")
  quit(save = "no", status = 0L)
}

if (overwrite) {
  unlink(output_files, force = TRUE)
}

# -------------------------------- utilities ----------------------------------

normalize_id <- function(x) {
  x <- as.character(x)
  x <- sub("^X", "", x)
  gsub("\\.", "-", x)
}

normalize_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  sub("\\.0$", "", x)
}

clean_probability <- function(p, eps = 1e-8) {
  p <- as.numeric(p)
  pmin(pmax(p, eps), 1 - eps)
}

safe_divide <- function(a, b) {
  if (!is.finite(b) || b == 0) return(NA_real_)
  as.numeric(a / b)
}

safe_roc_auc <- function(y, p) {
  y <- as.integer(y)
  p <- as.numeric(p)
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]
  p <- p[keep]

  if (length(unique(y)) < 2L) return(NA_real_)

  roc_obj <- tryCatch(
    pROC::roc(
      response = y,
      predictor = p,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(roc_obj)) return(NA_real_)
  as.numeric(pROC::auc(roc_obj))
}

safe_pr_auc <- function(y, p) {
  y <- as.integer(y)
  p <- as.numeric(p)
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]
  p <- p[keep]

  if (length(unique(y)) < 2L) return(NA_real_)

  pr_obj <- tryCatch(
    PRROC::pr.curve(
      scores.class0 = p[y == 1L],
      scores.class1 = p[y == 0L],
      curve = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(pr_obj)) return(NA_real_)
  as.numeric(pr_obj$auc.integral)
}

choose_youden_threshold <- function(y, p) {
  y <- as.integer(y)
  p <- as.numeric(p)
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]
  p <- p[keep]

  if (length(unique(y)) < 2L || length(unique(p)) < 2L) {
    return(0.5)
  }

  roc_obj <- tryCatch(
    pROC::roc(
      response = y,
      predictor = p,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(roc_obj)) return(0.5)

  threshold <- tryCatch(
    pROC::coords(
      roc_obj,
      x = "best",
      best.method = "youden",
      ret = "threshold",
      transpose = FALSE
    ),
    error = function(e) NULL
  )

  threshold <- suppressWarnings(as.numeric(threshold)[1L])

  if (!is.finite(threshold)) 0.5 else threshold
}

calculate_metrics <- function(y, p, threshold) {
  y <- as.integer(y)
  p <- as.numeric(p)
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]
  p <- p[keep]

  if (length(y) == 0L) {
    empty_names <- c(
      "n", "n_case", "n_control", "prevalence", "AUROC",
      "PR_AUC", "PR_baseline", "PR_AUC_lift", "Brier", "log_loss",
      "calibration_intercept", "calibration_slope", "sensitivity",
      "specificity", "balanced_accuracy", "PPV", "NPV"
    )
    out <- as.list(rep(NA_real_, length(empty_names)))
    names(out) <- empty_names
    return(out)
  }

  pred_class <- as.integer(p >= threshold)

  tp <- sum(pred_class == 1L & y == 1L)
  tn <- sum(pred_class == 0L & y == 0L)
  fp <- sum(pred_class == 1L & y == 0L)
  fn <- sum(pred_class == 0L & y == 1L)

  sensitivity <- safe_divide(tp, tp + fn)
  specificity <- safe_divide(tn, tn + fp)
  ppv <- safe_divide(tp, tp + fp)
  npv <- safe_divide(tn, tn + fn)
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)

  if (!is.finite(balanced_accuracy)) balanced_accuracy <- NA_real_

  p_clip <- clean_probability(p)
  brier <- mean((p - y)^2)
  log_loss <- -mean(y * log(p_clip) + (1 - y) * log(1 - p_clip))

  calibration_intercept <- NA_real_
  calibration_slope <- NA_real_
  lp <- qlogis(p_clip)

  sd_lp <- stats::sd(lp)

  if (length(unique(y)) >= 2L && is.finite(sd_lp) && sd_lp > 0) {
    intercept_fit <- tryCatch(
      suppressWarnings(
        stats::glm(
          y ~ 1 + offset(lp),
          family = stats::binomial()
        )
      ),
      error = function(e) NULL
    )

    slope_fit <- tryCatch(
      suppressWarnings(
        stats::glm(
          y ~ lp,
          family = stats::binomial()
        )
      ),
      error = function(e) NULL
    )

    if (!is.null(intercept_fit)) {
      calibration_intercept <- unname(stats::coef(intercept_fit)[1L])
    }

    if (!is.null(slope_fit) && length(stats::coef(slope_fit)) >= 2L) {
      calibration_slope <- unname(stats::coef(slope_fit)[2L])
    }
  }

  list(
    n = length(y),
    n_case = sum(y == 1L),
    n_control = sum(y == 0L),
    prevalence = mean(y),
    AUROC = safe_roc_auc(y, p),
    PR_AUC = safe_pr_auc(y, p),
    PR_baseline = mean(y),
    PR_AUC_lift = safe_pr_auc(y, p) - mean(y),
    Brier = brier,
    log_loss = log_loss,
    calibration_intercept = calibration_intercept,
    calibration_slope = calibration_slope,
    sensitivity = sensitivity,
    specificity = specificity,
    balanced_accuracy = balanced_accuracy,
    PPV = ppv,
    NPV = npv
  )
}

metrics_to_columns <- function(metrics, suffix) {
  out <- as.list(metrics)
  metric_names <- names(out)

  # Keep the explicit outer-split sample-count columns n_train and n_test
  # unique. The metric named "n" is the number of evaluable observations
  # after removing non-finite outcomes/predictions.
  names(out) <- ifelse(
    metric_names == "n",
    paste0("n_evaluable_", suffix),
    paste0(metric_names, "_", suffix)
  )

  out
}

most_frequent <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(x) == 0L) return(NA_character_)
  tab <- sort(table(as.character(x)), decreasing = TRUE)
  names(tab)[1L]
}

prepare_numeric_covariate <- function(train_x, test_x, name) {
  train_num <- suppressWarnings(as.numeric(as.character(train_x)))
  test_num <- suppressWarnings(as.numeric(as.character(test_x)))

  med <- stats::median(train_num[is.finite(train_num)], na.rm = TRUE)
  if (!is.finite(med)) med <- 0

  train_num[!is.finite(train_num)] <- med
  test_num[!is.finite(test_num)] <- med

  list(
    train = matrix(train_num, ncol = 1L, dimnames = list(NULL, name)),
    test = matrix(test_num, ncol = 1L, dimnames = list(NULL, name)),
    imputation_value = med
  )
}

prepare_sex_covariate <- function(train_x, test_x) {
  train_char <- as.character(train_x)
  test_char <- as.character(test_x)

  train_num <- suppressWarnings(as.numeric(train_char))
  test_num <- suppressWarnings(as.numeric(test_char))

  nonmissing_train <- !is.na(train_char) & nzchar(train_char)
  numeric_possible <- all(
    !nonmissing_train |
      is.finite(train_num)
  )

  if (numeric_possible) {
    return(
      prepare_numeric_covariate(
        train_x = train_num,
        test_x = test_num,
        name = "COV_Sex"
      )
    )
  }

  mode_level <- most_frequent(train_char)
  if (is.na(mode_level)) mode_level <- "Unknown"

  train_char[is.na(train_char) | !nzchar(train_char)] <- mode_level
  test_char[is.na(test_char) | !nzchar(test_char)] <- mode_level

  train_levels <- sort(unique(train_char))
  test_char[!test_char %in% train_levels] <- mode_level

  if (length(train_levels) <= 1L) {
    return(
      list(
        train = matrix(numeric(0), nrow = length(train_char), ncol = 0L),
        test = matrix(numeric(0), nrow = length(test_char), ncol = 0L),
        imputation_value = mode_level
      )
    )
  }

  reference <- train_levels[1L]
  dummy_levels <- train_levels[-1L]

  train_matrix <- sapply(
    dummy_levels,
    function(z) as.numeric(train_char == z)
  )

  test_matrix <- sapply(
    dummy_levels,
    function(z) as.numeric(test_char == z)
  )

  if (is.null(dim(train_matrix))) {
    train_matrix <- matrix(train_matrix, ncol = 1L)
    test_matrix <- matrix(test_matrix, ncol = 1L)
  }

  colnames(train_matrix) <- paste0("COV_Sex_", make.names(dummy_levels))
  colnames(test_matrix) <- colnames(train_matrix)

  list(
    train = train_matrix,
    test = test_matrix,
    imputation_value = paste0("mode=", mode_level, ";reference=", reference)
  )
}

make_covariate_matrices <- function(train_dt, test_dt) {
  age <- prepare_numeric_covariate(
    train_dt$AgeCalc,
    test_dt$AgeCalc,
    "COV_AgeCalc"
  )

  smoking <- prepare_numeric_covariate(
    train_dt$Non.smoker,
    test_dt$Non.smoker,
    "COV_NonSmoker"
  )

  bmi <- prepare_numeric_covariate(
    train_dt$BMI,
    test_dt$BMI,
    "COV_BMI"
  )

  sex <- prepare_sex_covariate(
    train_dt$Sex,
    test_dt$Sex
  )

  train_matrix <- do.call(
    cbind,
    list(age$train, sex$train, smoking$train, bmi$train)
  )

  test_matrix <- do.call(
    cbind,
    list(age$test, sex$test, smoking$test, bmi$test)
  )

  rownames(train_matrix) <- train_dt$ID
  rownames(test_matrix) <- test_dt$ID

  list(
    train = train_matrix,
    test = test_matrix,
    imputation = data.table(
      covariate = c("AgeCalc", "Sex", "Non.smoker", "BMI"),
      training_imputation = as.character(
        c(
          age$imputation_value,
          sex$imputation_value,
          smoking$imputation_value,
          bmi$imputation_value
        )
      )
    )
  )
}

make_family_inner_folds <- function(
  train_dt,
  requested_k,
  seed,
  fold_path
) {
  if (file.exists(fold_path)) {
    existing <- fread(fold_path)

    required <- c("splitID", "FID", "inner_fold")
    if (!all(required %in% names(existing))) {
      stop("Existing inner-fold file has incompatible columns: ", fold_path)
    }

    existing[, FID := as.character(FID)]
    observed_fids <- sort(unique(as.character(train_dt$FID)))

    if (!setequal(existing$FID, observed_fids)) {
      stop(
        "Existing inner-fold FIDs do not match this split's training FIDs: ",
        fold_path
      )
    }

    foldid <- existing$inner_fold[
      match(as.character(train_dt$FID), existing$FID)
    ]

    return(
      list(
        foldid = as.integer(foldid),
        family_folds = existing,
        n_folds = max(existing$inner_fold)
      )
    )
  }

  fam_dt <- train_dt[
    ,
    .(
      fam_AA = max(AA_only, na.rm = TRUE),
      n_subjects = .N,
      n_cases = sum(AA_only == 1L),
      n_controls = sum(AA_only == 0L)
    ),
    by = FID
  ]

  fam_dt[, FID := as.character(FID)]

  max_k <- min(
    requested_k,
    fam_dt[fam_AA == 1L, .N],
    fam_dt[fam_AA == 0L, .N]
  )

  if (!is.finite(max_k) || max_k < 3L) {
    stop("Too few case/control families for at least three inner folds.")
  }

  selected <- NULL

  for (k in seq.int(max_k, 3L, by = -1L)) {
    for (attempt in seq_len(200L)) {
      set.seed(seed + 1000L * k + attempt)
      candidate <- copy(fam_dt)
      candidate[, inner_fold := NA_integer_]

      for (cls in c(0L, 1L)) {
        idx <- which(candidate$fam_AA == cls)
        idx <- sample(idx, length(idx), replace = FALSE)
        candidate$inner_fold[idx] <- rep(seq_len(k), length.out = length(idx))
      }

      subject_fold <- candidate$inner_fold[
        match(as.character(train_dt$FID), candidate$FID)
      ]

      fold_class_ok <- all(
        vapply(
          seq_len(k),
          function(f) {
            length(unique(train_dt$AA_only[subject_fold == f])) == 2L
          },
          logical(1L)
        )
      )

      if (fold_class_ok) {
        selected <- candidate
        break
      }
    }

    if (!is.null(selected)) break
  }

  if (is.null(selected)) {
    stop("Could not construct family-grouped inner folds containing both classes.")
  }

  fold_out <- selected[
    ,
    .(
      splitID = split_id,
      FID,
      fam_AA,
      n_subjects,
      n_cases,
      n_controls,
      inner_fold
    )
  ]

  setorder(fold_out, inner_fold, FID)

  fwrite(
    fold_out,
    fold_path,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  foldid <- fold_out$inner_fold[
    match(as.character(train_dt$FID), fold_out$FID)
  ]

  list(
    foldid = as.integer(foldid),
    family_folds = fold_out,
    n_folds = max(fold_out$inner_fold)
  )
}

# ------------------------------- input data ----------------------------------

selected_dmr_file <- file.path(
  selected_dmr_dir,
  sprintf("24_dmrs_STRICT_split_%03d.tsv", split_id)
)

if (!file.exists(selected_dmr_file)) {
  stop("Split-specific MGCV DMR file does not exist: ", selected_dmr_file)
}

if (!file.exists(train_split_file)) {
  stop("Training-split file does not exist: ", train_split_file)
}

if (!file.exists(pheno_file_path)) {
  stop("Phenotype file does not exist: ", pheno_file_path)
}

selected_dmrs <- fread(selected_dmr_file)

required_dmr_columns <- c(
  "splitID",
  "region_id",
  "data_chunk_id",
  "chr",
  "region_start",
  "region_end"
)

missing_dmr_columns <- setdiff(
  required_dmr_columns,
  names(selected_dmrs)
)

if (length(missing_dmr_columns) > 0L) {
  stop(
    "Selected DMR file is missing: ",
    paste(missing_dmr_columns, collapse = ", ")
  )
}

selected_dmrs <- selected_dmrs[splitID == split_id]
selected_dmrs[, region_id := as.integer(region_id)]
selected_dmrs[, data_chunk_id := as.integer(data_chunk_id)]
selected_dmrs[, chr := normalize_chr(chr)]
selected_dmrs[, region_start := as.integer(region_start)]
selected_dmrs[, region_end := as.integer(region_end)]
selected_dmrs <- unique(selected_dmrs, by = "region_id")
setorder(selected_dmrs, data_chunk_id, chr, region_start, region_end, region_id)

if (max_dmrs > 0L && nrow(selected_dmrs) > max_dmrs) {
  selected_dmrs <- head(selected_dmrs, max_dmrs)
}

cat("Strict MGCV DMRs selected in this run:", nrow(selected_dmrs), "\n")

train_split_dt <- fread(train_split_file)
train_cols <- grep(
  "^train_FID_[0-9]+$",
  names(train_split_dt),
  value = TRUE
)

if (!"splitID" %in% names(train_split_dt) || length(train_cols) == 0L) {
  stop("Training-split file lacks splitID or train_FID_* columns.")
}

split_row <- train_split_dt[splitID == split_id]
if (nrow(split_row) != 1L) {
  stop("Expected exactly one row for splitID ", split_id, ".")
}

train_fids <- unique(
  trimws(
    as.character(
      unlist(split_row[, ..train_cols], use.names = FALSE)
    )
  )
)
train_fids <- train_fids[
  !is.na(train_fids) & nzchar(train_fids) & train_fids != "NA"
]

pheno_env <- new.env(parent = emptyenv())
load(pheno_file_path, envir = pheno_env, verbose = TRUE)

if (!exists("pheno_file", envir = pheno_env, inherits = FALSE)) {
  stop("Phenotype RData does not contain pheno_file.")
}

pheno_dt <- as.data.table(
  copy(get("pheno_file", envir = pheno_env, inherits = FALSE))
)
rm(pheno_env)

required_pheno <- c(
  "ID",
  "FID",
  "AA_only",
  "AgeCalc",
  "Sex",
  "Non.smoker",
  "BMI"
)

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

train_dt <- pheno_dt[FID %chin% train_fids]
test_dt <- pheno_dt[!FID %chin% train_fids]

setorder(train_dt, ID)
setorder(test_dt, ID)

if (nrow(train_dt) == 0L || nrow(test_dt) == 0L) {
  stop("The outer split has no training or no test subjects.")
}

if (length(unique(train_dt$AA_only)) < 2L ||
    length(unique(test_dt$AA_only)) < 2L) {
  stop("Training and test sets must each contain both AA classes.")
}

if (length(intersect(train_dt$FID, test_dt$FID)) > 0L) {
  stop("Family leakage detected between the outer training and test sets.")
}

cat(
  "Outer training subjects:", nrow(train_dt),
  "families:", uniqueN(train_dt$FID),
  "cases:", sum(train_dt$AA_only == 1L),
  "controls:", sum(train_dt$AA_only == 0L),
  "\n"
)
cat(
  "Outer test subjects:", nrow(test_dt),
  "families:", uniqueN(test_dt$FID),
  "cases:", sum(test_dt$AA_only == 1L),
  "controls:", sum(test_dt$AA_only == 0L),
  "\n"
)

covariates <- make_covariate_matrices(train_dt, test_dt)

inner_fold_obj <- make_family_inner_folds(
  train_dt = train_dt,
  requested_k = inner_folds_requested,
  seed = seed_base + split_id,
  fold_path = fold_file
)

inner_foldid <- inner_fold_obj$foldid

if (length(inner_foldid) != nrow(train_dt)) {
  stop("Internal fold assignment length does not match training subjects.")
}

# ------------------------ split-specific DMR PC features ---------------------

build_dmr_features <- function(
  dmrs,
  train_ids,
  test_ids,
  chunk_directory
) {
  train_feature_list <- list()
  test_feature_list <- list()
  qc_list <- list()

  all_ids <- c(train_ids, test_ids)
  all_ids <- unique(all_ids)

  add_qc <- function(info, status, message = NA_character_,
                     n_cpg_raw = NA_integer_,
                     n_cpg_training_eligible = NA_integer_,
                     n_cpg_variable = NA_integer_,
                     n_pc = NA_integer_,
                     n_train_missing_before = NA_integer_,
                     n_test_missing_before = NA_integer_) {
    qc_list[[length(qc_list) + 1L]] <<- data.table(
      splitID = split_id,
      region_id = as.integer(info$region_id),
      data_chunk_id = as.integer(info$data_chunk_id),
      chr = as.character(info$chr),
      region_start = as.integer(info$region_start),
      region_end = as.integer(info$region_end),
      feature_status = as.character(status),
      feature_message = as.character(message),
      n_cpg_raw = as.integer(n_cpg_raw),
      n_cpg_training_eligible = as.integer(n_cpg_training_eligible),
      n_cpg_variable = as.integer(n_cpg_variable),
      n_pc = as.integer(n_pc),
      n_train_missing_before = as.integer(n_train_missing_before),
      n_test_missing_before = as.integer(n_test_missing_before)
    )
  }

  if (nrow(dmrs) == 0L) {
    return(
      list(
        train = matrix(
          numeric(0),
          nrow = length(train_ids),
          ncol = 0L,
          dimnames = list(train_ids, NULL)
        ),
        test = matrix(
          numeric(0),
          nrow = length(test_ids),
          ncol = 0L,
          dimnames = list(test_ids, NULL)
        ),
        qc = data.table()
      )
    )
  }

  dmr_chunks <- split(dmrs, dmrs$data_chunk_id)

  for (chunk_name in names(dmr_chunks)) {
    chunk_dmrs <- as.data.table(dmr_chunks[[chunk_name]])
    chunk_id <- as.integer(chunk_name)
    chunk_file <- file.path(
      chunk_directory,
      sprintf("chunk_%04d.csv", chunk_id)
    )

    if (!file.exists(chunk_file)) {
      for (i in seq_len(nrow(chunk_dmrs))) {
        add_qc(
          chunk_dmrs[i],
          status = "MISSING_CHUNK_FILE",
          message = chunk_file
        )
      }
      next
    }

    header <- tryCatch(
      fread(chunk_file, nrows = 0L, showProgress = FALSE),
      error = function(e) e
    )

    if (inherits(header, "error")) {
      for (i in seq_len(nrow(chunk_dmrs))) {
        add_qc(
          chunk_dmrs[i],
          status = "CHUNK_HEADER_ERROR",
          message = conditionMessage(header)
        )
      }
      next
    }

    header_names <- names(header)
    coord_cols <- intersect(c("chr", "start", "end"), header_names)

    if (!all(c("chr", "start") %in% coord_cols)) {
      for (i in seq_len(nrow(chunk_dmrs))) {
        add_qc(
          chunk_dmrs[i],
          status = "CHUNK_COORDINATE_ERROR",
          message = "Chunk lacks chr/start columns."
        )
      }
      next
    }

    meth_cols <- grep("_meth$", header_names, value = TRUE)
    meth_map <- data.table(
      source_column = meth_cols,
      ID = normalize_id(sub("_meth$", "", meth_cols))
    )

    meth_map <- meth_map[ID %chin% all_ids]

    if (anyDuplicated(meth_map$ID)) {
      duplicated_ids <- unique(meth_map$ID[duplicated(meth_map$ID)])
      stop(
        "Duplicated normalized methylation IDs in chunk ",
        chunk_id,
        ": ",
        paste(head(duplicated_ids, 10L), collapse = ", ")
      )
    }

    if (nrow(meth_map) == 0L) {
      for (i in seq_len(nrow(chunk_dmrs))) {
        add_qc(
          chunk_dmrs[i],
          status = "NO_MATCHING_METHYLATION_SAMPLES",
          message = "No *_meth columns matched phenotype IDs."
        )
      }
      next
    }

    selected_columns <- c(coord_cols, meth_map$source_column)

    chunk_dt <- tryCatch(
      fread(
        chunk_file,
        select = selected_columns,
        showProgress = FALSE
      ),
      error = function(e) e
    )

    if (inherits(chunk_dt, "error")) {
      for (i in seq_len(nrow(chunk_dmrs))) {
        add_qc(
          chunk_dmrs[i],
          status = "CHUNK_READ_ERROR",
          message = conditionMessage(chunk_dt)
        )
      }
      next
    }

    chunk_dt[, chr := normalize_chr(chr)]
    chunk_dt[, start := suppressWarnings(as.integer(start))]

    for (i in seq_len(nrow(chunk_dmrs))) {
      info <- chunk_dmrs[i]

      region_dt <- chunk_dt[
        chr == info$chr &
          start >= info$region_start &
          start <= info$region_end
      ]

      if (nrow(region_dt) == 0L) {
        add_qc(
          info,
          status = "NO_MATCHING_CPGS",
          message = "No chunk rows matched the DMR coordinates.",
          n_cpg_raw = 0L
        )
        next
      }

      setorder(region_dt, start)
      region_dt <- unique(region_dt, by = c("chr", "start"))

      source_cols <- meth_map$source_column
      meth_source <- as.matrix(
        region_dt[, ..source_cols]
      )
      storage.mode(meth_source) <- "numeric"
      colnames(meth_source) <- meth_map$ID

      meth_all <- matrix(
        NA_real_,
        nrow = nrow(meth_source),
        ncol = length(all_ids),
        dimnames = list(NULL, all_ids)
      )

      destination_index <- match(colnames(meth_source), all_ids)
      meth_all[, destination_index] <- meth_source

      meth_all[
        !is.finite(meth_all) |
          meth_all < 0 |
          meth_all > 1
      ] <- NA_real_

      meth_asr <- asin(sqrt(meth_all))

      train_matrix <- meth_asr[, train_ids, drop = FALSE]
      test_matrix <- meth_asr[, test_ids, drop = FALSE]

      n_train_missing_before <- sum(is.na(train_matrix))
      n_test_missing_before <- sum(is.na(test_matrix))

      observed_fraction <- rowMeans(is.finite(train_matrix))
      keep_observed <- is.finite(observed_fraction) &
        observed_fraction >= min_train_cpg_observed_fraction

      train_matrix <- train_matrix[keep_observed, , drop = FALSE]
      test_matrix <- test_matrix[keep_observed, , drop = FALSE]

      if (nrow(train_matrix) == 0L) {
        add_qc(
          info,
          status = "NO_TRAINING_ELIGIBLE_CPGS",
          message = paste0(
            "No CpGs had training observed fraction >= ",
            min_train_cpg_observed_fraction,
            "."
          ),
          n_cpg_raw = nrow(region_dt),
          n_cpg_training_eligible = 0L,
          n_cpg_variable = 0L,
          n_pc = 0L,
          n_train_missing_before = n_train_missing_before,
          n_test_missing_before = n_test_missing_before
        )
        next
      }

      train_medians <- apply(
        train_matrix,
        1L,
        function(z) {
          out <- stats::median(z[is.finite(z)], na.rm = TRUE)
          if (is.finite(out)) out else NA_real_
        }
      )

      keep_median <- is.finite(train_medians)
      train_matrix <- train_matrix[keep_median, , drop = FALSE]
      test_matrix <- test_matrix[keep_median, , drop = FALSE]
      train_medians <- train_medians[keep_median]

      if (nrow(train_matrix) == 0L) {
        add_qc(
          info,
          status = "NO_FINITE_TRAINING_MEDIANS",
          n_cpg_raw = nrow(region_dt),
          n_cpg_training_eligible = 0L,
          n_cpg_variable = 0L,
          n_pc = 0L,
          n_train_missing_before = n_train_missing_before,
          n_test_missing_before = n_test_missing_before
        )
        next
      }

      for (r in seq_len(nrow(train_matrix))) {
        train_matrix[r, !is.finite(train_matrix[r, ])] <- train_medians[r]
        test_matrix[r, !is.finite(test_matrix[r, ])] <- train_medians[r]
      }

      training_sd <- apply(train_matrix, 1L, stats::sd)
      keep_variable <- is.finite(training_sd) &
        training_sd > variance_epsilon

      train_matrix <- train_matrix[keep_variable, , drop = FALSE]
      test_matrix <- test_matrix[keep_variable, , drop = FALSE]

      n_variable <- nrow(train_matrix)

      if (n_variable < min_variable_cpgs) {
        add_qc(
          info,
          status = "TOO_FEW_VARIABLE_CPGS",
          message = paste0(
            "Variable training CpGs = ",
            n_variable,
            "; required = ",
            min_variable_cpgs,
            "."
          ),
          n_cpg_raw = nrow(region_dt),
          n_cpg_training_eligible = length(train_medians),
          n_cpg_variable = n_variable,
          n_pc = 0L,
          n_train_missing_before = n_train_missing_before,
          n_test_missing_before = n_test_missing_before
        )
        next
      }

      n_pc <- min(
        pca_components_requested,
        n_variable,
        ncol(train_matrix) - 1L
      )

      if (n_pc < 1L) {
        add_qc(
          info,
          status = "PCA_DIMENSION_ERROR",
          n_cpg_raw = nrow(region_dt),
          n_cpg_training_eligible = length(train_medians),
          n_cpg_variable = n_variable,
          n_pc = 0L,
          n_train_missing_before = n_train_missing_before,
          n_test_missing_before = n_test_missing_before
        )
        next
      }

      pca_fit <- tryCatch(
        stats::prcomp(
          x = t(train_matrix),
          center = TRUE,
          scale. = TRUE,
          rank. = n_pc
        ),
        error = function(e) e
      )

      if (inherits(pca_fit, "error")) {
        add_qc(
          info,
          status = "PCA_ERROR",
          message = conditionMessage(pca_fit),
          n_cpg_raw = nrow(region_dt),
          n_cpg_training_eligible = length(train_medians),
          n_cpg_variable = n_variable,
          n_pc = 0L,
          n_train_missing_before = n_train_missing_before,
          n_test_missing_before = n_test_missing_before
        )
        next
      }

      n_pc_available <- min(
        n_pc,
        ncol(pca_fit$x),
        ncol(pca_fit$rotation)
      )

      if (!is.finite(n_pc_available) || n_pc_available < 1L) {
        add_qc(
          info,
          status = "PCA_NO_COMPONENTS",
          n_cpg_raw = nrow(region_dt),
          n_cpg_training_eligible = length(train_medians),
          n_cpg_variable = n_variable,
          n_pc = 0L,
          n_train_missing_before = n_train_missing_before,
          n_test_missing_before = n_test_missing_before
        )
        next
      }

      n_pc <- as.integer(n_pc_available)
      train_scores <- pca_fit$x[, seq_len(n_pc), drop = FALSE]
      test_scores <- tryCatch(
        stats::predict(
          pca_fit,
          newdata = t(test_matrix)
        )[, seq_len(n_pc), drop = FALSE],
        error = function(e) e
      )

      if (inherits(test_scores, "error")) {
        add_qc(
          info,
          status = "PCA_PROJECTION_ERROR",
          message = conditionMessage(test_scores),
          n_cpg_raw = nrow(region_dt),
          n_cpg_training_eligible = length(train_medians),
          n_cpg_variable = n_variable,
          n_pc = 0L,
          n_train_missing_before = n_train_missing_before,
          n_test_missing_before = n_test_missing_before
        )
        next
      }

      feature_names <- paste0(
        "DMRpc",
        seq_len(n_pc),
        "_",
        info$region_id
      )

      colnames(train_scores) <- feature_names
      colnames(test_scores) <- feature_names
      rownames(train_scores) <- train_ids
      rownames(test_scores) <- test_ids

      train_feature_list[[as.character(info$region_id)]] <- train_scores
      test_feature_list[[as.character(info$region_id)]] <- test_scores

      add_qc(
        info,
        status = "FEATURE_OK",
        n_cpg_raw = nrow(region_dt),
        n_cpg_training_eligible = length(train_medians),
        n_cpg_variable = n_variable,
        n_pc = n_pc,
        n_train_missing_before = n_train_missing_before,
        n_test_missing_before = n_test_missing_before
      )
    }

    rm(chunk_dt, header)
    invisible(gc(FALSE))
  }

  if (length(train_feature_list) > 0L) {
    train_features <- do.call(cbind, train_feature_list)
    test_features <- do.call(cbind, test_feature_list)
  } else {
    train_features <- matrix(
      numeric(0),
      nrow = length(train_ids),
      ncol = 0L,
      dimnames = list(train_ids, NULL)
    )

    test_features <- matrix(
      numeric(0),
      nrow = length(test_ids),
      ncol = 0L,
      dimnames = list(test_ids, NULL)
    )
  }

  if (length(qc_list) > 0L) {
    qc_dt <- rbindlist(qc_list, use.names = TRUE, fill = TRUE)
  } else {
    qc_dt <- data.table()
  }

  list(
    train = train_features,
    test = test_features,
    qc = qc_dt
  )
}

feature_obj <- build_dmr_features(
  dmrs = selected_dmrs,
  train_ids = train_dt$ID,
  test_ids = test_dt$ID,
  chunk_directory = chunk_dir
)

if (nrow(feature_obj$qc) > 0L) {
  feature_qc_out <- merge(
    feature_obj$qc,
    selected_dmrs,
    by = c(
      "splitID",
      "region_id",
      "data_chunk_id",
      "chr",
      "region_start",
      "region_end"
    ),
    all.x = TRUE,
    sort = FALSE
  )
} else {
  feature_qc_out <- data.table(
    splitID = integer(),
    region_id = integer(),
    data_chunk_id = integer(),
    chr = character(),
    region_start = integer(),
    region_end = integer(),
    feature_status = character(),
    feature_message = character(),
    n_cpg_raw = integer(),
    n_cpg_training_eligible = integer(),
    n_cpg_variable = integer(),
    n_pc = integer(),
    n_train_missing_before = integer(),
    n_test_missing_before = integer()
  )
}

fwrite(
  feature_qc_out,
  feature_qc_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

x_dmr_train <- feature_obj$train
x_dmr_test <- feature_obj$test

if (!identical(rownames(x_dmr_train), train_dt$ID) ||
    !identical(rownames(x_dmr_test), test_dt$ID)) {
  stop("DMR feature sample order does not match phenotype sample order.")
}

cat(
  "Usable MGCV DMRs:",
  if (nrow(feature_obj$qc) > 0L) {
    feature_obj$qc[feature_status == "FEATURE_OK", uniqueN(region_id)]
  } else {
    0L
  },
  "methylation features:", ncol(x_dmr_train),
  "\n"
)

# --------------------------- prediction model helpers ------------------------

fit_covariate_reference <- function(
  x_train,
  y_train,
  x_test,
  foldid
) {
  train_df <- as.data.frame(x_train)
  train_df$y <- as.integer(y_train)
  test_df <- as.data.frame(x_test)

  oof_prob <- rep(NA_real_, length(y_train))

  for (fold in sort(unique(foldid))) {
    fit_idx <- foldid != fold
    val_idx <- foldid == fold

    fit <- tryCatch(
      suppressWarnings(
        stats::glm(
          y ~ .,
          data = train_df[fit_idx, , drop = FALSE],
          family = stats::binomial()
        )
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      stop(
        "Covariate-only inner fold failed: ",
        conditionMessage(fit)
      )
    }

    oof_prob[val_idx] <- as.numeric(
      stats::predict(
        fit,
        newdata = train_df[val_idx, , drop = FALSE],
        type = "response"
      )
    )
  }

  final_fit <- tryCatch(
    suppressWarnings(
      stats::glm(
        y ~ .,
        data = train_df,
        family = stats::binomial()
      )
    ),
    error = function(e) e
  )

  if (inherits(final_fit, "error")) {
    stop("Final covariate-only model failed: ", conditionMessage(final_fit))
  }

  test_prob <- as.numeric(
    stats::predict(
      final_fit,
      newdata = test_df,
      type = "response"
    )
  )

  threshold <- choose_youden_threshold(y_train, oof_prob)

  coef_dt <- data.table(
    feature = names(stats::coef(final_fit)),
    coefficient = as.numeric(stats::coef(final_fit))
  )

  coef_dt[, feature_type := fifelse(
    feature == "(Intercept)",
    "intercept",
    "covariate"
  )]

  list(
    oof_prob = oof_prob,
    test_prob = test_prob,
    threshold = threshold,
    coefficients = coef_dt,
    lambda = NA_real_,
    n_features_input = ncol(x_train),
    n_features_nonzero = sum(
      coef_dt$feature != "(Intercept)" &
        is.finite(coef_dt$coefficient) &
        coef_dt$coefficient != 0
    ),
    n_methylation_nonzero = 0L,
    fit_status = "FIT_OK"
  )
}

fit_mgcv_lasso <- function(
  x_dmr_train,
  x_dmr_test,
  x_cov_train,
  x_cov_test,
  y_train,
  foldid,
  seed
) {
  if (ncol(x_dmr_train) == 0L) {
    return(NULL)
  }

  x_train <- cbind(x_dmr_train, x_cov_train)
  x_test <- cbind(x_dmr_test, x_cov_test)

  training_sd <- apply(x_train, 2L, stats::sd)
  keep <- is.finite(training_sd) & training_sd > variance_epsilon

  # Preserve all nonconstant predictors; a constant unpenalized covariate adds
  # no information and is safely removed.
  x_train <- x_train[, keep, drop = FALSE]
  x_test <- x_test[, keep, drop = FALSE]

  dmr_columns <- grepl("^DMRpc[0-9]+_", colnames(x_train))
  covariate_columns <- grepl("^COV_", colnames(x_train))

  if (!any(dmr_columns)) {
    return(NULL)
  }

  penalty_factor <- ifelse(dmr_columns, 1, 0)

  set.seed(seed)

  cvfit <- tryCatch(
    glmnet::cv.glmnet(
      x = x_train,
      y = as.integer(y_train),
      family = "binomial",
      alpha = 1,
      foldid = as.integer(foldid),
      type.measure = "auc",
      standardize = TRUE,
      intercept = TRUE,
      penalty.factor = penalty_factor,
      keep = TRUE,
      grouped = TRUE,
      maxit = 1000000,
      parallel = FALSE
    ),
    error = function(e) e
  )

  if (inherits(cvfit, "error")) {
    stop("MGCV LASSO failed: ", conditionMessage(cvfit))
  }


  lambda_selected <- switch(
  lambda_rule,
  "lambda.1se" = as.numeric(cvfit$lambda.1se),
  "lambda.min" = as.numeric(cvfit$lambda.min)
   )


  lambda_index <- which.min(abs(cvfit$lambda - lambda_selected))

  oof_prob <- as.numeric(cvfit$fit.preval[, lambda_index])

  # cv.glmnet normally stores response-scale predictions for AUC. This guard
  # also handles versions that store link-scale predictions.
  if (any(oof_prob < 0 | oof_prob > 1, na.rm = TRUE)) {
    oof_prob <- stats::plogis(oof_prob)
  }

  test_prob <- as.numeric(
    stats::predict(
      cvfit,
      newx = x_test,
      s = lambda_rule,
      type = "response"
    )
  )

  threshold <- choose_youden_threshold(y_train, oof_prob)

  coef_matrix <- stats::coef(cvfit, s = lambda_rule)
  coef_dt <- data.table(
    feature = rownames(coef_matrix),
    coefficient = as.numeric(coef_matrix)
  )

  coef_dt[, feature_type := fcase(
    feature == "(Intercept)", "intercept",
    grepl("^DMRpc[0-9]+_", feature), "methylation",
    grepl("^COV_", feature), "covariate",
    default = "other"
  )]

  nonzero <- coef_dt[
    feature != "(Intercept)" &
      is.finite(coefficient) &
      coefficient != 0
  ]

  list(
    oof_prob = oof_prob,
    test_prob = test_prob,
    threshold = threshold,
    coefficients = coef_dt,
    lambda = lambda_selected,
    n_features_input = ncol(x_train),
    n_features_nonzero = nrow(nonzero),
    n_methylation_nonzero = nonzero[feature_type == "methylation", .N],
    fit_status = "FIT_OK"
  )
}

# ------------------------------- fit models ----------------------------------

y_train <- train_dt$AA_only
y_test <- test_dt$AA_only

baseline_fit <- fit_covariate_reference(
  x_train = covariates$train,
  y_train = y_train,
  x_test = covariates$test,
  foldid = inner_foldid
)

mgcv_fit <- fit_mgcv_lasso(
  x_dmr_train = x_dmr_train,
  x_dmr_test = x_dmr_test,
  x_cov_train = covariates$train,
  x_cov_test = covariates$test,
  y_train = y_train,
  foldid = inner_foldid,
  seed = seed_base + 100000L + split_id
)

mgcv_fallback <- FALSE

if (is.null(mgcv_fit)) {
  mgcv_fallback <- TRUE
  mgcv_fit <- baseline_fit
  mgcv_fit$fit_status <- "COVARIATE_ONLY_FALLBACK_NO_USABLE_DMR_FEATURES"
}

# ----------------------------- collect outputs -------------------------------

make_model_outputs <- function(
  fit,
  model_name,
  model_type,
  fallback,
  n_dmr_selected,
  n_dmr_usable,
  n_methylation_features
) {
  inner_metrics <- calculate_metrics(
    y = y_train,
    p = fit$oof_prob,
    threshold = fit$threshold
  )

  test_metrics <- calculate_metrics(
    y = y_test,
    p = fit$test_prob,
    threshold = fit$threshold
  )

  performance <- data.table(
    splitID = split_id,
    DMR_method = "MGCV",
    model = model_name,
    model_type = model_type,
    alpha = if (model_type == "LASSO") 1 else NA_real_,
    lambda = fit$lambda,
    lambda_rule = if(model_type == "LASSO") lambda_rule else NA_character_, 
    threshold = fit$threshold,
    threshold_source = "family_grouped_inner_CV_Youden",
    inner_folds = inner_fold_obj$n_folds,
    fit_status = fit$fit_status,
    model_fallback = fallback,
    n_train = nrow(train_dt),
    n_test = nrow(test_dt),
    n_train_FIDs = uniqueN(train_dt$FID),
    n_test_FIDs = uniqueN(test_dt$FID),
    n_DMR_selected_input = n_dmr_selected,
    n_DMR_usable = n_dmr_usable,
    n_methylation_features_input = n_methylation_features,
    n_features_model_input = fit$n_features_input,
    n_features_nonzero = fit$n_features_nonzero,
    n_methylation_features_nonzero = fit$n_methylation_nonzero
  )

  performance <- cbind(
    performance,
    as.data.table(metrics_to_columns(inner_metrics, "inner_oof")),
    as.data.table(metrics_to_columns(test_metrics, "test"))
  )

  prediction_train <- data.table(
    splitID = split_id,
    DMR_method = "MGCV",
    model = model_name,
    ID = train_dt$ID,
    FID = train_dt$FID,
    AA_only = y_train,
    prediction_set = "train_inner_oof",
    inner_fold = inner_foldid,
    pred_prob = fit$oof_prob,
    threshold = fit$threshold,
    pred_class = as.integer(fit$oof_prob >= fit$threshold)
  )

  prediction_test <- data.table(
    splitID = split_id,
    DMR_method = "MGCV",
    model = model_name,
    ID = test_dt$ID,
    FID = test_dt$FID,
    AA_only = y_test,
    prediction_set = "test",
    inner_fold = NA_integer_,
    pred_prob = fit$test_prob,
    threshold = fit$threshold,
    pred_class = as.integer(fit$test_prob >= fit$threshold)
  )

  coefficients <- copy(fit$coefficients)
  coefficients[, `:=`(
    splitID = split_id,
    DMR_method = "MGCV",
    model = model_name,
    selected_nonzero = as.integer(
      is.finite(coefficient) & coefficient != 0
    )
  )]

  list(
    performance = performance,
    predictions = rbindlist(
      list(prediction_train, prediction_test),
      use.names = TRUE
    ),
    coefficients = coefficients
  )
}

n_dmr_usable <- if (nrow(feature_obj$qc) > 0L) {
  feature_obj$qc[feature_status == "FEATURE_OK", uniqueN(region_id)]
} else {
  0L
}

mgcv_output <- make_model_outputs(
  fit = mgcv_fit,
  model_name = "MGCV_PC12_LASSO",
  model_type = "LASSO",
  fallback = mgcv_fallback,
  n_dmr_selected = nrow(selected_dmrs),
  n_dmr_usable = n_dmr_usable,
  n_methylation_features = ncol(x_dmr_train)
)

baseline_output <- make_model_outputs(
  fit = baseline_fit,
  model_name = "Covariate_only",
  model_type = "unpenalized_logistic",
  fallback = FALSE,
  n_dmr_selected = 0L,
  n_dmr_usable = 0L,
  n_methylation_features = 0L
)

performance_out <- rbindlist(
  list(mgcv_output$performance, baseline_output$performance),
  use.names = TRUE,
  fill = TRUE
)

prediction_out <- rbindlist(
  list(mgcv_output$predictions, baseline_output$predictions),
  use.names = TRUE,
  fill = TRUE
)

coefficient_out <- rbindlist(
  list(mgcv_output$coefficients, baseline_output$coefficients),
  use.names = TRUE,
  fill = TRUE
)

setcolorder(
  coefficient_out,
  c(
    "splitID",
    "DMR_method",
    "model",
    "feature",
    "feature_type",
    "coefficient",
    "selected_nonzero"
  )
)

fwrite(
  performance_out,
  performance_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  prediction_out,
  prediction_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  coefficient_out,
  coefficient_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

split_qc <- data.table(
  splitID = split_id,
  DMR_method = "MGCV",
  selected_dmr_file = selected_dmr_file,
  n_selected_DMRs = nrow(selected_dmrs),
  n_usable_DMRs = n_dmr_usable,
  n_DMR_PC_features = ncol(x_dmr_train),
  n_feature_failures = if (nrow(feature_obj$qc) > 0L) {
    feature_obj$qc[feature_status != "FEATURE_OK", .N]
  } else {
    0L
  },
  n_train = nrow(train_dt),
  n_test = nrow(test_dt),
  n_train_FIDs = uniqueN(train_dt$FID),
  n_test_FIDs = uniqueN(test_dt$FID),
  n_train_cases = sum(y_train == 1L),
  n_train_controls = sum(y_train == 0L),
  n_test_cases = sum(y_test == 1L),
  n_test_controls = sum(y_test == 0L),
  inner_folds = inner_fold_obj$n_folds,
  family_overlap_count = length(intersect(train_dt$FID, test_dt$FID)),
  mgcv_model_fallback = mgcv_fallback,
  elapsed_minutes = as.numeric(
    difftime(Sys.time(), run_start, units = "mins")
  ),
  completed_at = format(
    Sys.time(),
    tz = "America/Toronto",
    usetz = TRUE
  )
)

fwrite(
  split_qc,
  split_qc_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\nFinished MGCV primary predictive analysis for split", split_id, "\n")
cat("Performance:", performance_file, "\n")
cat("Predictions:", prediction_file, "\n")
cat("Coefficients:", coefficient_file, "\n")
cat("Feature QC:", feature_qc_file, "\n")
cat("Split QC:", split_qc_file, "\n")
print(performance_out[, .(
  model,
  fit_status,
  n_DMR_selected_input,
  n_DMR_usable,
  n_methylation_features_input,
  AUROC_inner_oof,
  AUROC_test,
  PR_AUC_test,
  Brier_test,
  calibration_slope_test
)])

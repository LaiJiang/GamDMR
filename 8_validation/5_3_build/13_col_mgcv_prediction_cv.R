#!/usr/bin/env Rscript

# Collect split-specific MGCV prediction results from the intentional
# all-subject PCA sensitivity analysis.
#
# Expected upstream analysis:
#   10_run_mgcv_test_pca_all_subjects.R
#
# In that sensitivity analysis, PCA is fitted using the combined outer-training
# and outer-test subjects before the resulting scores are separated. This is an
# intentional leakage experiment and must not be reported as valid external-CV
# performance.
#
# Outputs include:
#   - combined split-level result files;
#   - descriptive performance summaries across splits;
#   - paired MGCV-versus-covariate-only comparisons;
#   - corrected repeated-CV confidence intervals;
#   - subject-averaged repeated test predictions and ROC/PR summaries;
#   - optional paired comparison with the proper training-only PCA analysis,
#     provided the two analyses use the same lambda rule.

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(PRROC)
  library(ggplot2)
})

data.table::setDTthreads(1L)
options(mc.cores = 1L)

# ------------------------------- configuration -------------------------------

PATH_wk <- path.expand(
  Sys.getenv("PATH_WK", "~/scratch/UQAC/meth/")
)

input_root <- Sys.getenv(
  "MGCV_PREDICTION_OUTPUT",
  file.path(
    PATH_wk,
    "results/15_revision/5_cv/5_prediction/1_mgcv_pca_all_subjects_test"
  )
)

performance_dir <- file.path(input_root, "performance")
prediction_dir <- file.path(input_root, "predictions")
coefficient_dir <- file.path(input_root, "coefficients")
feature_qc_dir <- file.path(input_root, "feature_qc")
split_qc_dir <- file.path(input_root, "split_qc")

output_dir <- file.path(input_root, "collected")
figure_dir <- file.path(output_dir, "figures")
qc_dir <- file.path(output_dir, "qc")

for (d in c(output_dir, figure_dir, qc_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

expected_splits <- suppressWarnings(
  as.integer(Sys.getenv("EXPECTED_SPLITS", "100"))
)

primary_model_name <- Sys.getenv(
  "PRIMARY_MODEL_NAME",
  "MGCV_PC12_LASSO"
)

expected_pca_scope <- Sys.getenv(
  "EXPECTED_PCA_SCOPE",
  "all_outer_split_subjects_train_plus_test_intentional_leakage"
)

compare_with_training_only <- Sys.getenv(
  "COMPARE_WITH_TRAINING_ONLY",
  "1"
) %in% c("1", "TRUE", "true", "T", "yes", "YES")

training_only_performance_file <- Sys.getenv(
  "TRAINING_ONLY_PERFORMANCE_FILE",
  file.path(
    PATH_wk,
    paste0(
      "results/15_revision/5_cv/5_prediction/1_mgcv_primary/",
      "collected/mgcv_primary_performance_all_splits.tsv"
    )
  )
)

allow_lambda_mismatch <- Sys.getenv(
  "ALLOW_LAMBDA_MISMATCH",
  "0"
) %in% c("1", "TRUE", "true", "T", "yes", "YES")

overwrite <- Sys.getenv("OVERWRITE", "0") %in%
  c("1", "TRUE", "true", "T", "yes", "YES")

if (is.na(expected_splits) || expected_splits < 1L) {
  stop("EXPECTED_SPLITS must be a positive integer.")
}

cat("Input root:", input_root, "\n")
cat("Expected splits:", expected_splits, "\n")
cat("Primary model:", primary_model_name, "\n")
cat("Expected PCA scope:", expected_pca_scope, "\n")

# -------------------------------- utilities ----------------------------------

safe_roc_auc <- function(y, p) {
  y <- as.integer(y)
  p <- as.numeric(p)
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]
  p <- p[keep]

  if (length(unique(y)) < 2L) return(NA_real_)

  obj <- tryCatch(
    pROC::roc(
      response = y,
      predictor = p,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(obj)) NA_real_ else as.numeric(pROC::auc(obj))
}

safe_pr_auc <- function(y, p) {
  y <- as.integer(y)
  p <- as.numeric(p)
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]
  p <- p[keep]

  if (length(unique(y)) < 2L) return(NA_real_)

  obj <- tryCatch(
    PRROC::pr.curve(
      scores.class0 = p[y == 1L],
      scores.class1 = p[y == 0L],
      curve = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(obj)) NA_real_ else as.numeric(obj$auc.integral)
}

aggregate_probability_metrics <- function(y, p) {
  y <- as.integer(y)
  p <- as.numeric(p)
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]
  p <- p[keep]

  if (length(y) == 0L) {
    return(
      data.table(
        n_subjects = 0L,
        n_cases = 0L,
        n_controls = 0L,
        prevalence = NA_real_,
        AUROC = NA_real_,
        PR_AUC = NA_real_,
        PR_baseline = NA_real_,
        PR_AUC_lift = NA_real_,
        Brier = NA_real_,
        log_loss = NA_real_,
        calibration_intercept = NA_real_,
        calibration_slope = NA_real_
      )
    )
  }

  p_clip <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  lp <- qlogis(p_clip)

  calibration_intercept <- NA_real_
  calibration_slope <- NA_real_
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

  pr_auc <- safe_pr_auc(y, p)

  data.table(
    n_subjects = length(y),
    n_cases = sum(y == 1L),
    n_controls = sum(y == 0L),
    prevalence = mean(y),
    AUROC = safe_roc_auc(y, p),
    PR_AUC = pr_auc,
    PR_baseline = mean(y),
    PR_AUC_lift = pr_auc - mean(y),
    Brier = mean((p - y)^2),
    log_loss = -mean(y * log(p_clip) + (1 - y) * log(1 - p_clip)),
    calibration_intercept = calibration_intercept,
    calibration_slope = calibration_slope
  )
}

summarize_numeric <- function(x, prefix) {
  x <- x[is.finite(x)]

  statistic_names <- c(
    "mean", "sd", "median", "q025", "q25", "q75", "q975", "min", "max"
  )

  if (length(x) == 0L) {
    out <- rep(NA_real_, length(statistic_names))
    names(out) <- statistic_names
  } else {
    out <- c(
      mean = mean(x),
      sd = stats::sd(x),
      median = stats::median(x),
      q025 = stats::quantile(x, 0.025, names = FALSE),
      q25 = stats::quantile(x, 0.25, names = FALSE),
      q75 = stats::quantile(x, 0.75, names = FALSE),
      q975 = stats::quantile(x, 0.975, names = FALSE),
      min = min(x),
      max = max(x)
    )
  }

  names(out) <- paste0(prefix, "_", names(out))
  as.list(out)
}

corrected_repeated_cv <- function(difference, test_train_ratio) {
  keep <- is.finite(difference) & is.finite(test_train_ratio)
  d <- difference[keep]
  ratio <- test_train_ratio[keep]
  n <- length(d)

  if (n < 3L) {
    return(
      data.table(
        n_splits = n,
        mean_difference = if (n > 0L) mean(d) else NA_real_,
        median_difference = if (n > 0L) median(d) else NA_real_,
        corrected_SE = NA_real_,
        corrected_CI_low = NA_real_,
        corrected_CI_high = NA_real_,
        corrected_t = NA_real_,
        corrected_df = NA_integer_,
        corrected_pvalue = NA_real_,
        mean_test_train_ratio = if (n > 0L) mean(ratio) else NA_real_,
        raw_positive_rate = if (n > 0L) mean(d > 0) else NA_real_
      )
    )
  }

  mean_d <- mean(d)
  variance_d <- stats::var(d)
  ratio_mean <- mean(ratio)
  corrected_se <- sqrt((1 / n + ratio_mean) * variance_d)
  df <- n - 1L

  if (!is.finite(corrected_se) || corrected_se == 0) {
    t_value <- NA_real_
    p_value <- NA_real_
    ci_low <- mean_d
    ci_high <- mean_d
  } else {
    t_value <- mean_d / corrected_se
    p_value <- 2 * stats::pt(-abs(t_value), df = df)
    critical <- stats::qt(0.975, df = df)
    ci_low <- mean_d - critical * corrected_se
    ci_high <- mean_d + critical * corrected_se
  }

  data.table(
    n_splits = n,
    mean_difference = mean_d,
    median_difference = median(d),
    corrected_SE = corrected_se,
    corrected_CI_low = ci_low,
    corrected_CI_high = ci_high,
    corrected_t = t_value,
    corrected_df = df,
    corrected_pvalue = p_value,
    mean_test_train_ratio = ratio_mean,
    raw_positive_rate = mean(d > 0)
  )
}

repair_duplicate_column_names <- function(dt, file_path, label) {
  current_names <- names(dt)

  if (!anyDuplicated(current_names)) return(dt)

  n_test_index <- which(current_names == "n_test")

  if (length(n_test_index) > 1L) {
    replacement_names <- paste0(
      "n_evaluable_test",
      ifelse(seq_along(n_test_index[-1L]) == 1L, "", paste0("_", seq_along(n_test_index[-1L])))
    )
    current_names[n_test_index[-1L]] <- replacement_names
  }

  if (anyDuplicated(current_names)) {
    duplicated_names <- unique(
      current_names[
        duplicated(current_names) |
          duplicated(current_names, fromLast = TRUE)
      ]
    )

    stop(
      label,
      " file still has duplicated columns after repair: ",
      basename(file_path),
      " [",
      paste(duplicated_names, collapse = ", "),
      "]"
    )
  }

  setnames(dt, current_names)
  dt
}

read_one_table <- function(file_path, label) {
  if (!file.exists(file_path)) {
    stop(label, " file does not exist: ", file_path)
  }

  if (file.info(file_path)$size == 0L) {
    stop(label, " file is empty: ", file_path)
  }

  dt <- fread(file_path, check.names = FALSE)
  repair_duplicate_column_names(dt, file_path, label)
}

read_and_bind <- function(files, label) {
  if (length(files) == 0L) {
    stop("No ", label, " files were found.")
  }

  table_list <- lapply(
    files,
    function(file_path) read_one_table(file_path, label)
  )

  out <- rbindlist(
    table_list,
    use.names = TRUE,
    fill = TRUE
  )

  rm(table_list)
  invisible(gc(FALSE))

  cat(label, "files:", length(files), "rows:", nrow(out), "\n")
  out
}

make_corrected_table <- function(
  paired,
  metrics,
  difference_prefix,
  comparison_label
) {
  rbindlist(
    lapply(
      metrics,
      function(metric) {
        difference_column <- paste0(difference_prefix, metric)
        d_values <- paired[[difference_column]]
        lower_is_better <- metric %in% c("Brier_test", "log_loss_test")

        result <- corrected_repeated_cv(
          difference = d_values,
          test_train_ratio = paired$test_train_ratio
        )

        result[, `:=`(
          comparison = comparison_label,
          metric = metric,
          difference_definition = difference_column,
          favorable_direction = if (lower_is_better) "negative" else "positive",
          win_rate = if (lower_is_better) {
            mean(d_values < 0, na.rm = TRUE)
          } else {
            mean(d_values > 0, na.rm = TRUE)
          }
        )]

        setcolorder(
          result,
          c(
            "comparison",
            "metric",
            "difference_definition",
            "favorable_direction",
            setdiff(
              names(result),
              c(
                "comparison",
                "metric",
                "difference_definition",
                "favorable_direction"
              )
            )
          )
        )

        result
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )
}

# -------------------------------- input files --------------------------------

performance_files <- sort(list.files(
  performance_dir,
  pattern = "^mgcv_primary_performance_split_[0-9]{3}\\.tsv$",
  full.names = TRUE
))

prediction_files <- sort(list.files(
  prediction_dir,
  pattern = "^mgcv_primary_predictions_split_[0-9]{3}\\.tsv$",
  full.names = TRUE
))

coefficient_files <- sort(list.files(
  coefficient_dir,
  pattern = "^mgcv_primary_coefficients_split_[0-9]{3}\\.tsv$",
  full.names = TRUE
))

feature_qc_files <- sort(list.files(
  feature_qc_dir,
  pattern = "^mgcv_primary_feature_qc_split_[0-9]{3}\\.tsv$",
  full.names = TRUE
))

split_qc_files <- sort(list.files(
  split_qc_dir,
  pattern = "^mgcv_primary_split_qc_[0-9]{3}\\.tsv$",
  full.names = TRUE
))

file_count_qc <- data.table(
  file_type = c(
    "performance", "prediction", "coefficient", "feature_qc", "split_qc"
  ),
  n_files = c(
    length(performance_files),
    length(prediction_files),
    length(coefficient_files),
    length(feature_qc_files),
    length(split_qc_files)
  ),
  expected_files = expected_splits
)

fwrite(
  file_count_qc,
  file.path(qc_dir, "mgcv_pca_all_subjects_file_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

print(file_count_qc)

if (length(performance_files) != expected_splits ||
    length(prediction_files) != expected_splits ||
    length(split_qc_files) != expected_splits) {
  stop(
    "Expected ", expected_splits,
    " performance, prediction, and split-QC files before collection."
  )
}

performance_all <- read_and_bind(performance_files, "Performance")
prediction_all <- read_and_bind(prediction_files, "Prediction")
coefficient_all <- read_and_bind(coefficient_files, "Coefficient")
feature_qc_all <- read_and_bind(feature_qc_files, "Feature QC")
split_qc_all <- read_and_bind(split_qc_files, "Split QC")

required_performance_columns <- c(
  "splitID", "model", "n_train", "n_test", "AUROC_test", "PR_AUC_test",
  "Brier_test", "log_loss_test"
)
missing_performance_columns <- setdiff(
  required_performance_columns,
  names(performance_all)
)

if (length(missing_performance_columns) > 0L) {
  stop(
    "Collected performance results lack: ",
    paste(missing_performance_columns, collapse = ", ")
  )
}

if (uniqueN(performance_all$splitID) != expected_splits ||
    uniqueN(prediction_all$splitID) != expected_splits) {
  stop("Not all expected splitIDs are represented in the collected results.")
}

expected_split_ids <- seq_len(expected_splits)
missing_performance_splits <- setdiff(
  expected_split_ids,
  sort(unique(performance_all$splitID))
)
missing_prediction_splits <- setdiff(
  expected_split_ids,
  sort(unique(prediction_all$splitID))
)

if (length(missing_performance_splits) > 0L ||
    length(missing_prediction_splits) > 0L) {
  stop(
    "Missing split IDs. Performance: ",
    paste(missing_performance_splits, collapse = ","),
    "; predictions: ",
    paste(missing_prediction_splits, collapse = ",")
  )
}

if (!"pca_scope" %in% names(performance_all)) {
  warning("pca_scope is absent from performance files; adding expected label.")
  performance_all[, pca_scope := expected_pca_scope]
}

if (!"pca_scope" %in% names(split_qc_all)) {
  split_qc_all[, pca_scope := expected_pca_scope]
}

observed_pca_scopes <- unique(
  performance_all[model == primary_model_name, as.character(pca_scope)]
)
observed_pca_scopes <- observed_pca_scopes[
  !is.na(observed_pca_scopes) & nzchar(observed_pca_scopes)
]

if (length(observed_pca_scopes) != 1L ||
    observed_pca_scopes != expected_pca_scope) {
  stop(
    "Unexpected PCA scope in primary model: ",
    paste(observed_pca_scopes, collapse = ", "),
    ". Expected: ",
    expected_pca_scope
  )
}

primary_rows <- performance_all[model == primary_model_name]
baseline_rows <- performance_all[model == "Covariate_only"]

if (nrow(primary_rows) != expected_splits) {
  stop(
    "Expected one ", primary_model_name,
    " row per split; observed ", nrow(primary_rows), "."
  )
}

if (nrow(baseline_rows) != expected_splits) {
  stop(
    "Expected one Covariate_only row per split; observed ",
    nrow(baseline_rows), "."
  )
}

if (anyDuplicated(primary_rows$splitID) || anyDuplicated(baseline_rows$splitID)) {
  stop("Duplicated primary or covariate-only performance rows by splitID.")
}

observed_lambda_rules <- if ("lambda_rule" %in% names(primary_rows)) {
  unique(na.omit(as.character(primary_rows$lambda_rule)))
} else {
  character()
}

if (length(observed_lambda_rules) > 1L) {
  stop(
    "More than one lambda rule is present across primary-model splits: ",
    paste(observed_lambda_rules, collapse = ", ")
  )
}

analysis_lambda_rule <- if (length(observed_lambda_rules) == 1L) {
  observed_lambda_rules
} else {
  NA_character_
}

cat("Collected PCA scope:", observed_pca_scopes, "\n")
cat("Collected lambda rule:", analysis_lambda_rule, "\n")

# -------------------------------- output files -------------------------------

combined_performance_file <- file.path(
  output_dir,
  "mgcv_pca_all_subjects_performance_all_splits.tsv"
)
combined_prediction_file <- file.path(
  output_dir,
  "mgcv_pca_all_subjects_predictions_all_splits.tsv"
)
combined_coefficient_file <- file.path(
  output_dir,
  "mgcv_pca_all_subjects_coefficients_all_splits.tsv"
)
combined_feature_qc_file <- file.path(
  output_dir,
  "mgcv_pca_all_subjects_feature_qc_all_splits.tsv"
)
combined_split_qc_file <- file.path(
  output_dir,
  "mgcv_pca_all_subjects_split_qc_all_splits.tsv"
)
performance_summary_file <- file.path(
  output_dir,
  "mgcv_pca_all_subjects_performance_summary_across_splits.tsv"
)
paired_difference_file <- file.path(
  output_dir,
  "mgcv_pca_all_subjects_vs_covariate_paired_differences.tsv"
)
corrected_comparison_file <- file.path(
  output_dir,
  "mgcv_pca_all_subjects_vs_covariate_corrected_repeated_cv.tsv"
)
subject_average_file <- file.path(
  output_dir,
  "mgcv_pca_all_subjects_subject_average_test_predictions.tsv"
)
subject_average_performance_file <- file.path(
  output_dir,
  "mgcv_pca_all_subjects_subject_average_performance.tsv"
)
roc_coordinates_file <- file.path(
  output_dir,
  "mgcv_pca_all_subjects_subject_average_ROC_coordinates.tsv"
)
pr_coordinates_file <- file.path(
  output_dir,
  "mgcv_pca_all_subjects_subject_average_PR_coordinates.tsv"
)
summary_text_file <- file.path(
  output_dir,
  "mgcv_pca_all_subjects_summary.txt"
)
training_comparison_status_file <- file.path(
  qc_dir,
  "training_only_comparison_status.tsv"
)
training_paired_file <- file.path(
  output_dir,
  "mgcv_pca_all_subjects_vs_training_only_paired_differences.tsv"
)
training_corrected_file <- file.path(
  output_dir,
  "mgcv_pca_all_subjects_vs_training_only_corrected_repeated_cv.tsv"
)

all_outputs <- c(
  combined_performance_file,
  combined_prediction_file,
  combined_coefficient_file,
  combined_feature_qc_file,
  combined_split_qc_file,
  performance_summary_file,
  paired_difference_file,
  corrected_comparison_file,
  subject_average_file,
  subject_average_performance_file,
  roc_coordinates_file,
  pr_coordinates_file,
  summary_text_file,
  training_comparison_status_file
)

if (!overwrite && any(file.exists(all_outputs))) {
  stop(
    "Collected output already exists. Set OVERWRITE=1 to replace it. First: ",
    all_outputs[file.exists(all_outputs)][1L]
  )
}

if (overwrite) {
  unlink(
    c(all_outputs, training_paired_file, training_corrected_file),
    force = TRUE
  )
}

setorder(performance_all, splitID, model)
setorder(prediction_all, splitID, model, prediction_set, ID)
setorder(coefficient_all, splitID, model, feature)
if ("splitID" %in% names(feature_qc_all)) setorder(feature_qc_all, splitID)
setorder(split_qc_all, splitID)

fwrite(performance_all, combined_performance_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(prediction_all, combined_prediction_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(coefficient_all, combined_coefficient_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(feature_qc_all, combined_feature_qc_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(split_qc_all, combined_split_qc_file, sep = "\t", quote = FALSE, na = "NA")

# ------------------------- descriptive split summaries -----------------------

summary_metrics <- c(
  "AUROC_inner_oof",
  "PR_AUC_inner_oof",
  "Brier_inner_oof",
  "AUROC_test",
  "PR_AUC_test",
  "PR_AUC_lift_test",
  "Brier_test",
  "log_loss_test",
  "calibration_intercept_test",
  "calibration_slope_test",
  "sensitivity_test",
  "specificity_test",
  "balanced_accuracy_test",
  "n_DMR_selected_input",
  "n_DMR_usable",
  "n_methylation_features_input",
  "n_methylation_features_nonzero"
)
summary_metrics <- summary_metrics[summary_metrics %in% names(performance_all)]

summary_group_columns <- c("DMR_method", "pca_scope", "lambda_rule", "model")
summary_group_columns <- summary_group_columns[
  summary_group_columns %in% names(performance_all)
]

performance_summary <- performance_all[
  ,
  {
    result <- list(
      n_splits = .N,
      n_fit_fallback = if ("model_fallback" %in% names(.SD)) {
        sum(model_fallback, na.rm = TRUE)
      } else {
        NA_integer_
      }
    )

    for (metric in summary_metrics) {
      result <- c(result, summarize_numeric(get(metric), metric))
    }

    result
  },
  by = summary_group_columns
]

fwrite(
  performance_summary,
  performance_summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# -------------------------- paired baseline comparison -----------------------

paired <- merge(
  primary_rows,
  baseline_rows,
  by = "splitID",
  suffixes = c("_PCA_all_subjects", "_Covariate"),
  all = FALSE
)

paired[, test_train_ratio := n_test_PCA_all_subjects / n_train_PCA_all_subjects]

comparison_metrics <- c(
  "AUROC_test",
  "PR_AUC_test",
  "PR_AUC_lift_test",
  "Brier_test",
  "log_loss_test",
  "calibration_slope_test",
  "balanced_accuracy_test"
)
comparison_metrics <- comparison_metrics[
  paste0(comparison_metrics, "_PCA_all_subjects") %in% names(paired) &
    paste0(comparison_metrics, "_Covariate") %in% names(paired)
]

for (metric in comparison_metrics) {
  paired[
    ,
    (paste0("difference_", metric)) :=
      get(paste0(metric, "_PCA_all_subjects")) -
      get(paste0(metric, "_Covariate"))
  ]
}

paired_output_columns <- c(
  "splitID",
  "test_train_ratio",
  paste0(comparison_metrics, "_PCA_all_subjects"),
  paste0(comparison_metrics, "_Covariate"),
  paste0("difference_", comparison_metrics),
  "n_DMR_selected_input_PCA_all_subjects",
  "n_DMR_usable_PCA_all_subjects",
  "n_methylation_features_input_PCA_all_subjects",
  "n_methylation_features_nonzero_PCA_all_subjects",
  "model_fallback_PCA_all_subjects"
)
paired_output_columns <- paired_output_columns[
  paired_output_columns %in% names(paired)
]

fwrite(
  paired[, ..paired_output_columns],
  paired_difference_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

corrected_results <- make_corrected_table(
  paired = paired,
  metrics = comparison_metrics,
  difference_prefix = "difference_",
  comparison_label = paste0(primary_model_name, "_PCA_all_subjects_minus_Covariate_only")
)

fwrite(
  corrected_results,
  corrected_comparison_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ---------------------- subject-averaged test predictions --------------------

prediction_test <- prediction_all[prediction_set == "test"]

outcome_check <- prediction_test[
  ,
  .(n_outcomes = uniqueN(AA_only)),
  by = .(model, ID, FID)
]

if (any(outcome_check$n_outcomes != 1L)) {
  stop("At least one subject has inconsistent outcomes across test appearances.")
}

subject_average <- prediction_test[
  ,
  .(
    AA_only = unique(AA_only)[1L],
    mean_pred_prob = mean(pred_prob, na.rm = TRUE),
    median_pred_prob = median(pred_prob, na.rm = TRUE),
    sd_pred_prob = stats::sd(pred_prob, na.rm = TRUE),
    n_test_appearances = .N,
    mean_split_threshold = mean(threshold, na.rm = TRUE)
  ),
  by = .(DMR_method, model, ID, FID)
]

setorder(subject_average, model, ID)

fwrite(
  subject_average,
  subject_average_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

subject_average_performance <- subject_average[
  ,
  aggregate_probability_metrics(
    y = AA_only,
    p = mean_pred_prob
  ),
  by = .(DMR_method, model)
]

subject_average_performance[, `:=`(
  pca_scope = expected_pca_scope,
  lambda_rule = fifelse(model == primary_model_name, analysis_lambda_rule, NA_character_)
)]

setcolorder(
  subject_average_performance,
  c(
    "DMR_method", "pca_scope", "lambda_rule", "model",
    setdiff(
      names(subject_average_performance),
      c("DMR_method", "pca_scope", "lambda_rule", "model")
    )
  )
)

fwrite(
  subject_average_performance,
  subject_average_performance_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------ curve outputs --------------------------------

roc_list <- list()
pr_list <- list()

for (model_name in unique(subject_average$model)) {
  dd <- subject_average[model == model_name]

  if (length(unique(dd$AA_only)) < 2L) next

  roc_obj <- pROC::roc(
    dd$AA_only,
    dd$mean_pred_prob,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  roc_list[[model_name]] <- data.table(
    model = model_name,
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    AUROC = as.numeric(pROC::auc(roc_obj))
  )

  pr_obj <- PRROC::pr.curve(
    scores.class0 = dd$mean_pred_prob[dd$AA_only == 1L],
    scores.class1 = dd$mean_pred_prob[dd$AA_only == 0L],
    curve = TRUE
  )

  pr_list[[model_name]] <- data.table(
    model = model_name,
    Recall = pr_obj$curve[, 1L],
    Precision = pr_obj$curve[, 2L],
    PR_AUC = as.numeric(pr_obj$auc.integral)
  )
}

roc_coordinates <- rbindlist(roc_list, use.names = TRUE, fill = TRUE)
pr_coordinates <- rbindlist(pr_list, use.names = TRUE, fill = TRUE)

fwrite(roc_coordinates, roc_coordinates_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(pr_coordinates, pr_coordinates_file, sep = "\t", quote = FALSE, na = "NA")

roc_coordinates[, label := paste0(model, " AUC=", sprintf("%.3f", AUROC))]
pr_coordinates[, label := paste0(model, " PR-AUC=", sprintf("%.3f", PR_AUC))]

p_roc <- ggplot(
  roc_coordinates,
  aes(x = FPR, y = TPR, color = label)
) +
  geom_line(linewidth = 1.1) +
  geom_abline(linetype = "dashed") +
  theme_bw(base_size = 13) +
  labs(
    title = "Intentional all-subject PCA: subject-averaged ROC",
    subtitle = "Sensitivity analysis with outer-test information in PCA",
    x = "False positive rate",
    y = "True positive rate",
    color = "Model"
  )

p_pr <- ggplot(
  pr_coordinates,
  aes(x = Recall, y = Precision, color = label)
) +
  geom_line(linewidth = 1.1) +
  theme_bw(base_size = 13) +
  labs(
    title = "Intentional all-subject PCA: subject-averaged precision-recall",
    subtitle = "Sensitivity analysis with outer-test information in PCA",
    x = "Recall",
    y = "Precision",
    color = "Model"
  )

ggsave(
  file.path(figure_dir, "mgcv_pca_all_subjects_subject_average_ROC.pdf"),
  p_roc,
  width = 8,
  height = 6
)

ggsave(
  file.path(figure_dir, "mgcv_pca_all_subjects_subject_average_ROC.png"),
  p_roc,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(figure_dir, "mgcv_pca_all_subjects_subject_average_PR.pdf"),
  p_pr,
  width = 8,
  height = 6
)

ggsave(
  file.path(figure_dir, "mgcv_pca_all_subjects_subject_average_PR.png"),
  p_pr,
  width = 8,
  height = 6,
  dpi = 300
)

# Split-level performance distribution figure.
plot_metrics <- c("AUROC_test", "PR_AUC_test", "Brier_test")
plot_metrics <- plot_metrics[plot_metrics %in% names(performance_all)]

plot_dt <- melt(
  performance_all[model %in% c(primary_model_name, "Covariate_only")],
  id.vars = c("splitID", "model"),
  measure.vars = plot_metrics,
  variable.name = "metric",
  value.name = "value"
)

p_split <- ggplot(
  plot_dt,
  aes(x = model, y = value, fill = model)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.25, size = 0.8) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Split-level performance: intentional all-subject PCA",
    subtitle = "Outer-test subjects contributed to PCA estimation",
    x = NULL,
    y = NULL
  )

ggsave(
  file.path(figure_dir, "mgcv_pca_all_subjects_split_performance.pdf"),
  p_split,
  width = 10,
  height = 5.5
)

ggsave(
  file.path(figure_dir, "mgcv_pca_all_subjects_split_performance.png"),
  p_split,
  width = 10,
  height = 5.5,
  dpi = 300
)

# ------------- optional comparison with proper training-only PCA -------------

training_comparison_status <- data.table(
  attempted = compare_with_training_only,
  reference_file = training_only_performance_file,
  reference_exists = file.exists(training_only_performance_file),
  current_lambda_rule = analysis_lambda_rule,
  reference_lambda_rule = NA_character_,
  lambda_rules_match = NA,
  comparison_completed = FALSE,
  status = "NOT_ATTEMPTED"
)

training_corrected_results <- NULL

if (compare_with_training_only) {
  if (!file.exists(training_only_performance_file)) {
    training_comparison_status[, status := "REFERENCE_FILE_NOT_FOUND"]
  } else {
    reference_perf <- read_one_table(
      training_only_performance_file,
      "Training-only reference performance"
    )

    reference_primary <- reference_perf[model == primary_model_name]

    if (nrow(reference_primary) == 0L) {
      training_comparison_status[, status := "REFERENCE_PRIMARY_MODEL_NOT_FOUND"]
    } else if (anyDuplicated(reference_primary$splitID)) {
      training_comparison_status[, status := "REFERENCE_HAS_DUPLICATED_SPLIT_ROWS"]
    } else {
      reference_lambda_rules <- if ("lambda_rule" %in% names(reference_primary)) {
        unique(na.omit(as.character(reference_primary$lambda_rule)))
      } else {
        character()
      }

      reference_lambda_rule <- if (length(reference_lambda_rules) == 1L) {
        reference_lambda_rules
      } else {
        NA_character_
      }

      lambda_match <- identical(
        as.character(analysis_lambda_rule),
        as.character(reference_lambda_rule)
      )

      training_comparison_status[, `:=`(
        reference_lambda_rule = reference_lambda_rule,
        lambda_rules_match = lambda_match
      )]

      if (!lambda_match && !allow_lambda_mismatch) {
        training_comparison_status[, status := "SKIPPED_LAMBDA_RULE_MISMATCH"]
      } else {
        reference_required <- c(
          "splitID", "n_train", "n_test", comparison_metrics
        )
        missing_reference <- setdiff(reference_required, names(reference_primary))

        if (length(missing_reference) > 0L) {
          training_comparison_status[, status := paste0(
            "REFERENCE_MISSING_COLUMNS:",
            paste(missing_reference, collapse = ",")
          )]
        } else {
          current_for_compare <- copy(primary_rows)
          reference_for_compare <- copy(reference_primary)

          pca_compare <- merge(
            current_for_compare,
            reference_for_compare,
            by = "splitID",
            suffixes = c("_PCA_all_subjects", "_PCA_training_only"),
            all = FALSE
          )

          pca_compare[, test_train_ratio :=
            n_test_PCA_all_subjects / n_train_PCA_all_subjects]

          for (metric in comparison_metrics) {
            pca_compare[
              ,
              (paste0("difference_", metric)) :=
                get(paste0(metric, "_PCA_all_subjects")) -
                get(paste0(metric, "_PCA_training_only"))
            ]
          }

          keep_columns <- c(
            "splitID",
            "test_train_ratio",
            paste0(comparison_metrics, "_PCA_all_subjects"),
            paste0(comparison_metrics, "_PCA_training_only"),
            paste0("difference_", comparison_metrics),
            "lambda_rule_PCA_all_subjects",
            "lambda_rule_PCA_training_only",
            "n_methylation_features_nonzero_PCA_all_subjects",
            "n_methylation_features_nonzero_PCA_training_only"
          )
          keep_columns <- keep_columns[keep_columns %in% names(pca_compare)]

          fwrite(
            pca_compare[, ..keep_columns],
            training_paired_file,
            sep = "\t",
            quote = FALSE,
            na = "NA"
          )

          training_corrected_results <- make_corrected_table(
            paired = pca_compare,
            metrics = comparison_metrics,
            difference_prefix = "difference_",
            comparison_label = paste0(
              primary_model_name,
              "_PCA_all_subjects_minus_PCA_training_only"
            )
          )

          fwrite(
            training_corrected_results,
            training_corrected_file,
            sep = "\t",
            quote = FALSE,
            na = "NA"
          )

          training_comparison_status[, `:=`(
            comparison_completed = TRUE,
            status = if (lambda_match) {
              "COMPLETED_SAME_LAMBDA_RULE"
            } else {
              "COMPLETED_WITH_LAMBDA_MISMATCH_ALLOWED"
            }
          )]

          # Paired delta figure. Positive favors all-subject PCA for AUROC/PR-AUC;
          # negative favors all-subject PCA for Brier.
          delta_metrics <- c("AUROC_test", "PR_AUC_test", "Brier_test")
          delta_metrics <- delta_metrics[
            paste0("difference_", delta_metrics) %in% names(pca_compare)
          ]

          delta_dt <- melt(
            pca_compare,
            id.vars = "splitID",
            measure.vars = paste0("difference_", delta_metrics),
            variable.name = "metric",
            value.name = "difference"
          )
          delta_dt[, metric := sub("^difference_", "", metric)]

          p_delta <- ggplot(
            delta_dt,
            aes(x = metric, y = difference)
          ) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.12, alpha = 0.3, size = 0.8) +
            theme_bw(base_size = 12) +
            labs(
              title = "All-subject PCA minus training-only PCA",
              subtitle = paste0(
                "Intentional leakage sensitivity; lambda rule = ",
                analysis_lambda_rule
              ),
              x = NULL,
              y = "Paired split-level difference"
            )

          ggsave(
            file.path(
              figure_dir,
              "mgcv_pca_all_subjects_vs_training_only_differences.pdf"
            ),
            p_delta,
            width = 8,
            height = 5.5
          )

          ggsave(
            file.path(
              figure_dir,
              "mgcv_pca_all_subjects_vs_training_only_differences.png"
            ),
            p_delta,
            width = 8,
            height = 5.5,
            dpi = 300
          )
        }
      }
    }
  }
}

fwrite(
  training_comparison_status,
  training_comparison_status_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ----------------------------- text summary ----------------------------------

primary_summary <- performance_all[
  model == primary_model_name,
  .(
    n_splits = .N,
    mean_AUROC_test = mean(AUROC_test, na.rm = TRUE),
    sd_AUROC_test = stats::sd(AUROC_test, na.rm = TRUE),
    median_AUROC_test = median(AUROC_test, na.rm = TRUE),
    mean_PR_AUC_test = mean(PR_AUC_test, na.rm = TRUE),
    sd_PR_AUC_test = stats::sd(PR_AUC_test, na.rm = TRUE),
    median_PR_AUC_test = median(PR_AUC_test, na.rm = TRUE),
    mean_Brier_test = mean(Brier_test, na.rm = TRUE),
    median_Brier_test = median(Brier_test, na.rm = TRUE),
    mean_nonzero_methylation_features = if (
      "n_methylation_features_nonzero" %in% names(performance_all)
    ) {
      mean(n_methylation_features_nonzero, na.rm = TRUE)
    } else {
      NA_real_
    },
    zero_methylation_feature_splits = if (
      "n_methylation_features_nonzero" %in% names(performance_all)
    ) {
      sum(n_methylation_features_nonzero == 0L, na.rm = TRUE)
    } else {
      NA_integer_
    }
  )
]

baseline_summary <- performance_all[
  model == "Covariate_only",
  .(
    mean_AUROC_test = mean(AUROC_test, na.rm = TRUE),
    mean_PR_AUC_test = mean(PR_AUC_test, na.rm = TRUE),
    mean_Brier_test = mean(Brier_test, na.rm = TRUE)
  )
]

summary_lines <- c(
  "MGCV prediction sensitivity analysis: PCA fitted on all outer-split subjects",
  "==========================================================================",
  paste0("Input root: ", input_root),
  paste0("PCA scope: ", expected_pca_scope),
  paste0("Lambda rule: ", analysis_lambda_rule),
  paste0("Splits collected: ", primary_summary$n_splits),
  "",
  "IMPORTANT INTERPRETATION:",
  paste(
    "The outer-test subjects contributed to PCA estimation. These results",
    "contain intentional information leakage and are for sensitivity/diagnostic",
    "purposes only. They are not unbiased estimates of predictive performance."
  ),
  "",
  paste0(
    primary_model_name,
    " mean test AUROC = ", sprintf("%.4f", primary_summary$mean_AUROC_test),
    " (SD ", sprintf("%.4f", primary_summary$sd_AUROC_test), ")"
  ),
  paste0(
    primary_model_name,
    " median test AUROC = ", sprintf("%.4f", primary_summary$median_AUROC_test)
  ),
  paste0(
    primary_model_name,
    " mean test PR-AUC = ", sprintf("%.4f", primary_summary$mean_PR_AUC_test),
    " (SD ", sprintf("%.4f", primary_summary$sd_PR_AUC_test), ")"
  ),
  paste0(
    primary_model_name,
    " mean test Brier score = ", sprintf("%.4f", primary_summary$mean_Brier_test)
  ),
  paste0(
    "Mean number of nonzero methylation PCs = ",
    sprintf("%.2f", primary_summary$mean_nonzero_methylation_features)
  ),
  paste0(
    "Splits with zero nonzero methylation PCs = ",
    primary_summary$zero_methylation_feature_splits
  ),
  "",
  paste0(
    "Covariate-only mean test AUROC = ",
    sprintf("%.4f", baseline_summary$mean_AUROC_test)
  ),
  paste0(
    "Covariate-only mean test PR-AUC = ",
    sprintf("%.4f", baseline_summary$mean_PR_AUC_test)
  ),
  paste0(
    "Covariate-only mean test Brier score = ",
    sprintf("%.4f", baseline_summary$mean_Brier_test)
  ),
  "",
  paste0(
    "Training-only PCA comparison status: ",
    training_comparison_status$status
  )
)

if (!is.null(training_corrected_results)) {
  auroc_delta <- training_corrected_results[metric == "AUROC_test"]
  pr_delta <- training_corrected_results[metric == "PR_AUC_test"]
  brier_delta <- training_corrected_results[metric == "Brier_test"]

  if (nrow(auroc_delta) == 1L) {
    summary_lines <- c(
      summary_lines,
      "",
      "All-subject PCA minus training-only PCA, corrected repeated-CV:",
      paste0(
        "AUROC difference = ", sprintf("%.4f", auroc_delta$mean_difference),
        "; 95% CI [", sprintf("%.4f", auroc_delta$corrected_CI_low),
        ", ", sprintf("%.4f", auroc_delta$corrected_CI_high), "]"
      )
    )
  }

  if (nrow(pr_delta) == 1L) {
    summary_lines <- c(
      summary_lines,
      paste0(
        "PR-AUC difference = ", sprintf("%.4f", pr_delta$mean_difference),
        "; 95% CI [", sprintf("%.4f", pr_delta$corrected_CI_low),
        ", ", sprintf("%.4f", pr_delta$corrected_CI_high), "]"
      )
    )
  }

  if (nrow(brier_delta) == 1L) {
    summary_lines <- c(
      summary_lines,
      paste0(
        "Brier difference = ", sprintf("%.4f", brier_delta$mean_difference),
        "; 95% CI [", sprintf("%.4f", brier_delta$corrected_CI_low),
        ", ", sprintf("%.4f", brier_delta$corrected_CI_high), "]",
        "; negative favors all-subject PCA"
      )
    )
  }
}

writeLines(summary_lines, summary_text_file)

# -------------------------------- final report -------------------------------

cat("\nMGCV all-subject PCA prediction collection finished.\n")
cat("Combined performance:", combined_performance_file, "\n")
cat("Combined predictions:", combined_prediction_file, "\n")
cat("Performance summary:", performance_summary_file, "\n")
cat("Paired differences versus covariates:", paired_difference_file, "\n")
cat("Corrected comparison versus covariates:", corrected_comparison_file, "\n")
cat("Subject-averaged performance:", subject_average_performance_file, "\n")
cat("Summary text:", summary_text_file, "\n")
cat("Training-only comparison status:", training_comparison_status$status, "\n")
cat("Figures:", figure_dir, "\n\n")

print(
  performance_all[
    ,
    .(
      n_splits = .N,
      mean_AUROC_test = mean(AUROC_test, na.rm = TRUE),
      median_AUROC_test = median(AUROC_test, na.rm = TRUE),
      mean_PR_AUC_test = mean(PR_AUC_test, na.rm = TRUE),
      median_PR_AUC_test = median(PR_AUC_test, na.rm = TRUE),
      mean_Brier_test = mean(Brier_test, na.rm = TRUE),
      median_Brier_test = median(Brier_test, na.rm = TRUE),
      mean_DMRs_usable = if ("n_DMR_usable" %in% names(performance_all)) {
        mean(n_DMR_usable, na.rm = TRUE)
      } else {
        NA_real_
      },
      mean_nonzero_methylation_features = if (
        "n_methylation_features_nonzero" %in% names(performance_all)
      ) {
        mean(n_methylation_features_nonzero, na.rm = TRUE)
      } else {
        NA_real_
      },
      fallback_splits = if ("model_fallback" %in% names(performance_all)) {
        sum(model_fallback, na.rm = TRUE)
      } else {
        NA_integer_
      }
    ),
    by = .(model, pca_scope, lambda_rule)
  ]
)

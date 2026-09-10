#!/usr/bin/env Rscript

# Collect the 100 split-specific outputs from
# 5_run_dmrcate_primary_prediction_cv.R.
#
# This script creates:
#   - combined performance, prediction, coefficient, feature-QC, and split-QC
#     files;
#   - descriptive performance summaries across splits;
#   - paired DMRcate-versus-covariate-only differences;
#   - corrected repeated-CV confidence intervals for paired differences;
#   - subject-level averages of repeated out-of-fold test predictions;
#   - summary ROC/PR curve coordinates and figures.

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(PRROC)
  library(ggplot2)
})

data.table::setDTthreads(1L)

# ------------------------------- configuration -------------------------------

PATH_wk <- path.expand(
  Sys.getenv("PATH_WK", "~/scratch/UQAC/meth/")
)

input_root <- Sys.getenv(
  "DMRCATE_PREDICTION_OUTPUT",
  file.path(
    PATH_wk,
    "results/15_revision/5_cv/5_prediction/3_dmrcate_primary"
  )
)

performance_dir <- file.path(input_root, "performance")
prediction_dir <- file.path(input_root, "predictions")
coefficient_dir <- file.path(input_root, "coefficients")
feature_qc_dir <- file.path(input_root, "feature_qc")
split_qc_dir <- file.path(input_root, "split_qc")

output_dir <- file.path(input_root, "collected")
figure_dir <- file.path(output_dir, "figures")

for (d in c(output_dir, figure_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

expected_splits <- suppressWarnings(
  as.integer(Sys.getenv("EXPECTED_SPLITS", "100"))
)

overwrite <- Sys.getenv("OVERWRITE", "0") %in%
  c("1", "TRUE", "true", "T", "yes", "YES")

if (is.na(expected_splits) || expected_splits < 1L) {
  stop("EXPECTED_SPLITS must be a positive integer.")
}

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
      y,
      p,
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

  data.table(
    n_subjects = length(y),
    n_cases = sum(y == 1L),
    n_controls = sum(y == 0L),
    prevalence = mean(y),
    AUROC = safe_roc_auc(y, p),
    PR_AUC = safe_pr_auc(y, p),
    PR_baseline = mean(y),
    PR_AUC_lift = safe_pr_auc(y, p) - mean(y),
    Brier = mean((p - y)^2),
    log_loss = -mean(y * log(p_clip) + (1 - y) * log(1 - p_clip)),
    calibration_intercept = calibration_intercept,
    calibration_slope = calibration_slope
  )
}

summarize_numeric <- function(x, prefix) {
  x <- x[is.finite(x)]

  statistic_names <- c(
    "mean", "sd", "median", "q025", "q25", "q75", "q975", "min"
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
      min = min(x)
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
        win_rate = if (n > 0L) mean(d > 0) else NA_real_
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
    win_rate = mean(d > 0)
  )
}

repair_duplicate_column_names <- function(dt, file_path, label) {
  current_names <- names(dt)

  if (!anyDuplicated(current_names)) {
    return(dt)
  }

  # Older prediction files contain two columns called n_test:
  #   1. total number of subjects in the outer test set;
  #   2. number of evaluable test predictions returned by calculate_metrics().
  # Preserve the first as n_test and rename the metric count.
  n_test_index <- which(current_names == "n_test")

  if (length(n_test_index) > 1L) {
    replacement_names <- c(
      "n_evaluable_test",
      if (length(n_test_index) > 2L) {
        paste0("n_evaluable_test_", seq_len(length(n_test_index) - 2L) + 1L)
      } else {
        character()
      }
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
      " file still has duplicated column names after repair: ",
      basename(file_path),
      " [",
      paste(duplicated_names, collapse = ", "),
      "]"
    )
  }

  setnames(dt, current_names)
  dt
}

read_and_bind <- function(files, label) {
  if (length(files) == 0L) {
    stop("No ", label, " files were found.")
  }

  table_list <- lapply(
    files,
    function(file_path) {
      dt <- fread(file_path)
      repair_duplicate_column_names(
        dt = dt,
        file_path = file_path,
        label = label
      )
    }
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

# -------------------------------- input files --------------------------------

performance_files <- list.files(
  performance_dir,
  pattern = "^dmrcate_primary_performance_split_[0-9]{3}\\.tsv$",
  full.names = TRUE
)

prediction_files <- list.files(
  prediction_dir,
  pattern = "^dmrcate_primary_predictions_split_[0-9]{3}\\.tsv$",
  full.names = TRUE
)

coefficient_files <- list.files(
  coefficient_dir,
  pattern = "^dmrcate_primary_coefficients_split_[0-9]{3}\\.tsv$",
  full.names = TRUE
)

feature_qc_files <- list.files(
  feature_qc_dir,
  pattern = "^dmrcate_primary_feature_qc_split_[0-9]{3}\\.tsv$",
  full.names = TRUE
)

split_qc_files <- list.files(
  split_qc_dir,
  pattern = "^dmrcate_primary_split_qc_[0-9]{3}\\.tsv$",
  full.names = TRUE
)

cat(
  "Found performance/prediction/coefficient/feature-QC/split-QC files:",
  length(performance_files),
  length(prediction_files),
  length(coefficient_files),
  length(feature_qc_files),
  length(split_qc_files),
  "\n"
)

if (length(performance_files) != expected_splits ||
    length(prediction_files) != expected_splits ||
    length(split_qc_files) != expected_splits) {
  stop(
    "Expected ", expected_splits,
    " complete performance, prediction, and split-QC files before collection."
  )
}

performance_all <- read_and_bind(performance_files, "Performance")
prediction_all <- read_and_bind(prediction_files, "Prediction")
coefficient_all <- read_and_bind(coefficient_files, "Coefficient")
feature_qc_all <- read_and_bind(feature_qc_files, "Feature QC")
split_qc_all <- read_and_bind(split_qc_files, "Split QC")

if (uniqueN(performance_all$splitID) != expected_splits ||
    uniqueN(prediction_all$splitID) != expected_splits) {
  stop("Not all expected splitIDs are represented in the collected results.")
}

# -------------------------------- output files -------------------------------

combined_performance_file <- file.path(
  output_dir,
  "dmrcate_primary_performance_all_splits.tsv"
)

combined_prediction_file <- file.path(
  output_dir,
  "dmrcate_primary_predictions_all_splits.tsv"
)

combined_coefficient_file <- file.path(
  output_dir,
  "dmrcate_primary_coefficients_all_splits.tsv"
)

combined_feature_qc_file <- file.path(
  output_dir,
  "dmrcate_primary_feature_qc_all_splits.tsv"
)

combined_split_qc_file <- file.path(
  output_dir,
  "dmrcate_primary_split_qc_all_splits.tsv"
)

performance_summary_file <- file.path(
  output_dir,
  "dmrcate_primary_performance_summary_across_splits.tsv"
)

paired_difference_file <- file.path(
  output_dir,
  "dmrcate_vs_covariate_paired_differences.tsv"
)

corrected_comparison_file <- file.path(
  output_dir,
  "dmrcate_vs_covariate_corrected_repeated_cv.tsv"
)

subject_average_file <- file.path(
  output_dir,
  "dmrcate_primary_subject_average_test_predictions.tsv"
)

subject_average_performance_file <- file.path(
  output_dir,
  "dmrcate_primary_subject_average_performance.tsv"
)

roc_coordinates_file <- file.path(
  output_dir,
  "dmrcate_primary_subject_average_ROC_coordinates.tsv"
)

pr_coordinates_file <- file.path(
  output_dir,
  "dmrcate_primary_subject_average_PR_coordinates.tsv"
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
  pr_coordinates_file
)

if (!overwrite && any(file.exists(all_outputs))) {
  stop(
    "Collected output already exists. Set OVERWRITE=1 to replace it. First: ",
    all_outputs[file.exists(all_outputs)][1L]
  )
}

if (overwrite) unlink(all_outputs, force = TRUE)

setorder(performance_all, splitID, model)
setorder(prediction_all, splitID, model, prediction_set, ID)
setorder(coefficient_all, splitID, model, feature)
setorder(split_qc_all, splitID)

fwrite(performance_all, combined_performance_file, sep = "\t", na = "NA")
fwrite(prediction_all, combined_prediction_file, sep = "\t", na = "NA")
fwrite(coefficient_all, combined_coefficient_file, sep = "\t", na = "NA")
fwrite(feature_qc_all, combined_feature_qc_file, sep = "\t", na = "NA")
fwrite(split_qc_all, combined_split_qc_file, sep = "\t", na = "NA")

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

summary_list <- performance_all[
  ,
  {
    result <- list(
      n_splits = .N,
      n_fit_fallback = sum(model_fallback, na.rm = TRUE)
    )

    for (metric in summary_metrics) {
      result <- c(
        result,
        summarize_numeric(get(metric), metric)
      )
    }

    result
  },
  by = .(
    DMR_method,
    model
  )
]

fwrite(
  summary_list,
  performance_summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# -------------------------- paired baseline comparison -----------------------

dmrcate_perf <- performance_all[model == "DMRcate_PC12_LASSO"]
baseline_perf <- performance_all[model == "Covariate_only"]

paired <- merge(
  dmrcate_perf,
  baseline_perf,
  by = "splitID",
  suffixes = c("_DMRcate", "_Covariate"),
  all = FALSE
)

paired[, test_train_ratio := n_test_DMRcate / n_train_DMRcate]

comparison_metrics <- c(
  "AUROC_test",
  "PR_AUC_test",
  "PR_AUC_lift_test",
  "Brier_test",
  "log_loss_test",
  "calibration_slope_test",
  "balanced_accuracy_test"
)

for (metric in comparison_metrics) {
  paired[
    ,
    (paste0("difference_", metric)) :=
      get(paste0(metric, "_DMRcate")) -
      get(paste0(metric, "_Covariate"))
  ]
}

paired_output_columns <- c(
  "splitID",
  "test_train_ratio",
  paste0(comparison_metrics, "_DMRcate"),
  paste0(comparison_metrics, "_Covariate"),
  paste0("difference_", comparison_metrics),
  "n_DMR_selected_input_DMRcate",
  "n_DMR_usable_DMRcate",
  "n_methylation_features_input_DMRcate",
  "n_methylation_features_nonzero_DMRcate",
  "model_fallback_DMRcate"
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

corrected_results <- rbindlist(
  lapply(
    comparison_metrics,
    function(metric) {
      result <- corrected_repeated_cv(
        difference = paired[[paste0("difference_", metric)]],
        test_train_ratio = paired$test_train_ratio
      )

      lower_is_better <- metric %in% c("Brier_test", "log_loss_test")
      d_values <- paired[[paste0("difference_", metric)]]

      result[, `:=`(
        comparison = "DMRcate_PC12_LASSO_minus_Covariate_only",
        metric = metric,
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
          "favorable_direction",
          setdiff(names(result), c("comparison", "metric", "favorable_direction"))
        )
      )

      result
    }
  ),
  use.names = TRUE,
  fill = TRUE
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
  .(
    n_outcomes = uniqueN(AA_only)
  ),
  by = .(
    model,
    ID,
    FID
  )
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
  by = .(
    DMR_method,
    model,
    ID,
    FID
  )
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
  by = .(
    DMR_method,
    model
  )
]

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

roc_coordinates <- rbindlist(roc_list)
pr_coordinates <- rbindlist(pr_list)

fwrite(roc_coordinates, roc_coordinates_file, sep = "\t", na = "NA")
fwrite(pr_coordinates, pr_coordinates_file, sep = "\t", na = "NA")

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
    title = "Subject-averaged repeated test predictions: ROC",
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
    title = "Subject-averaged repeated test predictions: precision-recall",
    x = "Recall",
    y = "Precision",
    color = "Model"
  )

ggsave(
  file.path(figure_dir, "dmrcate_subject_average_ROC.pdf"),
  p_roc,
  width = 8,
  height = 6
)

ggsave(
  file.path(figure_dir, "dmrcate_subject_average_ROC.png"),
  p_roc,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(figure_dir, "dmrcate_subject_average_PR.pdf"),
  p_pr,
  width = 8,
  height = 6
)

ggsave(
  file.path(figure_dir, "dmrcate_subject_average_PR.png"),
  p_pr,
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------------- final report -------------------------------

cat("\nDMRcate primary prediction collection finished.\n")
cat("Combined performance:", combined_performance_file, "\n")
cat("Combined predictions:", combined_prediction_file, "\n")
cat("Performance summary:", performance_summary_file, "\n")
cat("Paired differences:", paired_difference_file, "\n")
cat("Corrected repeated-CV comparison:", corrected_comparison_file, "\n")
cat("Subject-averaged performance:", subject_average_performance_file, "\n")
cat("Figures:", figure_dir, "\n")

print(
  performance_all[
    ,
    .(
      n_splits = .N,
      median_AUROC_test = median(AUROC_test, na.rm = TRUE),
      median_PR_AUC_test = median(PR_AUC_test, na.rm = TRUE),
      median_Brier_test = median(Brier_test, na.rm = TRUE),
      mean_DMRs_usable = mean(n_DMR_usable, na.rm = TRUE),
      fallback_splits = sum(model_fallback, na.rm = TRUE)
    ),
    by = model
  ]
)

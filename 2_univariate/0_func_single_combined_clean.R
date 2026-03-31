suppressPackageStartupMessages({
  library(lmerTest)
  library(glmmLasso)
  library(lme4)
})

lambda_grid <- c(5, 10, 20, 40, 50, 60, 80, 100)

covariates <- c(
  "AgeCalc", "Sex", "Non.smoker",
  "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc",
  "sv1", "sv2", "sv3", "sv4", "sv5", "BMI")

features_model1 <- c("meth_response", covariates)
features_model2 <- c("AA_only", covariates)
features_model3 <- c("AA_only", covariates)

fit_glmmlasso_with_bic <- function(formula, data, family, lambda_grid) {
  bic_values <- rep(Inf, length(lambda_grid))

  for (idx in seq_along(lambda_grid)) {
    fit <- tryCatch(
      glmmLasso(
        fix = formula,
        rnd = list(FID = ~ 1),
        data = data,
        family = family,
        lambda = lambda_grid[idx]
      ),
      error = function(e) NULL
    )

    if (!is.null(fit) && !is.null(fit$bic) && is.finite(fit$bic)) {
      bic_values[idx] <- fit$bic
    }
  }

  if (all(!is.finite(bic_values))) {
    stop("All glmmLasso fits failed during lambda tuning.")
  }

  best_idx <- which.min(bic_values)
  best_lambda <- lambda_grid[best_idx]

  final_fit <- glmmLasso(
    fix = formula,
    rnd = list(FID = ~ 1),
    data = data,
    family = family,
    lambda = best_lambda,
    final.re = TRUE
  )

  list(
    fit = final_fit,
    optimal_lambda = best_lambda
  )
}

extract_glmmlasso_results <- function(fit_obj, features, coef_names) {
  coefficients <- coef(fit_obj)
  coefficient_template <- stats::setNames(rep(0, length(coef_names)), coef_names)
  coefficient_template[names(coefficients)] <- coefficients

  fixed_effects_table <- coef(summary(fit_obj))
  coef_pvals <- stats::setNames(rep(1, length(features)), features)

  matched_rows <- match(features, rownames(fixed_effects_table))
  valid_rows <- !is.na(matched_rows)
  coef_pvals[valid_rows] <- fixed_effects_table[matched_rows[valid_rows], "p.value"]
  coef_pvals[is.na(coef_pvals)] <- 1

  list(
    coefficients = coefficient_template,
    FID_sdv = as.numeric(fit_obj$StdDev),
    coef_pvals = coef_pvals
  )
}

default_model_result <- function(coef_names, pval_names) {
  list(
    coefficients = stats::setNames(rep(0, length(coef_names)), coef_names),
    FID_sdv = 0,
    coef_pvals = stats::setNames(rep(1, length(pval_names)), pval_names)
  )
}

coef_names_model1 <- c("(Intercept)", features_model1)
coef_names_model2 <- c("(Intercept)", features_model2)
coef_names_model3 <- c("(Intercept)", features_model3)

model1_result <- default_model_result(coef_names_model1, features_model1)
model2_result <- default_model_result(coef_names_model2, features_model2)
model3_result <- default_model_result(coef_names_model3, features_model3)

optimal_lambda_model1 <- 0
optimal_lambda_model3 <- 0

model1_formula <- as.formula(
  paste("AA_only ~ meth_response +", paste(covariates, collapse = " + "))
)

model3_formula <- as.formula(
  paste("meth_response ~ AA_only +", paste(covariates, collapse = " + "))
)

# Model 1: penalized logistic mixed model + stage-2 GLMM refit
model1_result <- tryCatch({
  model1_fit <- fit_glmmlasso_with_bic(
    formula = model1_formula,
    data = pheno_file_feed,
    family = binomial(link = "logit"),
    lambda_grid = lambda_grid
  )
  optimal_lambda_model1 <- model1_fit$optimal_lambda

  res <- extract_glmmlasso_results(
    fit_obj = model1_fit$fit,
    features = features_model1,
    coef_names = coef_names_model1
  )

  selected_features <- names(res$coefficients)[res$coefficients != 0]
  selected_features <- setdiff(selected_features, "(Intercept)")

  if (length(selected_features) > 0) {
    stage2_formula <- as.formula(
      paste0("AA_only ~ ", paste(selected_features, collapse = " + "), " + (1|FID)")
    )

    stage2_res <- tryCatch({
      glmod <- glmer(
        formula = stage2_formula,
        data = pheno_file_feed,
        family = binomial(link = "logit")
      )

      smry <- summary(glmod)$coefficients
      pvals <- smry[, "Pr(>|z|)"]
      ests <- smry[, "Estimate"]

      selected_in_summary <- intersect(selected_features, rownames(smry))
      res$coef_pvals[selected_in_summary] <- pvals[selected_in_summary]
      res$coefficients[selected_in_summary] <- ests[selected_in_summary]
      res$FID_sdv <- as.numeric(attr(VarCorr(glmod)$FID, "stddev"))

      res
    }, error = function(e) res)

    stage2_res
  } else {
    res
  }
}, error = function(e) {
  warning("Model 1 failed for CpG ", i_cpg, ": ", e$message)
  default_model_result(coef_names_model1, features_model1)
})

# Model 2: linear mixed model
model2_result <- tryCatch({
  lmer_fit <- lmer(
    meth_response ~ AA_only + AgeCalc + Sex + Non.smoker +
      EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc +
      sv1 + sv2 + sv3 + sv4 + sv5 + (1 | FID),
    data = pheno_file_feed
  )

  smry <- summary(lmer_fit)$coefficients
  ests <- smry[, "Estimate"]
  pvals <- smry[, "Pr(>|t|)"]
  pvals[is.na(pvals)] <- 1

  res <- default_model_result(coef_names_model2, features_model2)
  matched_coef <- intersect(names(res$coefficients), rownames(smry))
  matched_pval <- intersect(names(res$coef_pvals), rownames(smry))

  res$coefficients[matched_coef] <- ests[matched_coef]
  res$coef_pvals[matched_pval] <- pvals[matched_pval]
  res$FID_sdv <- as.numeric(attr(VarCorr(lmer_fit)$FID, "stddev"))
  res
}, error = function(e) {
  warning("Model 2 failed for CpG ", i_cpg, ": ", e$message)
  default_model_result(coef_names_model2, features_model2)
})

# Model 3: penalized Gaussian mixed model + stage-2 LMM refit
model3_result <- tryCatch({
  model3_fit <- fit_glmmlasso_with_bic(
    formula = model3_formula,
    data = pheno_file_feed,
    family = gaussian(link = "identity"),
    lambda_grid = lambda_grid
  )
  optimal_lambda_model3 <- model3_fit$optimal_lambda

  res <- extract_glmmlasso_results(
    fit_obj = model3_fit$fit,
    features = features_model3,
    coef_names = coef_names_model3
  )

  selected_features <- names(res$coefficients)[res$coefficients != 0]
  selected_features <- setdiff(selected_features, "(Intercept)")

  if (length(selected_features) > 0) {
    stage2_formula <- as.formula(
      paste0("meth_response ~ ", paste(selected_features, collapse = " + "), " + (1|FID)")
    )

    stage2_res <- tryCatch({
      lmod <- lmer(
        formula = stage2_formula,
        data = pheno_file_feed,
        REML = FALSE
      )

      smry <- summary(lmod)$coefficients
      pvals <- smry[, "Pr(>|t|)"]
      ests <- smry[, "Estimate"]

      selected_in_summary <- intersect(selected_features, rownames(smry))
      res$coef_pvals[selected_in_summary] <- pvals[selected_in_summary]
      res$coefficients[selected_in_summary] <- ests[selected_in_summary]
      res$FID_sdv <- as.numeric(attr(VarCorr(lmod)$FID, "stddev"))

      res
    }, error = function(e) res)

    stage2_res
  } else {
    res
  }
}, error = function(e) {
  warning("Model 3 failed for CpG ", i_cpg, ": ", e$message)
  default_model_result(coef_names_model3, features_model3)
})

model1_coefficients <- model1_result$coefficients
model1_FID_sdv <- model1_result$FID_sdv
model1_coef_pvals <- model1_result$coef_pvals

model2_coefficients <- model2_result$coefficients
model2_FID_sdv <- model2_result$FID_sdv
model2_coef_pvals <- model2_result$coef_pvals

model3_coefficients <- model3_result$coefficients
model3_FID_sdv <- model3_result$FID_sdv
model3_coef_pvals <- model3_result$coef_pvals

output_complete_results <- c(
  i_cpg, meth_cpg_info,
  model1_coef_pvals, model1_coefficients, model1_FID_sdv,
  model2_coef_pvals, model2_coefficients, model2_FID_sdv,
  model3_coef_pvals, model3_coefficients, model3_FID_sdv
)

output_simple_results <- c(
  i_cpg,
  meth_cpg_info,
  model1_coefficients["meth_response"],
  model1_coef_pvals["meth_response"],
  model1_FID_sdv,
  model2_coefficients["AA_only"],
  model2_coef_pvals["AA_only"],
  model2_FID_sdv,
  model3_coefficients["AA_only"],
  model3_coef_pvals["AA_only"],
  model3_FID_sdv,
  optimal_lambda_model1,
  optimal_lambda_model3
)

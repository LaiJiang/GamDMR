
# This script runs an MGCV model on a single region of interest from one data chunk
#adding effect size estimation to the 6_run_region.R script 

i_region <- 57
min_cpgs <- 10

# Define paths
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
PATH_results <- paste0(PATH_wk, "results/10_sanity/")
PATH_scr6 <- paste0(PATH_wk, "scr/6_beluga/")
PATH_scr11 <- paste0(PATH_wk, "scr/11_mgcv/")
PATH_output <- paste0(PATH_wk, "results/11_mgcv/")

# Ensure output directory exists
if (!dir.exists(PATH_output)) dir.create(PATH_output, recursive = TRUE)

# Load libraries
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(tidyr)
library(mgcv)

# Load region info
region_file <- fread(paste0(PATH_scr11, "dat/region_file_1_chunk.csv"))
i_region_info <- region_file[i_region, ]
i_chunk_id <- i_region_info$data_chunk_id

# Load the corresponding data chunk
i_chunk <- fread(paste0(PATH_wk, "dat/chunk_", sprintf("%04d", i_chunk_id), ".csv"))

# Extract region of interest
i_region_data <- i_chunk %>%
  filter(chr == i_region_info$chr,
         start >= i_region_info$region_start,
         start <= i_region_info$region_end)

if (nrow(i_region_data) > min_cpgs) {
  setDT(i_region_data)

  # Convert _meth columns to numeric
  meth_cols <- grep("_meth$", names(i_region_data), value = TRUE)
  i_region_data[, (meth_cols) := lapply(.SD, as.numeric), .SDcols = meth_cols]

  # Melt to long format
  meth_long <- melt(i_region_data,
                    id.vars = "start",
                    measure.vars = patterns("_meth$"),
                    variable.name = "sample",
                    value.name = "meth",
                    na.rm = TRUE)
  meth_long[, sample := str_replace(as.character(sample), "_meth", "")]

  # Load phenotype data
  load(file = paste0(PATH_scr6, "4_1_data.RData"), verbose = TRUE)
  pheno_file <- as_tibble(pheno_file)

  # Join phenotype data
  meth_long <- meth_long %>%
    left_join(pheno_file %>% select(ID, AA_only), by = c("sample" = "ID")) %>%
    na.omit()

  i_pheno_file <- pheno_file %>% rename(sample = ID) %>% select(-AA_only)

  i_meth_long_cov <- meth_long %>%
    left_join(i_pheno_file, by = "sample") %>%
    na.omit()

  if (nrow(i_meth_long_cov) < min_cpgs) {
    stop("Too few CpGs after full data join to proceed.")
  }

  # Arcsin transform
  i_meth_long_cov$meth_arcsin <- asin(sqrt(i_meth_long_cov$meth))

  # Fix type bugs
  for (v in colnames(i_meth_long_cov)) {
    if (is.matrix(i_meth_long_cov[[v]]) || is.array(i_meth_long_cov[[v]])) {
      i_meth_long_cov[[v]] <- as.numeric(i_meth_long_cov[[v]])
    }
  }

  # Define basis dimension
  k_cpg <- max(4, min(10, floor(length(unique(i_meth_long_cov$start)) / 20)))

  # Fit GAM model
  gam_model <- mgcv::gam(meth_arcsin ~ 
    s(start, bs = "cs", k = k_cpg) +
    s(start, by = AA_only, bs = "cs", k = k_cpg) +
    AA_only + AgeCalc + Sex + Non.smoker +
    EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc +
    sv1 + sv2 + sv3 + sv4 + sv5 +
    s(FID, bs = "re"),
    data = i_meth_long_cov,
    method = "REML")

  # Extract results
  feature_list <- c("(Intercept)", "AA_only", "AgeCalc", "Sex", "Non.smoker",
                    "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc", 
                    "sv1", "sv2", "sv3", "sv4", "sv5")
  feature_smooth_terms <- c("s(start)", "s(start):AA_only", "s(FID)")

  param_coef <- summary(gam_model)$p.table
  smooth_terms <- summary(gam_model)$s.table



  present_features <- intersect(feature_list, rownames(param_coef))
  result_pvalue_fixed_effect <- rep(NA, length(feature_list))
  names(result_pvalue_fixed_effect) <- feature_list
  result_pvalue_fixed_effect[present_features] <- param_coef[present_features, "Pr(>|t|)"]

  present_smooth <- intersect(feature_smooth_terms, rownames(smooth_terms))
  result_smooth_terms <- rep(NA, length(feature_smooth_terms))
  names(result_smooth_terms) <- feature_smooth_terms
  result_smooth_terms[present_smooth] <- smooth_terms[present_smooth, "p-value"]



  #now get the effect size estimations
  # Create new data with AA_only = 1 and 0
newdata_AA1 <- i_meth_long_cov
newdata_AA1$AA_only <- 1
newdata_AA0 <- i_meth_long_cov
newdata_AA0$AA_only <- 0

# Predict
pred_AA1 <- predict(gam_model, newdata = newdata_AA1, type = "response")
pred_AA0 <- predict(gam_model, newdata = newdata_AA0, type = "response")

# Effect size
delta <- pred_AA1 - pred_AA0
max_diff <- max(abs(delta))
mean_diff <- mean(abs(delta))

  # Model fit statistics
  result_stats <- c(result_pvalue_fixed_effect,
                    result_smooth_terms,
                    "R2" = summary(gam_model)$r.sq,
                    "AIC" = AIC(gam_model),
                    "Deviance_explained" = summary(gam_model)$dev.expl,
                    "REML" = gam_model$gcv.ubre,
                    "N_cpgs" = length(unique(i_meth_long_cov$start)),
                    "N_samples" = length(unique(i_meth_long_cov$sample)),
                    "max_diff" = max_diff,
                    "mean_diff" = mean_diff)

  # Save output
  result_stats_df <- data.frame(term = names(result_stats), value = as.numeric(result_stats))
  fwrite(result_stats_df, file = paste0(PATH_output, "pval_region_", i_region, ".csv"))
}

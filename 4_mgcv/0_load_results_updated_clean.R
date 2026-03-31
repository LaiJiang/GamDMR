suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

load_mgcv_batch <- function(path_dat, include_bmi = TRUE) {
  all_files <- list.files(path_dat, pattern = "\\.txt$", full.names = TRUE)
  if (length(all_files) == 0) return(data.table())

  feature_fixed <- c(
    "(Intercept)", "AA_only", "AgeCalc", "Sex", "Non.smoker",
    "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc",
    "sv1", "sv2", "sv3", "sv4", "sv5"  )
    
  if (include_bmi) feature_fixed <- c(feature_fixed, "BMI")

  feature_smooth_terms <- c("s(start)", "s(start):AA_only", "s(FID)")
  feature_edf <- paste0("edf_", feature_smooth_terms)

  features_column_names <- c(
    feature_fixed,
    feature_smooth_terms,
    feature_edf,
    "R2",
    "AIC",
    "Deviance_explained",
    "REML",
    "N_cpgs",
    "N_samples",
    "data_chunk_id",
    "region_id",
    "max_diff",
    "mean_diff"
  )

  results_col <- rbindlist(lapply(all_files, function(i_file) {
    i_file_dat <- fread(i_file, header = TRUE, sep = "\t")
    if (ncol(i_file_dat) != length(features_column_names)) {
      stop("Unexpected number of columns in file: ", i_file)
    }
    setnames(i_file_dat, features_column_names)
    i_file_dat$job_id <- as.numeric(str_extract(i_file, "(?<=_job_)\\d+(?=\\.txt)"))
    i_file_dat
  }), fill = TRUE)

  results_col
}

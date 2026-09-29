
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

load_mgcv_batch <- function(path_dat, include_bmi = TRUE) {
  all_files <- list.files(path_dat, pattern = "\\.txt$", full.names = TRUE)
  if (length(all_files) == 0L) return(data.table())

  fixed <- c(
    "(Intercept)", "AA_only", "AgeCalc", "Sex", "Non.smoker",
    "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc",
    "sv1", "sv2", "sv3", "sv4", "sv5"
  )
  if (include_bmi) fixed <- c(fixed, "BMI")

  smooth <- c("s(start)", "s(start):AA_only", "s(FID)")
  edf <- paste0("edf_", smooth)
  cols <- c(
    fixed, smooth, edf, "R2", "AIC", "Deviance_explained", "REML",
    "N_cpgs", "N_samples", "data_chunk_id", "region_id", "max_diff", "mean_diff"
  )

  rbindlist(lapply(all_files, function(f) {
    x <- fread(f, header = TRUE, sep = "\t")
    if (ncol(x) != length(cols)) stop("Unexpected column count: ", f)
    setnames(x, cols)
    x[, job_id := as.integer(str_extract(basename(f), "(?<=_job_)\\d+"))]
    x
  }), use.names = TRUE, fill = TRUE)
}

select_strict_gam_dmrs <- function(x, r2_threshold = 0.50) {
  stopifnot("s(start):AA_only" %in% names(x))
  stopifnot("edf_s(start):AA_only" %in% names(x))
  x <- copy(as.data.table(x))
  x[, FDR := p.adjust(`s(start):AA_only`, method = "BH")]
  x[
    `s(start):AA_only` < 1e-5 &
      FDR < 0.05 &
      `edf_s(start):AA_only` >= 0.5 &
      R2 >= r2_threshold &
      mean_diff >= 0.05 &
      N_cpgs >= 10
  ]
}

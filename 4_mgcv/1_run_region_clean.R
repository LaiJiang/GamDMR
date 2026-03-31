suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(mgcv)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 1_run_region_clean.R <job_id>")
}
job_id <- as.integer(args[1])
if (is.na(job_id)) {
  stop("job_id must be an integer.")
}

min_cpgs <- 10
n_jobs <- 300

base_dir <- Sys.getenv("METH_BASE_DIR", unset = getwd())
data_dir <- file.path(base_dir, "data", "meth_split")
script_dir <- file.path(base_dir, "4_mgcv")
output_dir <- file.path(base_dir, "results", "4_mgcv", "B1")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

region_file_path <- file.path(script_dir, "dat", "region_file_1_chunk.csv")
pheno_file_path <- file.path(script_dir, "dat", "18_pheno_BMI.RData")

if (!file.exists(region_file_path)) stop("Missing region file: ", region_file_path)
if (!file.exists(pheno_file_path)) stop("Missing phenotype file: ", pheno_file_path)

region_file <- fread(region_file_path)
region_index_list <- split(
  seq_len(nrow(region_file)),
  cut(seq_len(nrow(region_file)), breaks = n_jobs, labels = FALSE)
)

if (job_id < 1 || job_id > length(region_index_list) || is.null(region_index_list[[job_id]])) {
  stop("job_id is outside the valid range.")
}

load(pheno_file_path)
pheno_file <- as_tibble(pheno_file)

region_ids <- region_index_list[[job_id]]
output_file <- file.path(output_dir, sprintf("results_job_%03d.txt", job_id))

for (i_region in region_ids) {
  i_region_info <- region_file[i_region, ]
  i_chunk_id <- as.integer(i_region_info$data_chunk_id)

  chunk_path <- file.path(data_dir, sprintf("chunk_%04d.csv", i_chunk_id))
  if (!file.exists(chunk_path)) {
    warning("Skipping region ", i_region, ": missing chunk file ", chunk_path)
    next
  }

  i_chunk <- fread(chunk_path)

  i_region_data <- i_chunk %>%
    filter(
      chr == i_region_info$chr,
      start >= i_region_info$region_start,
      start <= i_region_info$region_end
    )

  if (nrow(i_region_data) < min_cpgs) {
    next
  }

  setDT(i_region_data)
  meth_cols <- grep("_meth$", names(i_region_data), value = TRUE)
  if (length(meth_cols) == 0) {
    warning("Skipping region ", i_region, ": no methylation columns found.")
    next
  }

  i_region_data[, (meth_cols) := lapply(.SD, as.numeric), .SDcols = meth_cols]

  meth_long <- melt(
    i_region_data,
    id.vars = "start",
    measure.vars = patterns("_meth$"),
    variable.name = "sample",
    value.name = "meth",
    na.rm = TRUE
  )
  meth_long[, sample := str_replace(as.character(sample), "_meth", "")]

  meth_long <- meth_long %>%
    left_join(pheno_file %>% select(ID, AA_only), by = c("sample" = "ID")) %>%
    na.omit()

  covariates <- pheno_file %>%
    rename(sample = ID) %>%
    select(-AA_only)

  i_meth_long_cov <- meth_long %>%
    left_join(covariates, by = "sample") %>%
    na.omit()

  if (nrow(i_meth_long_cov) < min_cpgs) {
    next
  }

  i_meth_long_cov$meth_arcsin <- asin(sqrt(i_meth_long_cov$meth))

  for (v in names(i_meth_long_cov)) {
    if (is.matrix(i_meth_long_cov[[v]]) || is.array(i_meth_long_cov[[v]])) {
      i_meth_long_cov[[v]] <- as.numeric(i_meth_long_cov[[v]])
    }
  }

  k_cpg <- max(5, min(10, floor(length(unique(i_meth_long_cov$start)) / 20)))

  gam_model <- tryCatch(
    gam(
      meth_arcsin ~
        s(start, bs = "cs", k = k_cpg) +
        s(start, by = AA_only, bs = "cs", k = k_cpg) +
        AA_only + AgeCalc + Sex + Non.smoker +
        EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc +
        sv1 + sv2 + sv3 + sv4 + sv5 + BMI +
        s(FID, bs = "re"),
      data = i_meth_long_cov,
      method = "REML"
    ),
    error = function(e) NULL
  )

  if (is.null(gam_model)) {
    warning("Skipping region ", i_region, ": model fit failed.")
    next
  }

  feature_list <- c(
    "(Intercept)", "AA_only", "AgeCalc", "Sex", "Non.smoker",
    "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc",
    "sv1", "sv2", "sv3", "sv4", "sv5", "BMI"
  )
  smooth_terms <- c("s(start)", "s(start):AA_only", "s(FID)")

  gam_summary <- summary(gam_model)
  param_coef <- gam_summary$p.table
  smooth_table <- gam_summary$s.table

  result_pvalue_fixed_effect <- setNames(rep(1, length(feature_list)), feature_list)
  present_features <- intersect(feature_list, rownames(param_coef))
  result_pvalue_fixed_effect[present_features] <- param_coef[present_features, "Pr(>|t|)"]

  result_smooth_terms <- setNames(rep(1, length(smooth_terms)), smooth_terms)
  present_smooth <- intersect(smooth_terms, rownames(smooth_table))
  result_smooth_terms[present_smooth] <- smooth_table[present_smooth, "p-value"]

  edf_smooth_terms <- setNames(rep(0, length(smooth_terms)), paste0("edf_", smooth_terms))
  edf_smooth_terms[paste0("edf_", present_smooth)] <- smooth_table[present_smooth, "edf"]

  newdata_AA1 <- i_meth_long_cov
  newdata_AA1$AA_only <- 1
  newdata_AA0 <- i_meth_long_cov
  newdata_AA0$AA_only <- 0

  pred_AA1 <- predict(gam_model, newdata = newdata_AA1, type = "response")
  pred_AA0 <- predict(gam_model, newdata = newdata_AA0, type = "response")
  delta <- pred_AA1 - pred_AA0

  result_stats <- c(
    result_pvalue_fixed_effect,
    result_smooth_terms,
    edf_smooth_terms,
    R2 = gam_summary$r.sq,
    AIC = AIC(gam_model),
    Deviance_explained = gam_summary$dev.expl,
    REML = gam_model$gcv.ubre,
    N_cpgs = length(unique(i_meth_long_cov$start)),
    N_samples = length(unique(i_meth_long_cov$sample)),
    data_chunk_id = i_chunk_id,
    region_id = i_region,
    max_diff = max(abs(delta)),
    mean_diff = mean(abs(delta))
  )

  fwrite(
    as.data.frame(as.list(result_stats)),
    file = output_file,
    append = TRUE,
    sep = "\t",
    col.names = !file.exists(output_file)
  )
}


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(mgcv)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) stop("Usage: Rscript 1_run_region_clean.R <job_id>")
job_id <- as.integer(args[1])
if (is.na(job_id)) stop("job_id must be an integer.")

min_cpgs <- 10L
n_jobs <- 300L
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
if (!"region_id" %in% names(region_file)) region_file[, region_id := .I]
region_index_list <- split(
  seq_len(nrow(region_file)),
  cut(seq_len(nrow(region_file)), breaks = n_jobs, labels = FALSE)
)
if (job_id < 1L || job_id > length(region_index_list) || is.null(region_index_list[[job_id]])) {
  stop("job_id is outside the valid range.")
}

env <- new.env(parent = emptyenv())
load(pheno_file_path, envir = env)
if (!exists("pheno_file", envir = env, inherits = FALSE)) stop("pheno_file is missing.")
pheno_file <- as_tibble(get("pheno_file", envir = env))
rm(env)
pheno_file <- pheno_file %>%
  mutate(ID = as.character(ID), FID = as.character(FID))

region_ids <- region_index_list[[job_id]]
output_file <- file.path(output_dir, sprintf("results_job_%03d.txt", job_id))

for (i_region in region_ids) {
  info <- region_file[i_region]
  chunk_path <- file.path(data_dir, sprintf("chunk_%04d.csv", as.integer(info$data_chunk_id)))
  if (!file.exists(chunk_path)) next

  region_dt <- fread(chunk_path)[
    chr == info$chr &
      start >= info$region_start &
      start <= info$region_end
  ]
  if (uniqueN(region_dt$start) < min_cpgs) next

  meth_cols <- grep("_meth$", names(region_dt), value = TRUE)
  if (length(meth_cols) == 0L) next
  region_dt[, (meth_cols) := lapply(.SD, as.numeric), .SDcols = meth_cols]

  meth_long <- melt(
    region_dt,
    id.vars = "start",
    measure.vars = meth_cols,
    variable.name = "sample",
    value.name = "meth",
    na.rm = TRUE
  )
  meth_long[, sample := str_remove(as.character(sample), "_meth$")]

  dat <- meth_long %>%
    inner_join(pheno_file %>% rename(sample = ID), by = "sample")

  model_vars <- c(
    "meth", "start", "FID", "AA_only", "AgeCalc", "Sex", "Non.smoker",
    "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc",
    "sv1", "sv2", "sv3", "sv4", "sv5", "BMI"
  )
  dat <- dat[complete.cases(dat[, model_vars]), , drop = FALSE]
  if (nrow(dat) == 0L) next

  for (v in names(dat)) {
    if (is.matrix(dat[[v]]) || is.array(dat[[v]])) dat[[v]] <- as.numeric(dat[[v]])
  }

  dat$FID <- droplevels(factor(dat$FID))
  dat$AA_only <- as.numeric(dat$AA_only)
  n_cpgs <- data.table::uniqueN(dat$start)
  if (n_cpgs < min_cpgs || length(unique(dat$AA_only)) < 2L || nlevels(dat$FID) < 2L) next

  dat$meth_arcsin <- asin(sqrt(dat$meth))
  k_cpg <- max(5L, min(10L, floor(n_cpgs / 20L)))

  fit <- tryCatch(
    mgcv::gam(
      meth_arcsin ~
        s(start, bs = "cs", k = k_cpg) +
        s(start, by = AA_only, bs = "cs", k = k_cpg) +
        AA_only + AgeCalc + Sex + Non.smoker +
        EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc +
        sv1 + sv2 + sv3 + sv4 + sv5 + BMI +
        s(FID, bs = "re"),
      data = dat,
      method = "REML"
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) next

  sm <- summary(fit)
  feature_list <- c(
    "(Intercept)", "AA_only", "AgeCalc", "Sex", "Non.smoker",
    "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc",
    "sv1", "sv2", "sv3", "sv4", "sv5", "BMI"
  )
  smooth_names <- c("s(start)", "s(start):AA_only", "s(FID)")
  p_fixed <- setNames(rep(1, length(feature_list)), feature_list)
  p_smooth <- setNames(rep(1, length(smooth_names)), smooth_names)
  edf_smooth <- setNames(rep(0, length(smooth_names)), paste0("edf_", smooth_names))

  present_fixed <- intersect(feature_list, rownames(sm$p.table))
  present_smooth <- intersect(smooth_names, rownames(sm$s.table))
  p_fixed[present_fixed] <- sm$p.table[present_fixed, "Pr(>|t|)"]
  p_smooth[present_smooth] <- sm$s.table[present_smooth, "p-value"]
  edf_smooth[paste0("edf_", present_smooth)] <- sm$s.table[present_smooth, "edf"]

  d1 <- dat
  d0 <- dat
  d1$AA_only <- 1
  d0$AA_only <- 0
  delta <- predict(fit, newdata = d1, type = "response") -
    predict(fit, newdata = d0, type = "response")

  out <- c(
    p_fixed,
    p_smooth,
    edf_smooth,
    R2 = sm$r.sq,
    AIC = AIC(fit),
    Deviance_explained = sm$dev.expl,
    REML = fit$gcv.ubre,
    N_cpgs = n_cpgs,
    N_samples = data.table::uniqueN(dat$sample),
    data_chunk_id = as.integer(info$data_chunk_id),
    region_id = as.integer(info$region_id),
    max_diff = max(abs(delta), na.rm = TRUE),
    mean_diff = mean(abs(delta), na.rm = TRUE)
  )

  fwrite(
    as.data.frame(as.list(out)),
    output_file,
    append = file.exists(output_file),
    sep = "\t",
    col.names = !file.exists(output_file)
  )
}

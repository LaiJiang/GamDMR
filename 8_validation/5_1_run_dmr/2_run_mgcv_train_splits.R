
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(mgcv)
})

min_cpgs <- 10L
N_jobs <- as.integer(Sys.getenv("N_JOBS", "990"))
PATH_wk <- path.expand(Sys.getenv("METH_BASE_DIR", "~/scratch/UQAC/meth/"))
PATH_scr11 <- file.path(PATH_wk, "scr", "11_mgcv")
PATH_output <- file.path(PATH_wk, "results", "15_revision", "5_cv", "1_mgcv")
PATH_train_splits <- file.path(PATH_scr11, "dat", "100_family_training_splits.csv")
dir.create(PATH_output, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) stop("Pass SLURM_ARRAY_TASK_ID.")
job_id <- as.integer(args[1])
if (is.na(job_id)) stop("Invalid job ID.")

region_file <- fread(file.path(PATH_scr11, "dat", "region_file_1_chunk.csv"))
if (!"region_id" %in% names(region_file)) region_file[, region_id := .I]
if ("n_cpgs" %in% names(region_file)) region_file <- region_file[n_cpgs >= min_cpgs]

train_df <- fread(PATH_train_splits)
train_cols <- grep("^train_FID_[0-9]+$", names(train_df), value = TRUE)
train_cols <- train_cols[order(as.integer(sub("^train_FID_", "", train_cols)))]

env <- new.env(parent = emptyenv())
load(file.path(PATH_scr11, "dat", "18_pheno_BMI.RData"), envir = env)
if (!exists("pheno_file", envir = env, inherits = FALSE)) stop("pheno_file is missing.")
pheno <- as.data.table(get("pheno_file", envir = env))
rm(env)
pheno[, `:=`(ID = as.character(ID), FID = as.character(FID))]

job_groups <- split(
  seq_len(nrow(region_file)),
  cut(seq_len(nrow(region_file)), breaks = min(N_jobs, nrow(region_file)), labels = FALSE)
)
if (job_id < 1L || job_id > length(job_groups)) quit(save = "no", status = 0L)
rows <- job_groups[[job_id]]
output_file <- file.path(PATH_output, sprintf("2_run_mgcv_training_splits_results_%d.txt", job_id))

feature_list <- c(
  "(Intercept)", "AA_only", "AgeCalc", "Sex", "Non.smoker",
  "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc",
  "sv1", "sv2", "sv3", "sv4", "sv5", "BMI"
)
smooth_names <- c("s(start)", "s(start):AA_only", "s(FID)")

for (ri in rows) {
  info <- region_file[ri]
  chunk_file <- file.path(PATH_wk, "data", "meth_split", sprintf("chunk_%04d.csv", info$data_chunk_id))
  if (!file.exists(chunk_file)) next

  region_dt <- fread(chunk_file)[
    chr == info$chr &
      start >= info$region_start &
      start <= info$region_end
  ]
  if (uniqueN(region_dt$start) < min_cpgs) next

  all_meth_cols <- grep("_meth$", names(region_dt), value = TRUE)
  if (!length(all_meth_cols)) next

  for (sid in train_df$splitID) {
    raw_fids <- unlist(train_df[splitID == sid, ..train_cols], use.names = FALSE)
    train_fids <- unique(as.character(raw_fids))
    train_fids <- train_fids[!is.na(train_fids) & nzchar(train_fids)]

    ph_train <- pheno[FID %in% train_fids]
    train_ids <- ph_train$ID
    meth_cols <- intersect(paste0(train_ids, "_meth"), all_meth_cols)
    if (!length(meth_cols)) next

    x <- copy(region_dt)
    x[, (meth_cols) := lapply(.SD, as.numeric), .SDcols = meth_cols]
    long <- melt(
      x,
      id.vars = "start",
      measure.vars = meth_cols,
      variable.name = "sample",
      value.name = "meth",
      na.rm = TRUE
    )
    long[, sample := str_remove(as.character(sample), "_meth$")]
    dat <- merge(long, ph_train, by.x = "sample", by.y = "ID", all = FALSE)

    model_vars <- c(
      "meth", "start", "FID", "AA_only", "AgeCalc", "Sex", "Non.smoker",
      "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc",
      "sv1", "sv2", "sv3", "sv4", "sv5", "BMI"
    )
    dat <- dat[complete.cases(dat[, ..model_vars])]
    if (!nrow(dat)) next

    for (v in names(dat)) {
      if (is.matrix(dat[[v]]) || is.array(dat[[v]])) dat[[v]] <- as.numeric(dat[[v]])
    }
    dat[, FID := droplevels(factor(FID))]
    dat[, AA_only := as.numeric(AA_only)]

    n_cpgs <- uniqueN(dat$start)
    n_samples <- uniqueN(dat$sample)
    n_fids <- nlevels(dat$FID)
    if (n_cpgs < min_cpgs || uniqueN(dat$AA_only) < 2L || n_fids < 2L) next

    dat[, meth_arcsin := asin(sqrt(meth))]
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
    p_fixed <- setNames(rep(1, length(feature_list)), feature_list)
    p_smooth <- setNames(rep(1, length(smooth_names)), smooth_names)
    edf_smooth <- setNames(rep(0, length(smooth_names)), smooth_names)

    present_fixed <- intersect(feature_list, rownames(sm$p.table))
    present_smooth <- intersect(smooth_names, rownames(sm$s.table))
    p_fixed[present_fixed] <- sm$p.table[present_fixed, "Pr(>|t|)"]
    p_smooth[present_smooth] <- sm$s.table[present_smooth, "p-value"]
    edf_smooth[present_smooth] <- sm$s.table[present_smooth, "edf"]

    d1 <- copy(dat)
    d0 <- copy(dat)
    d1[, AA_only := 1]
    d0[, AA_only := 0]
    delta <- predict(fit, newdata = d1, type = "response") -
      predict(fit, newdata = d0, type = "response")

    row <- c(
      splitID = sid,
      p_fixed,
      p_smooth,
      edf_smooth,
      R2 = sm$r.sq,
      AIC = AIC(fit),
      Deviance_explained = sm$dev.expl,
      REML = fit$gcv.ubre,
      N_cpgs = n_cpgs,
      N_samples = n_samples,
      N_train_FIDs = n_fids,
      data_chunk_id = as.integer(info$data_chunk_id),
      region_id = as.integer(info$region_id),
      max_diff = max(abs(delta), na.rm = TRUE),
      mean_diff = mean(abs(delta), na.rm = TRUE)
    )

    fwrite(
      as.data.frame(t(as.numeric(row))),
      output_file,
      append = file.exists(output_file),
      sep = "\t",
      col.names = FALSE
    )
  }
}

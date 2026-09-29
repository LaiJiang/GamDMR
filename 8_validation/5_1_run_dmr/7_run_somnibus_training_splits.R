
#!/usr/bin/env Rscript

user_lib <- path.expand("~/scratch/R/library")
.libPaths(unique(c(user_lib, .libPaths())))

suppressPackageStartupMessages({
  library(data.table)
  library(SOMNiBUS)
})

PATH_wk <- path.expand(Sys.getenv("METH_BASE_DIR", "~/scratch/UQAC/meth/"))
bundle_dir <- file.path(PATH_wk, "results", "15_revision", "5_cv", "2_somnibus", "prep", "bundles")
train_split_file <- file.path(PATH_wk, "scr", "11_mgcv", "dat", "100_family_training_splits.csv")
pheno_path <- file.path(PATH_wk, "scr", "11_mgcv", "dat", "18_pheno_BMI.RData")
output_dir <- file.path(PATH_wk, "results", "15_revision", "5_cv", "2_somnibus", "training_splits")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

min_cpgs <- 50L
max_cpgs <- 2000L
p0_value <- 0.003
p1_value <- 0.9
epsilon_value <- 1e-6
epsilon_lambda_value <- 0.001
max_step <- 200L

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Pass SLURM_ARRAY_TASK_ID.")
job_id <- as.integer(args[1])
if (is.na(job_id)) stop("Invalid job ID.")

bundle_file <- file.path(bundle_dir, sprintf("somnibus_bundle_job_%04d.rds", job_id))
if (!file.exists(bundle_file)) stop("Missing bundle: ", bundle_file)
bundle <- readRDS(bundle_file)
if (!all(c("sample_map", "manifest_ready", "regions") %in% names(bundle))) stop("Invalid bundle.")

env <- new.env(parent = emptyenv())
load(pheno_path, envir = env)
if (!exists("pheno_file", envir = env, inherits = FALSE)) stop("pheno_file is missing.")
pheno <- as.data.table(get("pheno_file", envir = env))
rm(env)

covs <- c(
  "AA_only", "AgeCalc", "Sex", "Non.smoker", "BMI",
  "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc",
  "sv1", "sv2", "sv3", "sv4", "sv5"
)
need_pheno <- c("ID", "FID", covs)
missing <- setdiff(need_pheno, names(pheno))
if (length(missing)) stop("Missing phenotype columns: ", paste(missing, collapse = ", "))
pheno <- pheno[, ..need_pheno]
pheno[, `:=`(ID = as.character(ID), FID = as.character(FID))]

train_df <- fread(train_split_file)
train_cols <- grep("^train_FID_[0-9]+$", names(train_df), value = TRUE)
train_cols <- train_cols[order(as.integer(sub("^train_FID_", "", train_cols)))]

manifest <- as.data.table(copy(bundle$manifest_ready))
if (nrow(manifest)) {
  manifest[, omnibus_region_id := as.character(omnibus_region_id)]
  setkey(manifest, omnibus_region_id)
}

result_file <- file.path(output_dir, sprintf("somnibus_training_splits_job_%04d.tsv", job_id))
rows <- list()

for (region_name in names(bundle$regions)) {
  region_all <- as.data.table(copy(bundle$regions[[region_name]]))
  required <- c("Meth_Counts", "Total_Counts", "Position", "ID")
  if (!all(required %in% names(region_all))) next
  region_all <- region_all[, ..required]
  region_all[, ID := as.character(ID)]
  region_all[, `:=`(
    Meth_Counts = as.numeric(Meth_Counts),
    Total_Counts = as.numeric(Total_Counts),
    Position = as.integer(Position)
  )]
  region_all <- region_all[
    is.finite(Meth_Counts) & is.finite(Total_Counts) &
      Total_Counts != 0 & Meth_Counts >= 0 & Meth_Counts <= Total_Counts
  ]

  meta <- if (nrow(manifest)) manifest[J(region_name), nomatch = 0L] else data.table()
  if (!nrow(meta)) {
    meta <- data.table(
      omnibus_region_id = region_name,
      data_chunk_id = NA_integer_,
      parent_region_id = NA_integer_,
      chr = NA_character_,
      region_start = NA_integer_,
      region_end = NA_integer_
    )
  }

  for (sid in train_df$splitID) {
    raw_fids <- unlist(train_df[splitID == sid, ..train_cols], use.names = FALSE)
    train_fids <- unique(as.character(raw_fids))
    train_fids <- train_fids[!is.na(train_fids) & nzchar(train_fids)]

    ph <- pheno[FID %in% train_fids]
    ph <- ph[complete.cases(ph[, ..need_pheno])]
    dat <- merge(region_all, ph, by = "ID", all = FALSE)
    if (!nrow(dat)) next

    n_cpg <- uniqueN(dat$Position)
    n_samples <- uniqueN(dat$ID)
    n_fids <- uniqueN(dat$FID)
    n_case <- uniqueN(dat[AA_only == 1L, ID])
    n_ctrl <- uniqueN(dat[AA_only == 0L, ID])

    if (n_cpg < min_cpgs || n_cpg > max_cpgs || n_samples < 2L ||
        n_fids < 2L || n_case == 0L || n_ctrl == 0L) next

    model_dat <- as.data.frame(dat[, c(
      "Meth_Counts", "Total_Counts", "Position", "ID", covs
    ), with = FALSE])
    k <- max(5L, min(10L, floor(n_cpg / 20L)))

    fit <- tryCatch(
      SOMNiBUS::binomRegMethModel(
        data = model_dat,
        n.k = rep(k, length(covs) + 1L),
        p0 = p0_value,
        p1 = p1_value,
        Quasi = TRUE,
        epsilon = epsilon_value,
        epsilon.lambda = epsilon_lambda_value,
        maxStep = max_step,
        binom.link = "logit",
        method = "REML",
        covs = covs,
        RanEff = TRUE,
        reml.scale = -2,
        scale = -2,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(fit) || is.null(fit$reg.out) || !"AA_only" %in% rownames(fit$reg.out)) next

    stat_col <- intersect(c("Chi.sq", "F"), colnames(fit$reg.out))
    rows[[length(rows) + 1L]] <- data.table(
      prep_job_id = job_id,
      omnibus_region_id = as.character(region_name),
      data_chunk_id = as.integer(meta$data_chunk_id[1L]),
      parent_region_id = as.integer(meta$parent_region_id[1L]),
      chr = as.character(meta$chr[1L]),
      region_start = as.integer(meta$region_start[1L]),
      region_end = as.integer(meta$region_end[1L]),
      splitID = as.integer(sid),
      status = "FIT_OK",
      n_train_FIDs_requested = length(train_fids),
      n_train_FIDs_observed = n_fids,
      n_train_samples = n_samples,
      n_train_cases = n_case,
      n_train_controls = n_ctrl,
      n_train_cpgs = n_cpg,
      n_train_rows = nrow(dat),
      n_k = k,
      aa_EDF = if ("EDF" %in% colnames(fit$reg.out)) as.numeric(fit$reg.out["AA_only", "EDF"]) else NA_real_,
      aa_statistic = if (length(stat_col)) as.numeric(fit$reg.out["AA_only", stat_col[1L]]) else NA_real_,
      aa_statistic_type = if (length(stat_col)) stat_col[1L] else NA_character_,
      aa_pvalue = if ("p-value" %in% colnames(fit$reg.out)) as.numeric(fit$reg.out["AA_only", "p-value"]) else NA_real_
    )
  }
}

if (length(rows)) {
  fwrite(rbindlist(rows, fill = TRUE), result_file, sep = "\t")
} else {
  fwrite(data.table(
    prep_job_id = integer(), omnibus_region_id = character(),
    data_chunk_id = integer(), parent_region_id = integer(), chr = character(),
    region_start = integer(), region_end = integer(), splitID = integer(),
    status = character(), aa_pvalue = numeric()
  ), result_file, sep = "\t")
}

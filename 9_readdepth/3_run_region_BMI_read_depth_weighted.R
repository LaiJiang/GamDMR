#!/usr/bin/env Rscript

# =============================================================================
# 22_run_region_BMI_read_depth_weighted.R
#
# Read-depth-weighted sensitivity analysis for the existing GAM-DMR model.
#
# This script mirrors 20_run_region_BMI.R as closely as possible, with one
# deliberate change:
#
#   mgcv::gam(..., weights = depth_weight, method = "REML")
#
# where depth_weight is proportional to total read depth and normalized to have
# mean 1 within each fitted region:
#
#   depth_weight = read_depth / mean(read_depth)
#
# For ASR-transformed binomial proportions,
# Var[asin(sqrt(p_hat))] is approximately proportional to 1 / read_depth, so a
# weight proportional to read depth is a natural inverse-variance sensitivity
# analysis. Normalizing the weights preserves relative precision while avoiding
# an arbitrary change in the overall weight scale.
#
# Input format is the same as the original script:
#   ~/scratch/UQAC/meth/data/meth_split/chunk_XXXX.csv
#
# Expected paired columns:
#   <sample>_meth
#   <sample>_tot
#
# Primary model:
#   meth_arcsin ~
#     s(start, bs="cs", k=k_cpg) +
#     s(start, by=AA_only, bs="cs", k=k_cpg) +
#     AA_only + AgeCalc + Sex + Non.smoker +
#     EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc +
#     sv1 + sv2 + sv3 + sv4 + sv5 + BMI +
#     s(FID, bs="re")
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(mgcv)
})

data.table::setDTthreads(1L)
options(mc.cores = 1L)

# ------------------------------- configuration -------------------------------

min_cpgs <- as.integer(Sys.getenv("MIN_CPGS", "10"))
min_read_depth <- as.numeric(Sys.getenv("MIN_READ_DEPTH", "5"))

N_jobs <- as.integer(Sys.getenv("N_GAM_JOBS", "900"))

PATH_wk <- path.expand(
  Sys.getenv("PATH_WK", "~/scratch/UQAC/meth/")
)

PATH_scr11 <- file.path(PATH_wk, "scr/11_mgcv")

region_file_path <- path.expand(
  Sys.getenv(
    "REGION_FILE",
    file.path(PATH_scr11, "dat/region_file_1_chunk.csv")
  )
)

pheno_file_path <- path.expand(
  Sys.getenv(
    "PHENO_FILE",
    file.path(PATH_scr11, "dat/18_pheno_BMI.RData")
  )
)

PATH_output <- path.expand(
  Sys.getenv(
    "WEIGHTED_GAM_OUTPUT",
    file.path(
      PATH_wk,
      "results/15_revision/7_read_depth_weighted_GAM"
    )
  )
)

overwrite <- Sys.getenv(
  "OVERWRITE", "0"
) %in% c("1", "TRUE", "true", "T", "yes", "YES")

if (!is.finite(min_cpgs) || min_cpgs < 5L) {
  stop("MIN_CPGS must be >= 5.")
}

if (!is.finite(min_read_depth) || min_read_depth < 1) {
  stop("MIN_READ_DEPTH must be >= 1.")
}

if (!is.finite(N_jobs) || N_jobs < 1L) {
  stop("N_GAM_JOBS must be >= 1.")
}

# ------------------------------ array task ID --------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1L) {
  stop(
    "No command-line arguments supplied. ",
    "Pass the SLURM_ARRAY_TASK_ID."
  )
}

i_job_id <- as.integer(args[1L])

if (
  !is.finite(i_job_id) ||
    i_job_id < 1L ||
    i_job_id > N_jobs
) {
  stop(
    "Job ID must be between 1 and ",
    N_jobs,
    "."
  )
}

cat(
  "Weighted GAM-DMR job ID: ",
  i_job_id,
  " / ",
  N_jobs,
  "\n",
  sep = ""
)

# -------------------------------- directories --------------------------------

dir.create(
  PATH_output,
  recursive = TRUE,
  showWarnings = FALSE
)

job_result_file <- file.path(
  PATH_output,
  sprintf(
    "weighted_GAM_results_job_%04d.tsv",
    i_job_id
  )
)

job_qc_file <- file.path(
  PATH_output,
  sprintf(
    "weighted_GAM_qc_job_%04d.tsv",
    i_job_id
  )
)

if (
  !overwrite &&
    file.exists(job_result_file) &&
    file.exists(job_qc_file)
) {
  cat(
    "Outputs already exist for job ",
    i_job_id,
    "; exiting. Set OVERWRITE=1 to rerun.\n",
    sep = ""
  )
  quit(
    save = "no",
    status = 0L
  )
}

if (overwrite) {
  unlink(
    c(
      job_result_file,
      job_qc_file
    ),
    force = TRUE
  )
}

# ------------------------------- load metadata -------------------------------

if (!file.exists(region_file_path)) {
  stop(
    "Region file does not exist: ",
    region_file_path
  )
}

region_file <- fread(region_file_path)

required_region_cols <- c(
  "chr",
  "region_start",
  "region_end",
  "data_chunk_id"
)

missing_region_cols <- setdiff(
  required_region_cols,
  names(region_file)
)

if (length(missing_region_cols)) {
  stop(
    "Region file is missing: ",
    paste(
      missing_region_cols,
      collapse = ", "
    )
  )
}

if (!file.exists(pheno_file_path)) {
  stop(
    "Phenotype RData does not exist: ",
    pheno_file_path
  )
}

pheno_env <- new.env(
  parent = emptyenv()
)

load(
  pheno_file_path,
  envir = pheno_env
)

if (!exists(
  "pheno_file",
  envir = pheno_env,
  inherits = FALSE
)) {
  stop(
    "Object pheno_file was not found in ",
    pheno_file_path
  )
}

pheno_dt <- as.data.table(
  copy(
    get(
      "pheno_file",
      envir = pheno_env
    )
  )
)

rm(pheno_env)

required_pheno_cols <- c(
  "ID",
  "AA_only",
  "AgeCalc",
  "Sex",
  "Non.smoker",
  "EOSINOpc",
  "LYMPHOpc",
  "MONOpc",
  "NEUTROpc",
  "sv1",
  "sv2",
  "sv3",
  "sv4",
  "sv5",
  "BMI",
  "FID"
)

missing_pheno_cols <- setdiff(
  required_pheno_cols,
  names(pheno_dt)
)

if (length(missing_pheno_cols)) {
  stop(
    "pheno_file is missing required columns: ",
    paste(
      missing_pheno_cols,
      collapse = ", "
    )
  )
}

pheno_dt[
  ,
  ID := as.character(ID)
]

pheno_dt <- unique(
  pheno_dt,
  by = "ID"
)

# ----------------------------- region assignment -----------------------------

region_indices <- seq_len(
  nrow(region_file)
)

# Round-robin assignment is robust even when N_jobs is changed.
job_assignment <- (
  (
    region_indices - 1L
  ) %% N_jobs
) + 1L

i_region_id_vector <- region_indices[
  job_assignment == i_job_id
]

cat(
  "Regions assigned to this job: ",
  length(i_region_id_vector),
  "\n",
  sep = ""
)

# ------------------------------- model formula --------------------------------

gam_formula <- meth_arcsin ~
  s(start, bs = "cs", k = k_cpg) +
  s(start, by = AA_only, bs = "cs", k = k_cpg) +
  AA_only + AgeCalc + Sex + Non.smoker +
  EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc +
  sv1 + sv2 + sv3 + sv4 + sv5 + BMI +
  s(FID, bs = "re")

feature_list <- c(
  "(Intercept)",
  "AA_only",
  "AgeCalc",
  "Sex",
  "Non.smoker",
  "EOSINOpc",
  "LYMPHOpc",
  "MONOpc",
  "NEUTROpc",
  "sv1",
  "sv2",
  "sv3",
  "sv4",
  "sv5",
  "BMI"
)

feature_output_names <- c(
  "Intercept",
  "AA_only",
  "AgeCalc",
  "Sex",
  "Non_smoker",
  "EOSINOpc",
  "LYMPHOpc",
  "MONOpc",
  "NEUTROpc",
  "sv1",
  "sv2",
  "sv3",
  "sv4",
  "sv5",
  "BMI"
)

feature_smooth_terms <- c(
  "s(start)",
  "s(start):AA_only",
  "s(FID)"
)

# ---------------------------------- helpers ----------------------------------

append_tsv <- function(dt, file) {
  fwrite(
    dt,
    file = file,
    append = file.exists(file),
    sep = "\t",
    col.names = !file.exists(file),
    quote = FALSE,
    na = "NA"
  )
}

safe_num <- function(x) {
  suppressWarnings(
    as.numeric(x)
  )
}

# -------------------------------- main loop ----------------------------------

for (i_region in i_region_id_vector) {

  region_start_time <- Sys.time()

  i_region_info <- region_file[i_region]
  i_chunk_id <- as.integer(
    i_region_info$data_chunk_id
  )

  cat(
    "Region ",
    i_region,
    " | chunk ",
    i_chunk_id,
    "\n",
    sep = ""
  )

  chunk_file <- file.path(
    PATH_wk,
    "data/meth_split",
    sprintf(
      "chunk_%04d.csv",
      i_chunk_id
    )
  )

  if (!file.exists(chunk_file)) {
    append_tsv(
      data.table(
        region_id = i_region,
        data_chunk_id = i_chunk_id,
        status = "MISSING_CHUNK",
        message = chunk_file
      ),
      job_qc_file
    )
    next
  }

  i_chunk <- tryCatch(
    fread(
      chunk_file,
      showProgress = FALSE
    ),
    error = function(e) e
  )

  if (inherits(i_chunk, "error")) {
    append_tsv(
      data.table(
        region_id = i_region,
        data_chunk_id = i_chunk_id,
        status = "CHUNK_READ_ERROR",
        message = conditionMessage(i_chunk)
      ),
      job_qc_file
    )
    next
  }

  i_region_data <- i_chunk[
    chr == i_region_info$chr &
      start >= i_region_info$region_start &
      start <= i_region_info$region_end
  ]

  if (nrow(i_region_data) < min_cpgs) {
    append_tsv(
      data.table(
        region_id = i_region,
        data_chunk_id = i_chunk_id,
        status = "TOO_FEW_CPGS_BEFORE_LONG",
        n_cpgs = nrow(i_region_data),
        message = NA_character_
      ),
      job_qc_file
    )
    rm(
      i_chunk,
      i_region_data
    )
    next
  }

  # -------------------------- paired meth/depth columns ----------------------

  meth_cols <- grep(
    "_meth$",
    names(i_region_data),
    value = TRUE
  )

  if (!length(meth_cols)) {
    append_tsv(
      data.table(
        region_id = i_region,
        data_chunk_id = i_chunk_id,
        status = "NO_METH_COLUMNS",
        message = NA_character_
      ),
      job_qc_file
    )
    rm(
      i_chunk,
      i_region_data
    )
    next
  }

  sample_ids <- sub(
    "_meth$",
    "",
    meth_cols
  )

  tot_cols_expected <- paste0(
    sample_ids,
    "_tot"
  )

  has_tot <- tot_cols_expected %in%
    names(i_region_data)

  if (!all(has_tot)) {
    cat(
      "  Missing paired *_tot columns for ",
      sum(!has_tot),
      " samples; these samples are omitted.\n",
      sep = ""
    )
  }

  meth_cols <- meth_cols[
    has_tot
  ]

  sample_ids <- sample_ids[
    has_tot
  ]

  tot_cols <- tot_cols_expected[
    has_tot
  ]

  if (!length(meth_cols)) {
    append_tsv(
      data.table(
        region_id = i_region,
        data_chunk_id = i_chunk_id,
        status = "NO_PAIRED_DEPTH_COLUMNS",
        message = NA_character_
      ),
      job_qc_file
    )
    rm(
      i_chunk,
      i_region_data
    )
    next
  }

  i_region_data[
    ,
    (meth_cols) := lapply(
      .SD,
      safe_num
    ),
    .SDcols = meth_cols
  ]

  i_region_data[
    ,
    (tot_cols) := lapply(
      .SD,
      safe_num
    ),
    .SDcols = tot_cols
  ]

  # Long methylation table.
  meth_long <- melt(
    i_region_data,
    id.vars = "start",
    measure.vars = meth_cols,
    variable.name = "meth_col",
    value.name = "meth",
    variable.factor = FALSE,
    na.rm = FALSE
  )

  meth_long[
    ,
    sample := sub(
      "_meth$",
      "",
      meth_col
    )
  ]

  meth_long[
    ,
    meth_col := NULL
  ]

  # Long total-read-depth table.
  depth_long <- melt(
    i_region_data,
    id.vars = "start",
    measure.vars = tot_cols,
    variable.name = "depth_col",
    value.name = "read_depth",
    variable.factor = FALSE,
    na.rm = FALSE
  )

  depth_long[
    ,
    sample := sub(
      "_tot$",
      "",
      depth_col
    )
  ]

  depth_long[
    ,
    depth_col := NULL
  ]

  setkey(
    meth_long,
    start,
    sample
  )

  setkey(
    depth_long,
    start,
    sample
  )

  long_dt <- merge(
    meth_long,
    depth_long,
    by = c(
      "start",
      "sample"
    ),
    all = FALSE,
    sort = FALSE
  )

  # Join the same phenotype/covariate table used in the original GAM-DMR.
  i_meth_long_cov <- merge(
    long_dt,
    pheno_dt,
    by.x = "sample",
    by.y = "ID",
    all = FALSE,
    sort = FALSE
  )

  required_complete_cols <- c(
    "start",
    "sample",
    "meth",
    "read_depth",
    "AA_only",
    "AgeCalc",
    "Sex",
    "Non.smoker",
    "EOSINOpc",
    "LYMPHOpc",
    "MONOpc",
    "NEUTROpc",
    "sv1",
    "sv2",
    "sv3",
    "sv4",
    "sv5",
    "BMI",
    "FID"
  )

  # Explicitly mirror the stated depth>=5 QC while also ensuring valid beta
  # values for the ASR transformation.
  n_before_depth_filter <- nrow(
    i_meth_long_cov
  )

  i_meth_long_cov <- i_meth_long_cov[
    complete.cases(
      i_meth_long_cov[
        ,
        ..required_complete_cols
      ]
    ) &
      is.finite(meth) &
      meth >= 0 &
      meth <= 1 &
      is.finite(read_depth) &
      read_depth >= min_read_depth
  ]

  n_after_depth_filter <- nrow(
    i_meth_long_cov
  )

  n_cpg_unique <- uniqueN(
    i_meth_long_cov$start
  )

  n_samples_unique <- uniqueN(
    i_meth_long_cov$sample
  )

  if (
    n_cpg_unique < min_cpgs ||
      n_after_depth_filter < min_cpgs
  ) {
    append_tsv(
      data.table(
        region_id = i_region,
        data_chunk_id = i_chunk_id,
        status = "TOO_FEW_OBSERVATIONS_AFTER_QC",
        n_cpgs = n_cpg_unique,
        n_samples = n_samples_unique,
        n_obs_before_depth_filter = n_before_depth_filter,
        n_obs_after_depth_filter = n_after_depth_filter,
        message = NA_character_
      ),
      job_qc_file
    )

    rm(
      i_chunk,
      i_region_data,
      meth_long,
      depth_long,
      long_dt,
      i_meth_long_cov
    )

    next
  }

  # ----------------------------- ASR + weights -------------------------------

  i_meth_long_cov[
    ,
    meth_arcsin := asin(
      sqrt(meth)
    )
  ]

  mean_region_depth <- mean(
    i_meth_long_cov$read_depth
  )

  # Relative inverse-variance weights. Mean weight = 1 in each region.
  i_meth_long_cov[
    ,
    depth_weight :=
      read_depth /
      mean_region_depth
  ]

  # Fix rare matrix/array columns inherited from upstream objects.
  for (v in names(i_meth_long_cov)) {
    if (
      is.matrix(
        i_meth_long_cov[[v]]
      ) ||
        is.array(
          i_meth_long_cov[[v]]
        )
    ) {
      i_meth_long_cov[[v]] <- as.numeric(
        i_meth_long_cov[[v]]
      )
    }
  }

  # Same basis-dimension rule as the original script.
  k_cpg <- max(
    5,
    min(
      10,
      floor(
        n_cpg_unique /
          20
      )
    )
  )

  # ------------------------------ weighted GAM -------------------------------

  warning_messages <- character()

  gam_model <- tryCatch(
    withCallingHandlers(
      mgcv::gam(
        gam_formula,
        data = i_meth_long_cov,
        weights = depth_weight,
        method = "REML"
      ),
      warning = function(w) {
        warning_messages <<- c(
          warning_messages,
          conditionMessage(w)
        )
        invokeRestart(
          "muffleWarning"
        )
      }
    ),
    error = function(e) e
  )

  if (inherits(gam_model, "error")) {
    append_tsv(
      data.table(
        region_id = i_region,
        data_chunk_id = i_chunk_id,
        status = "GAM_ERROR",
        n_cpgs = n_cpg_unique,
        n_samples = n_samples_unique,
        n_obs_after_depth_filter = n_after_depth_filter,
        message = conditionMessage(gam_model)
      ),
      job_qc_file
    )

    rm(
      i_chunk,
      i_region_data,
      meth_long,
      depth_long,
      long_dt,
      i_meth_long_cov,
      gam_model
    )

    next
  }

  sm <- summary(
    gam_model
  )

  param_coef <- sm$p.table
  smooth_terms <- sm$s.table

  # --------------------------- fixed-effect p-values -------------------------

  fixed_p <- rep(
    1,
    length(feature_list)
  )

  names(fixed_p) <- feature_output_names

  present_features <- intersect(
    feature_list,
    rownames(param_coef)
  )

  if (length(present_features)) {
    for (jj in seq_along(feature_list)) {
      old_name <- feature_list[jj]

      if (old_name %in% present_features) {
        fixed_p[
          jj
        ] <- param_coef[
          old_name,
          "Pr(>|t|)"
        ]
      }
    }
  }

  # ---------------------------- smooth statistics ----------------------------

  smooth_p <- rep(
    1,
    length(feature_smooth_terms)
  )

  names(smooth_p) <- c(
    "p_s_start",
    "p_s_start_AA",
    "p_s_FID"
  )

  smooth_edf <- rep(
    0,
    length(feature_smooth_terms)
  )

  names(smooth_edf) <- c(
    "edf_s_start",
    "edf_s_start_AA",
    "edf_s_FID"
  )

  present_smooth <- intersect(
    feature_smooth_terms,
    rownames(smooth_terms)
  )

  for (jj in seq_along(feature_smooth_terms)) {
    nm <- feature_smooth_terms[jj]

    if (nm %in% present_smooth) {
      smooth_p[
        jj
      ] <- smooth_terms[
        nm,
        "p-value"
      ]

      smooth_edf[
        jj
      ] <- smooth_terms[
        nm,
        "edf"
      ]
    }
  }

  # -------------------------- effect-size estimates --------------------------

  newdata_AA1 <- copy(
    i_meth_long_cov
  )

  newdata_AA0 <- copy(
    i_meth_long_cov
  )

  newdata_AA1[
    ,
    AA_only := 1
  ]

  newdata_AA0[
    ,
    AA_only := 0
  ]

  pred_AA1 <- predict(
    gam_model,
    newdata = newdata_AA1,
    type = "response"
  )

  pred_AA0 <- predict(
    gam_model,
    newdata = newdata_AA0,
    type = "response"
  )

  delta <- pred_AA1 -
    pred_AA0

  max_diff <- max(
    abs(delta),
    na.rm = TRUE
  )

  mean_diff <- mean(
    abs(delta),
    na.rm = TRUE
  )

  # ------------------------------- save result -------------------------------

  result_dt <- as.data.table(
    as.list(
      c(
        setNames(
          fixed_p,
          paste0(
            "p_",
            names(fixed_p)
          )
        ),
        smooth_p,
        smooth_edf
      )
    )
  )

  result_dt[
    ,
    `:=`(
      R2 = sm$r.sq,
      AIC = AIC(gam_model),
      Deviance_explained = sm$dev.expl,
      REML = gam_model$gcv.ubre,
      N_cpgs = n_cpg_unique,
      N_samples = n_samples_unique,
      N_observations = n_after_depth_filter,
      data_chunk_id = i_chunk_id,
      region_id = i_region,
      max_diff = max_diff,
      mean_diff = mean_diff,
      mean_read_depth = mean_region_depth,
      median_read_depth = median(
        i_meth_long_cov$read_depth
      ),
      q25_read_depth = as.numeric(
        quantile(
          i_meth_long_cov$read_depth,
          0.25,
          names = FALSE
        )
      ),
      q75_read_depth = as.numeric(
        quantile(
          i_meth_long_cov$read_depth,
          0.75,
          names = FALSE
        )
      ),
      min_read_depth = min(
        i_meth_long_cov$read_depth
      ),
      max_read_depth = max(
        i_meth_long_cov$read_depth
      ),
      weight_mean = mean(
        i_meth_long_cov$depth_weight
      ),
      weight_sd = sd(
        i_meth_long_cov$depth_weight
      ),
      weight_min = min(
        i_meth_long_cov$depth_weight
      ),
      weight_max = max(
        i_meth_long_cov$depth_weight
      ),
      k_cpg = k_cpg
    )
  ]

  append_tsv(
    result_dt,
    job_result_file
  )

  append_tsv(
    data.table(
      region_id = i_region,
      data_chunk_id = i_chunk_id,
      status = "COMPLETED",
      n_cpgs = n_cpg_unique,
      n_samples = n_samples_unique,
      n_obs_before_depth_filter = n_before_depth_filter,
      n_obs_after_depth_filter = n_after_depth_filter,
      n_warnings = length(
        warning_messages
      ),
      warning_message = if (
        length(
          warning_messages
        )
      ) {
        paste(
          unique(
            warning_messages
          ),
          collapse = " | "
        )
      } else {
        NA_character_
      },
      elapsed_seconds = as.numeric(
        difftime(
          Sys.time(),
          region_start_time,
          units = "secs"
        )
      ),
      message = NA_character_
    ),
    job_qc_file
  )

  rm(
    i_chunk,
    i_region_data,
    meth_long,
    depth_long,
    long_dt,
    i_meth_long_cov,
    newdata_AA1,
    newdata_AA0,
    gam_model
  )

  invisible(
    gc(FALSE)
  )
}

cat(
  "\nWeighted GAM-DMR job completed: ",
  i_job_id,
  "\n",
  sep = ""
)

###############################################################################
# 2_run_mgcv.R
#
# Leakage-free GAM-DMR discovery for ONE outer split.
#
# Usage:
#
# Rscript 2_run_mgcv.R SPLIT_ID JOB_ID N_JOBS
#
# Example:
#
# Rscript 2_run_mgcv.R 1 1 990
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
})

###############################################################################
# Parameters
###############################################################################

min_cpgs <- 10L

PATH_wk <- path.expand("~/scratch/UQAC/meth/")

PATH_scr11 <- file.path(
  PATH_wk,
  "scr/11_mgcv"
)

PATH_validation <- file.path(
  PATH_wk,
  "results/15_revision/5_cv_leakfree"
)

PATH_output <- file.path(
  PATH_validation,
  "1_mgcv_raw"
)

dir.create(
  PATH_output,
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################################
# Command line arguments
###############################################################################

args <- commandArgs(
  trailingOnly = TRUE
)

if (length(args) < 3L) {

  stop(
    "Usage: Rscript 2_run_mgcv.R ",
    "SPLIT_ID JOB_ID N_JOBS"
  )
}

split_id <- as.integer(args[1])
job_id   <- as.integer(args[2])
N_jobs   <- as.integer(args[3])

if (
  anyNA(
    c(
      split_id,
      job_id,
      N_jobs
    )
  )
) {
  stop("Invalid command-line arguments.")
}

cat(
  "Split:", split_id,
  " Job:", job_id,
  " N_jobs:", N_jobs,
  "\n"
)

###############################################################################
# IMPORTANT:
#
# Use ALL predefined regions.
#
# DO NOT use:
#
# region_GAM_M12.csv
# full-cohort GAM DMRs
# M1_manhattan
# M2_manhattan
###############################################################################

region_file <- fread(
  file.path(
    PATH_scr11,
    "dat/region_file_1_chunk.csv"
  )
)

required_region_cols <- c(
  "data_chunk_id",
  "chr",
  "region_start",
  "region_end",
  "region_id"
)

missing_region_cols <- setdiff(
  required_region_cols,
  names(region_file)
)

if (length(missing_region_cols) > 0L) {

  stop(
    "Missing region columns: ",
    paste(
      missing_region_cols,
      collapse = ", "
    )
  )
}

###############################################################################
# Give every predefined region a permanent row identifier
###############################################################################

region_file[
  ,
  region_row_id := .I
]

###############################################################################
# Split regions over SLURM jobs
###############################################################################

N_regions <- nrow(region_file)

start_index <- floor(
  (job_id - 1) * N_regions / N_jobs
) + 1L

end_index <- floor(
  job_id * N_regions / N_jobs
)

if (start_index > end_index) {

  cat("No regions assigned to this job.\n")
  quit(save = "no")
}

region_indices <- seq.int(
  start_index,
  end_index
)

cat(
  "Regions:",
  start_index,
  "to",
  end_index,
  "\n"
)

###############################################################################
# Load outer splits
###############################################################################

split_file <- fread(
  file.path(
    PATH_validation,
    "outer_family_membership.csv"
  )
)

train_FIDs <- split_file[
  splitID == split_id &
  set == "train",
  unique(as.character(FID))
]

test_FIDs <- split_file[
  splitID == split_id &
  set == "test",
  unique(as.character(FID))
]

stopifnot(
  length(
    intersect(
      train_FIDs,
      test_FIDs
    )
  ) == 0L
)

###############################################################################
# Load phenotype
###############################################################################

pheno_env <- new.env()

load(
  file.path(
    PATH_scr11,
    "dat/18_pheno_BMI.RData"
  ),
  envir = pheno_env
)

if (
  !exists(
    "pheno_file",
    envir = pheno_env
  )
) {
  stop(
    "pheno_file missing from 18_pheno_BMI.RData"
  )
}

pheno <- as.data.table(
  pheno_env$pheno_file
)

rm(pheno_env)

###############################################################################
# Only variables needed by GAM
###############################################################################

pheno_vars <- c(
  "ID",
  "FID",
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

missing_pheno_vars <- setdiff(
  pheno_vars,
  names(pheno)
)

if (length(missing_pheno_vars) > 0L) {

  stop(
    "Missing phenotype variables: ",
    paste(
      missing_pheno_vars,
      collapse = ", "
    )
  )
}

pheno <- pheno[
  ,
  ..pheno_vars
]

pheno[
  ,
  ID := as.character(ID)
]

pheno[
  ,
  FID := as.character(FID)
]

###############################################################################
# CRITICAL:
#
# Remove all test-family subjects before touching the methylation matrix.
###############################################################################

pheno_train <- pheno[
  FID %in% train_FIDs
]

if (
  any(
    pheno_train$FID %in% test_FIDs
  )
) {
  stop("TEST FAMILY LEAKAGE DETECTED")
}

setnames(
  pheno_train,
  "ID",
  "sample"
)

###############################################################################
# Function for empty/error rows
###############################################################################

make_result_row <- function(
  region_info,
  status,
  error_message = NA_character_
) {

  data.table(
    splitID = split_id,
    jobID = job_id,
    region_row_id = region_info$region_row_id,
    region_id = region_info$region_id,
    data_chunk_id = region_info$data_chunk_id,
    chr = region_info$chr,
    region_start = region_info$region_start,
    region_end = region_info$region_end,

    p_AA = NA_real_,
    p_smooth_AA = NA_real_,
    edf_smooth_AA = NA_real_,

    R2 = NA_real_,
    AIC = NA_real_,
    deviance_explained = NA_real_,

    max_diff = NA_real_,
    mean_diff = NA_real_,

    N_cpgs = NA_integer_,
    N_samples = NA_integer_,
    N_train_FIDs = NA_integer_,

    status = status,
    error_message = error_message
  )
}

###############################################################################
# Process regions
###############################################################################

result_list <- vector(
  "list",
  length(region_indices)
)

current_chunk_id <- NA_integer_
current_chunk <- NULL

counter <- 0L

for (i_region in region_indices) {

  counter <- counter + 1L

  region_info <- region_file[
    i_region
  ]

  chunk_id <- region_info$data_chunk_id

  ###########################################################################
  # Load methylation chunk only when chunk changes
  ###########################################################################

  if (
    is.na(current_chunk_id) ||
    chunk_id != current_chunk_id
  ) {

    chunk_path <- file.path(
      PATH_wk,
      "data/meth_split",
      paste0(
        "chunk_",
        sprintf(
          "%04d",
          chunk_id
        ),
        ".csv"
      )
    )

    if (!file.exists(chunk_path)) {

      result_list[[counter]] <-
        make_result_row(
          region_info,
          "chunk_missing",
          chunk_path
        )

      next
    }

    current_chunk <- fread(
      chunk_path
    )

    current_chunk_id <- chunk_id
  }

  ###########################################################################
  # Extract predefined region
  ###########################################################################

  region_dat <- current_chunk[
    chr == region_info$chr &
    start >= region_info$region_start &
    start <= region_info$region_end
  ]

  n_raw_cpgs <- uniqueN(
    region_dat$start
  )

  if (
    n_raw_cpgs < min_cpgs
  ) {

    rr <- make_result_row(
      region_info,
      "fewer_than_10_CpGs"
    )

    rr[, N_cpgs := n_raw_cpgs]

    result_list[[counter]] <- rr

    next
  }

  ###########################################################################
  # CRITICAL:
  #
  # Use methylation columns belonging to TRAINING subjects only.
  #
  # Test-family methylation never enters the GAM fitting dataset.
  ###########################################################################

  train_meth_cols <- intersect(
    paste0(
      pheno_train$sample,
      "_meth"
    ),
    names(region_dat)
  )

  if (
    length(train_meth_cols) == 0L
  ) {

    result_list[[counter]] <-
      make_result_row(
        region_info,
        "no_training_samples"
      )

    next
  }

  ###########################################################################
  # Convert training methylation columns
  ###########################################################################

  region_dat[
    ,
    (train_meth_cols) :=
      lapply(
        .SD,
        as.numeric
      ),
    .SDcols = train_meth_cols
  ]

  ###########################################################################
  # Long format
  ###########################################################################

  meth_long <- melt(
    region_dat,
    id.vars = "start",
    measure.vars = train_meth_cols,
    variable.name = "sample",
    value.name = "meth",
    na.rm = TRUE
  )

  meth_long[
    ,
    sample :=
      sub(
        "_meth$",
        "",
        as.character(sample)
      )
  ]

  ###########################################################################
  # Merge TRAINING phenotype only
  ###########################################################################

  dat <- merge(
    meth_long,
    pheno_train,
    by = "sample",
    all = FALSE
  )

  ###########################################################################
  # Convert matrix/array phenotype variables to numeric vectors
  ###########################################################################

  for (v in names(dat)) {

    if (
      is.matrix(dat[[v]]) ||
      is.array(dat[[v]])
    ) {

      dat[
        ,
        (v) := as.numeric(dat[[v]])
      ]
    }
  }

  ###########################################################################
  # Complete cases only for model variables
  ###########################################################################

  model_vars <- c(
    "meth",
    "start",
    "FID",
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

  dat <- dat[
    complete.cases(
      dat[
        ,
        ..model_vars
      ]
    )
  ]

  ###########################################################################
  # Training-only eligibility
  ###########################################################################

  n_cpgs <- uniqueN(
    dat$start
  )

  n_samples <- uniqueN(
    dat$sample
  )

  n_FIDs <- uniqueN(
    dat$FID
  )

  if (
    n_cpgs < min_cpgs
  ) {

    rr <- make_result_row(
      region_info,
      "fewer_than_10_training_CpGs"
    )

    rr[, N_cpgs := n_cpgs]
    rr[, N_samples := n_samples]
    rr[, N_train_FIDs := n_FIDs]

    result_list[[counter]] <- rr

    next
  }

  if (
    uniqueN(dat$AA_only) < 2L
  ) {

    result_list[[counter]] <-
      make_result_row(
        region_info,
        "only_one_AA_class"
      )

    next
  }

  if (
    n_FIDs < 2L
  ) {

    result_list[[counter]] <-
      make_result_row(
        region_info,
        "fewer_than_2_families"
      )

    next
  }

  ###########################################################################
  # Prepare variables
  ###########################################################################

  dat[
    ,
    FID := droplevels(
      factor(FID)
    )
  ]

  dat[
    ,
    AA_only := as.numeric(
      AA_only
    )
  ]

  dat[
    ,
    meth_arcsin :=
      asin(
        sqrt(meth)
      )
  ]

  ###########################################################################
  # Same CpG-adaptive basis rule as primary GAM-DMR
  ###########################################################################

  k_cpg <- max(
    5L,
    min(
      10L,
      floor(
        n_cpgs / 20L
      )
    )
  )

  ###########################################################################
  # Fit GAM using TRAINING subjects only
  ###########################################################################

  gam_model <- tryCatch(

    mgcv::gam(

      meth_arcsin ~

        s(
          start,
          bs = "cs",
          k = k_cpg
        ) +

        s(
          start,
          by = AA_only,
          bs = "cs",
          k = k_cpg
        ) +

        AA_only +

        AgeCalc +
        Sex +
        Non.smoker +

        EOSINOpc +
        LYMPHOpc +
        MONOpc +
        NEUTROpc +

        sv1 +
        sv2 +
        sv3 +
        sv4 +
        sv5 +

        BMI +

        s(
          FID,
          bs = "re"
        ),

      data = dat,

      method = "REML"
    ),

    error = function(e) e
  )

  if (
    inherits(
      gam_model,
      "error"
    )
  ) {

    result_list[[counter]] <-
      make_result_row(
        region_info,
        "gam_error",
        conditionMessage(gam_model)
      )

    next
  }

  ###########################################################################
  # Extract GAM statistics
  ###########################################################################

  sm <- summary(
    gam_model
  )

  p_AA <- NA_real_

  if (
    "AA_only" %in%
    rownames(sm$p.table)
  ) {

    p_AA <-
      sm$p.table[
        "AA_only",
        "Pr(>|t|)"
      ]
  }

  smooth_name <- "s(start):AA_only"

  p_smooth_AA <- NA_real_
  edf_smooth_AA <- NA_real_

  if (
    smooth_name %in%
    rownames(sm$s.table)
  ) {

    p_smooth_AA <-
      sm$s.table[
        smooth_name,
        "p-value"
      ]

    edf_smooth_AA <-
      sm$s.table[
        smooth_name,
        "edf"
      ]
  }

  ###########################################################################
  # Training-only AA effect estimate
  ###########################################################################

  new_AA1 <- copy(dat)
  new_AA0 <- copy(dat)

  new_AA1[
    ,
    AA_only := 1
  ]

  new_AA0[
    ,
    AA_only := 0
  ]

  pred_AA1 <- predict(
    gam_model,
    newdata = new_AA1,
    type = "response"
  )

  pred_AA0 <- predict(
    gam_model,
    newdata = new_AA0,
    type = "response"
  )

  delta <- pred_AA1 - pred_AA0

  ###########################################################################
  # Save statistics
  ###########################################################################

  result_list[[counter]] <-
    data.table(

      splitID = split_id,
      jobID = job_id,

      region_row_id =
        region_info$region_row_id,

      region_id =
        region_info$region_id,

      data_chunk_id =
        region_info$data_chunk_id,

      chr =
        region_info$chr,

      region_start =
        region_info$region_start,

      region_end =
        region_info$region_end,

      p_AA =
        p_AA,

      p_smooth_AA =
        p_smooth_AA,

      edf_smooth_AA =
        edf_smooth_AA,

      R2 =
        sm$r.sq,

      AIC =
        AIC(gam_model),

      deviance_explained =
        sm$dev.expl,

      max_diff =
        max(
          abs(delta),
          na.rm = TRUE
        ),

      mean_diff =
        mean(
          abs(delta),
          na.rm = TRUE
        ),

      N_cpgs =
        n_cpgs,

      N_samples =
        n_samples,

      N_train_FIDs =
        n_FIDs,

      status =
        "ok",

      error_message =
        NA_character_
    )

  if (
    counter %% 10L == 0L
  ) {

    cat(
      "Completed",
      counter,
      "of",
      length(region_indices),
      "\n"
    )
  }
}

###############################################################################
# Save one file per split/job
###############################################################################

result_table <- rbindlist(
  result_list,
  fill = TRUE
)

output_file <- file.path(
  PATH_output,
  sprintf(
    "gam_split_%03d_job_%04d.tsv",
    split_id,
    job_id
  )
)

fwrite(
  result_table,
  output_file,
  sep = "\t"
)

cat(
  "\nSaved:",
  output_file,
  "\n"
)

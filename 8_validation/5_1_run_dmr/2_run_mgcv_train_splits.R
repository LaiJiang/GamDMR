# This script runs the MGCV model for each selected region and each
# family-level training split. Each model is fitted using training families only.

min_cpgs <- 10L

# Number of SLURM array jobs. This must match the SLURM_ARRAY_TASK_ID range.
N_jobs <- 990L

# Define paths
# PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
PATH_wk <- "~/scratch/UQAC/meth/"

PATH_scr6 <- paste0(PATH_wk, "scr/6_whole/")
PATH_scr11 <- paste0(PATH_wk, "scr/11_mgcv/")
PATH_output <- paste0(PATH_wk, "results/15_revision/5_cv/1_mgcv/")

# Input file containing:
# splitID, train_FID_1, ..., train_FID_135
PATH_train_splits <- paste0(
  PATH_scr11,
  "dat/100_family_training_splits.csv"
)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0L) {
  stop(
    "No command-line argument supplied. ",
    "Please pass the SLURM_ARRAY_TASK_ID."
  )
}

i_job_id <- as.integer(args[1])

if (is.na(i_job_id)) {
  stop("The supplied SLURM_ARRAY_TASK_ID is not an integer.")
}

cat("Job ID:", i_job_id, "\n")

# Ensure output directory exists
if (!dir.exists(PATH_output)) {
  dir.create(PATH_output, recursive = TRUE)
}

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(stringr)
  library(tidyr)
  library(mgcv)
})

# Load region information
region_file <- fread(
  paste0(PATH_scr11, "dat/region_GAM_M12.csv")
)

required_region_columns <- c(
  "data_chunk_id",
  "chr",
  "region_start",
  "region_end",
  "region_id"
)

missing_region_columns <- setdiff(
  required_region_columns,
  names(region_file)
)

if (length(missing_region_columns) > 0L) {
  stop(
    "region_file is missing: ",
    paste(missing_region_columns, collapse = ", ")
  )
}

# Load the 100 family-level training splits once
train_FID_df <- fread(PATH_train_splits)

if (!"splitID" %in% names(train_FID_df)) {
  stop("train_FID_df must contain a splitID column.")
}

train_FID_cols <- grep(
  "^train_FID_[0-9]+$",
  names(train_FID_df),
  value = TRUE
)

if (length(train_FID_cols) == 0L) {
  stop("No train_FID_* columns were found in train_FID_df.")
}

# Put the training-FID columns in numeric order
train_FID_col_number <- as.integer(
  sub("^train_FID_", "", train_FID_cols)
)

train_FID_cols <- train_FID_cols[
  order(train_FID_col_number)
]

train_FID_df[, splitID := as.integer(splitID)]

if (anyNA(train_FID_df$splitID)) {
  stop("At least one splitID is missing or non-numeric.")
}

if (anyDuplicated(train_FID_df$splitID)) {
  stop("splitID values must be unique.")
}

cat(
  "Number of splits:", nrow(train_FID_df), "\n",
  "Training-FID columns per split:", length(train_FID_cols), "\n"
)

# Divide regions among SLURM array jobs
vec_list <- split(
  seq_len(nrow(region_file)),
  cut(
    seq_len(nrow(region_file)),
    breaks = N_jobs,
    labels = FALSE
  )
)

if (i_job_id < 1L || i_job_id > length(vec_list)) {
  stop(
    "Job ID ", i_job_id,
    " is outside the valid range 1 to ",
    length(vec_list), "."
  )
}

i_region_id_vector <- vec_list[[i_job_id]]

# Keep output separate from results produced by the original script
output_file <- paste0(
  PATH_output,
  "2_run_mgcv_training_splits_results_",
  i_job_id,
  ".txt"
)

# The file is headerless to preserve the format of the original script.
# Column order:
# splitID, fixed-effect p-values, smooth-term p-values, smooth-term EDFs,
# R2, AIC, deviance explained, REML, N_cpgs, N_samples,
# N_train_FIDs, data_chunk_id, region_id, max_diff, mean_diff.

for (i_region in i_region_id_vector) {

  i_region_info <- region_file[i_region]

  i_chunk_id <- i_region_info$data_chunk_id
  i_region_id <- i_region_info$region_id

  cat(
    "\nRegion row:", i_region,
    " region_id:", i_region_id,
    " chunk:", i_chunk_id,
    "\n"
  )

  # Load the corresponding methylation-data chunk
  i_chunk <- fread(
    paste0(
      PATH_wk,
      "data/meth_split/chunk_",
      sprintf("%04d", i_chunk_id),
      ".csv"
    )
  )

  # Extract the current region
  i_region_data <- i_chunk[
    chr == i_region_info$chr &
      start >= i_region_info$region_start &
      start <= i_region_info$region_end
  ]

  rm(i_chunk)

  if (nrow(i_region_data) < min_cpgs) {
    cat(
      "Skipping region", i_region_id,
      "because it contains fewer than",
      min_cpgs, "CpG rows before reshaping.\n"
    )
    next
  }

  setDT(i_region_data)

  # Convert sample methylation columns to numeric
  meth_cols <- grep(
    "_meth$",
    names(i_region_data),
    value = TRUE
  )

  if (length(meth_cols) == 0L) {
    cat(
      "Skipping region", i_region_id,
      "because no *_meth columns were found.\n"
    )
    next
  }

  i_region_data[
    ,
    (meth_cols) := lapply(.SD, as.numeric),
    .SDcols = meth_cols
  ]

  # Convert the regional methylation matrix to long format
  meth_long <- melt(
    i_region_data,
    id.vars = "start",
    measure.vars = meth_cols,
    variable.name = "sample",
    value.name = "meth",
    na.rm = TRUE
  )

  meth_long[
    ,
    sample := str_remove(
      as.character(sample),
      "_meth$"
    )
  ]

  rm(i_region_data)

  # Load phenotype data separately for the current region, as requested.
  # Loading into a temporary environment prevents the large unused matrices
  # in the RData file from remaining in the global environment.
  pheno_env <- new.env()

  load(
    file = paste0(
      PATH_scr11,
      "dat/18_pheno_BMI.RData"
    ),
    envir = pheno_env,
    verbose = TRUE
  )

  if (!exists("pheno_file", envir = pheno_env, inherits = FALSE)) {
    stop("18_pheno_BMI.RData does not contain pheno_file.")
  }

  pheno_file <- as_tibble(
    get("pheno_file", envir = pheno_env)
  )

  rm(pheno_env)
  invisible(gc(verbose = FALSE))

  # Convert FID to character before matching it to train_FID_df
  pheno_file <- pheno_file %>%
    mutate(
      ID = as.character(ID),
      FID = as.character(FID)
    )

  # Run one training-only MGCV model for every splitID
  for (i_split_id in train_FID_df$splitID) {

    split_row <- train_FID_df[
      splitID == i_split_id,
      ..train_FID_cols
    ]

    train_FIDs <- unlist(
      split_row,
      use.names = FALSE
    )

    train_FIDs <- unique(
      as.character(
        train_FIDs[
          !is.na(train_FIDs) &
            nzchar(as.character(train_FIDs))
        ]
      )
    )

    if (length(train_FIDs) == 0L) {
      cat(
        "Skipping split", i_split_id,
        "for region", i_region_id,
        "because no training FIDs were found.\n"
      )
      next
    }

    # Restrict phenotype data to training families only
    pheno_train <- pheno_file %>%
      filter(FID %in% train_FIDs) %>%
      rename(sample = ID)

    if (nrow(pheno_train) == 0L) {
      cat(
        "Skipping split", i_split_id,
        "for region", i_region_id,
        "because no phenotype rows matched the training FIDs.\n"
      )
      next
    }

    # Inner join guarantees that test-family samples are excluded
    i_meth_long_cov <- meth_long %>%
      inner_join(
        pheno_train,
        by = "sample"
      ) %>%
      drop_na()

    if (nrow(i_meth_long_cov) == 0L) {
      cat(
        "Skipping split", i_split_id,
        "for region", i_region_id,
        "because no complete training observations remain.\n"
      )
      next
    }

    # Convert any one-column matrix or array variables to numeric vectors
    for (v in names(i_meth_long_cov)) {
      if (
        is.matrix(i_meth_long_cov[[v]]) ||
          is.array(i_meth_long_cov[[v]])
      ) {
        i_meth_long_cov[[v]] <- as.numeric(
          i_meth_long_cov[[v]]
        )
      }
    }

    # Restore FID as a factor for s(FID, bs = "re")
    i_meth_long_cov$FID <- droplevels(
      factor(i_meth_long_cov$FID)
    )

    i_meth_long_cov$AA_only <- as.numeric(
      i_meth_long_cov$AA_only
    )

    n_cpgs_train <- data.table::uniqueN(
      i_meth_long_cov$start
    )

    n_samples_train <- data.table::uniqueN(
      i_meth_long_cov$sample
    )

    n_train_FIDs_observed <- data.table::uniqueN(
      i_meth_long_cov$FID
    )

    if (n_cpgs_train < min_cpgs) {
      cat(
        "Skipping split", i_split_id,
        "for region", i_region_id,
        "because only", n_cpgs_train,
        "unique CpGs remain after the training-data join.\n"
      )
      next
    }

    if (length(unique(i_meth_long_cov$AA_only)) < 2L) {
      cat(
        "Skipping split", i_split_id,
        "for region", i_region_id,
        "because AA_only has fewer than two values in training data.\n"
      )
      next
    }

    if (n_train_FIDs_observed < 2L) {
      cat(
        "Skipping split", i_split_id,
        "for region", i_region_id,
        "because fewer than two training FIDs remain.\n"
      )
      next
    }

    # Arcsine-square-root transformation
    i_meth_long_cov$meth_arcsin <- asin(
      sqrt(i_meth_long_cov$meth)
    )

    # Define the spline basis dimension from the number of training CpGs
    k_cpg <- max(
      5L,
      min(
        10L,
        floor(n_cpgs_train / 20L)
      )
    )

    # Fit the model using training families only.
    # tryCatch allows other splits and regions to continue if one fit fails.
    gam_model <- tryCatch(
      mgcv::gam(
        meth_arcsin ~
          s(start, bs = "cs", k = k_cpg) +
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
          s(FID, bs = "re"),
        data = i_meth_long_cov,
        method = "REML"
      ),
      error = function(e) {
        cat(
          "MGCV failed for region", i_region_id,
          "split", i_split_id, ":",
          conditionMessage(e), "\n"
        )
        NULL
      }
    )

    if (is.null(gam_model)) {
      next
    }

    gam_summary <- summary(gam_model)

    # Extract parametric and smooth-term statistics
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

    feature_smooth_terms <- c(
      "s(start)",
      "s(start):AA_only",
      "s(FID)"
    )

    param_coef <- gam_summary$p.table
    smooth_terms <- gam_summary$s.table

    present_features <- intersect(
      feature_list,
      rownames(param_coef)
    )

    result_pvalue_fixed_effect <- rep(
      1,
      length(feature_list)
    )

    names(result_pvalue_fixed_effect) <- feature_list

    result_pvalue_fixed_effect[present_features] <-
      param_coef[present_features, "Pr(>|t|)"]

    present_smooth <- intersect(
      feature_smooth_terms,
      rownames(smooth_terms)
    )

    result_smooth_terms <- rep(
      1,
      length(feature_smooth_terms)
    )

    names(result_smooth_terms) <- feature_smooth_terms

    result_smooth_terms[present_smooth] <-
      smooth_terms[present_smooth, "p-value"]

    edf_smooth_terms <- rep(
      0,
      length(feature_smooth_terms)
    )

    names(edf_smooth_terms) <- feature_smooth_terms

    edf_smooth_terms[present_smooth] <-
      smooth_terms[present_smooth, "edf"]

    # Estimate the AA effect using the training observations
    newdata_AA1 <- i_meth_long_cov
    newdata_AA1$AA_only <- 1

    newdata_AA0 <- i_meth_long_cov
    newdata_AA0$AA_only <- 0

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

    delta <- pred_AA1 - pred_AA0

    max_diff <- max(
      abs(delta),
      na.rm = TRUE
    )

    mean_diff <- mean(
      abs(delta),
      na.rm = TRUE
    )

    # Add splitID to every output row
    result_stats <- c(
      "splitID" = i_split_id,
      result_pvalue_fixed_effect,
      result_smooth_terms,
      edf_smooth_terms,
      "R2" = gam_summary$r.sq,
      "AIC" = AIC(gam_model),
      "Deviance_explained" = gam_summary$dev.expl,
      "REML" = gam_model$gcv.ubre,
      "N_cpgs" = n_cpgs_train,
      "N_samples" = n_samples_train,
      "N_train_FIDs" = n_train_FIDs_observed,
      "data_chunk_id" = i_chunk_id,
      "region_id" = i_region_id,
      "max_diff" = max_diff,
      "mean_diff" = mean_diff
    )

    result_stats_df <- as.data.frame(
      t(as.numeric(result_stats))
    )

    # Preserve the original headerless, append-only output structure
    fwrite(
      result_stats_df,
      file = output_file,
      append = TRUE,
      sep = "\t",
      col.names = FALSE
    )

    cat(
      "Completed region", i_region_id,
      "split", i_split_id,
      "with", n_train_FIDs_observed,
      "training FIDs and", n_samples_train,
      "training samples.\n"
    )

    rm(
      gam_model,
      gam_summary,
      i_meth_long_cov,
      pheno_train,
      newdata_AA1,
      newdata_AA0,
      pred_AA1,
      pred_AA0,
      delta
    )
  }

  rm(
    meth_long,
    pheno_file
  )

  invisible(gc(verbose = FALSE))
}

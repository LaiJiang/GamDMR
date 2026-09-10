#!/usr/bin/env Rscript

# =============================================================================
# 20_run_read_depth_variability_chunks.R
#
# Purpose
# -------
# Quantify read-depth heterogeneity for the methylation dataset used by GAM-DMR.
#
# Input format (same as existing methylation chunks):
#   ~/scratch/UQAC/meth/data/meth_split/chunk_XXXX.csv
#
# Expected wide columns:
#   chr, start, ...
#   <sample>_meth : methylation proportion
#   <sample>_tot  : total read depth
#
# By default, the analysis is restricted to the final analytic cohort contained
# in pheno_file from:
#   ~/scratch/UQAC/meth/scr/11_mgcv/dat/18_pheno_BMI.RData
#
# Two read-depth sets are summarized:
#   1) all_depth:
#      all finite non-negative *_tot values;
#   2) gam_observed:
#      *_tot values only when the paired *_meth value is finite and in [0,1].
#      This second set most closely reflects observations that can contribute to
#      the GAM-DMR methylation-proportion analysis.
#
# One SLURM array task processes a subset of methylation chunks.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

if (!requireNamespace("matrixStats", quietly = TRUE)) {
  stop(
    "Package 'matrixStats' is required. It is normally available in the ",
    "StdEnv/2023 + r-bundle-bioconductor/3.20 environment."
  )
}

data.table::setDTthreads(1L)
options(mc.cores = 1L)

# ------------------------------- configuration -------------------------------

PATH_wk <- path.expand(
  Sys.getenv("PATH_WK", "~/scratch/UQAC/meth/")
)

chunk_dir <- path.expand(
  Sys.getenv(
    "METH_CHUNK_DIR",
    file.path(PATH_wk, "data/meth_split")
  )
)

pheno_file_path <- path.expand(
  Sys.getenv(
    "PHENO_FILE",
    file.path(PATH_wk, "scr/11_mgcv/dat/18_pheno_BMI.RData")
  )
)

output_root <- path.expand(
  Sys.getenv(
    "READ_DEPTH_OUTPUT",
    file.path(
      PATH_wk,
      "results/15_revision/6_read_depth_variability"
    )
  )
)

n_jobs <- as.integer(
  Sys.getenv("N_READ_DEPTH_JOBS", "300")
)

max_hist_depth <- as.integer(
  Sys.getenv("MAX_HIST_DEPTH", "200")
)

analytic_only <- Sys.getenv(
  "ANALYTIC_ONLY", "1"
) %in% c("1", "TRUE", "true", "T", "yes", "YES")

overwrite <- Sys.getenv(
  "OVERWRITE", "0"
) %in% c("1", "TRUE", "true", "T", "yes", "YES")

if (is.na(n_jobs) || n_jobs < 1L) {
  stop("N_READ_DEPTH_JOBS must be >= 1.")
}

if (is.na(max_hist_depth) || max_hist_depth < 20L) {
  stop("MAX_HIST_DEPTH must be >= 20.")
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1L) {
  stop("Pass array task ID as the first argument.")
}

job_id <- as.integer(args[1L])

if (is.na(job_id) || job_id < 1L || job_id > n_jobs) {
  stop(
    "Job ID must be between 1 and ",
    n_jobs,
    "."
  )
}

# -------------------------------- directories --------------------------------

partial_root <- file.path(output_root, "partial")
cpg_dir <- file.path(partial_root, "cpg")
subject_dir <- file.path(partial_root, "subject")
hist_dir <- file.path(partial_root, "subject_hist")
chr_dir <- file.path(partial_root, "subject_chr")
qc_dir <- file.path(partial_root, "qc")

for (d in c(
  output_root,
  partial_root,
  cpg_dir,
  subject_dir,
  hist_dir,
  chr_dir,
  qc_dir
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

cpg_file <- file.path(
  cpg_dir,
  sprintf(
    "read_depth_cpg_job_%03d.rds",
    job_id
  )
)

subject_file <- file.path(
  subject_dir,
  sprintf(
    "read_depth_subject_partial_job_%03d.tsv",
    job_id
  )
)

hist_file <- file.path(
  hist_dir,
  sprintf(
    "read_depth_subject_hist_job_%03d.rds",
    job_id
  )
)

chr_file <- file.path(
  chr_dir,
  sprintf(
    "read_depth_subject_chr_partial_job_%03d.tsv",
    job_id
  )
)

qc_file <- file.path(
  qc_dir,
  sprintf(
    "read_depth_qc_job_%03d.tsv",
    job_id
  )
)

outs <- c(
  cpg_file,
  subject_file,
  hist_file,
  chr_file,
  qc_file
)

if (!overwrite && all(file.exists(outs))) {
  cat(
    "All outputs already exist for job ",
    job_id,
    "; exiting.\n",
    sep = ""
  )
  quit(
    save = "no",
    status = 0L
  )
}

if (overwrite) {
  unlink(
    outs,
    force = TRUE
  )
}

run_start <- Sys.time()

# ---------------------------------- helpers ----------------------------------

normalize_id <- function(x) {
  x <- as.character(x)
  x <- sub("^X", "", x)
  x <- gsub("\\.", "-", x)
  x
}

normalize_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- sub("\\.0$", "", x)
  x
}

safe_cv <- function(sd_value, mean_value) {
  out <- rep(
    NA_real_,
    length(mean_value)
  )

  ok <- is.finite(sd_value) &
    is.finite(mean_value) &
    mean_value > 0

  out[ok] <- sd_value[ok] /
    mean_value[ok]

  out
}

safe_fraction <- function(num, den) {
  out <- rep(
    NA_real_,
    length(den)
  )

  ok <- is.finite(den) &
    den > 0

  out[ok] <- num[ok] /
    den[ok]

  out
}

safe_col_sums <- function(x) {
  colSums(
    x,
    na.rm = TRUE
  )
}

depth_to_hist <- function(v, max_depth) {
  v <- as.numeric(v)
  v <- v[
    is.finite(v) &
      v >= 0
  ]

  if (!length(v)) {
    return(
      integer(max_depth + 2L)
    )
  }

  # Read depth should be integer. Round only for histogram bookkeeping;
  # all means/SDs use the original numeric depth values.
  vi <- as.integer(round(v))
  vi <- pmin(
    vi,
    max_depth + 1L
  )

  tabulate(
    vi + 1L,
    nbins = max_depth + 2L
  )
}

# -------------------------------- phenotype ----------------------------------

if (!file.exists(pheno_file_path)) {
  stop(
    "Phenotype file does not exist: ",
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
    "pheno_file was not found in ",
    pheno_file_path
  )
}

ph <- as.data.table(
  copy(
    get(
      "pheno_file",
      envir = pheno_env
    )
  )
)

rm(pheno_env)

if (!"ID" %in% names(ph)) {
  stop("pheno_file lacks ID.")
}

ph[
  ,
  ID := normalize_id(ID)
]

ph <- unique(
  ph,
  by = "ID"
)

analytic_ids <- ph$ID

if (!length(analytic_ids)) {
  stop("No analytic phenotype IDs found.")
}

cat(
  "Analytic phenotype IDs: ",
  length(analytic_ids),
  "\n",
  sep = ""
)

# ------------------------------- chunk assignment -----------------------------

chunk_files <- list.files(
  chunk_dir,
  pattern = "^chunk_[0-9]+\\.csv$",
  full.names = TRUE
)

if (!length(chunk_files)) {
  stop(
    "No chunk_*.csv files found under ",
    chunk_dir
  )
}

chunk_ids <- suppressWarnings(
  as.integer(
    sub(
      "^chunk_([0-9]+)\\.csv$",
      "\\1",
      basename(chunk_files)
    )
  )
)

ord <- order(
  chunk_ids,
  basename(chunk_files)
)

chunk_files <- chunk_files[ord]
chunk_ids <- chunk_ids[ord]

assignment <- (
  (
    seq_along(chunk_files) - 1L
  ) %% n_jobs
) + 1L

my_idx <- which(
  assignment == job_id
)

my_files <- chunk_files[my_idx]
my_chunk_ids <- chunk_ids[my_idx]

cat(
  "Job ",
  job_id,
  " of ",
  n_jobs,
  ": assigned ",
  length(my_files),
  " chunks.\n",
  sep = ""
)

# ------------------------- subject-level accumulators -------------------------

subject_acc <- data.table(
  ID = analytic_ids,
  n_cpg_rows = 0,
  n_all_depth = 0,
  sum_depth_all = 0,
  sumsq_depth_all = 0,
  n_depth_lt5_all = 0,
  n_depth_ge5_all = 0,
  n_depth_ge10_all = 0,
  n_depth_ge20_all = 0,
  n_depth_ge30_all = 0,
  n_depth_ge50_all = 0,
  n_gam_observed = 0,
  sum_depth_gam = 0,
  sumsq_depth_gam = 0,
  n_depth_5_9_gam = 0,
  n_depth_10_19_gam = 0,
  n_depth_20_29_gam = 0,
  n_depth_30_49_gam = 0,
  n_depth_ge50_gam = 0
)

n_subjects <- nrow(
  subject_acc
)

hist_all <- matrix(
  0L,
  nrow = n_subjects,
  ncol = max_hist_depth + 2L,
  dimnames = list(
    analytic_ids,
    c(
      as.character(
        0:max_hist_depth
      ),
      paste0(
        ">",
        max_hist_depth
      )
    )
  )
)

hist_gam <- hist_all

cpg_list <- vector(
  "list",
  length(my_files)
)

chr_list <- list()

qc_list <- vector(
  "list",
  length(my_files)
)

# -------------------------------- main loop ----------------------------------

for (ii in seq_along(my_files)) {
  ff <- my_files[ii]
  cid <- my_chunk_ids[ii]

  cat(
    "\nChunk ",
    cid,
    " [",
    ii,
    "/",
    length(my_files),
    "]\n",
    sep = ""
  )

  header <- tryCatch(
    fread(
      ff,
      nrows = 0L,
      showProgress = FALSE
    ),
    error = function(e) e
  )

  if (inherits(header, "error")) {
    qc_list[[ii]] <- data.table(
      job_id = job_id,
      chunk_id = cid,
      chunk_file = ff,
      n_cpg_rows = NA_integer_,
      n_sample_pairs = NA_integer_,
      n_noninteger_depth_values = NA_real_,
      status = "HEADER_ERROR",
      message = conditionMessage(header)
    )
    next
  }

  hn <- names(header)

  if (!all(c("chr", "start") %in% hn)) {
    qc_list[[ii]] <- data.table(
      job_id = job_id,
      chunk_id = cid,
      chunk_file = ff,
      n_cpg_rows = NA_integer_,
      n_sample_pairs = NA_integer_,
      n_noninteger_depth_values = NA_real_,
      status = "MISSING_COORDINATES",
      message = "chr/start columns are required."
    )
    next
  }

  tot_cols <- grep(
    "_tot$",
    hn,
    value = TRUE
  )

  meth_cols <- grep(
    "_meth$",
    hn,
    value = TRUE
  )

  if (!length(tot_cols)) {
    qc_list[[ii]] <- data.table(
      job_id = job_id,
      chunk_id = cid,
      chunk_file = ff,
      n_cpg_rows = NA_integer_,
      n_sample_pairs = 0L,
      n_noninteger_depth_values = NA_real_,
      status = "NO_TOT_COLUMNS",
      message = "No *_tot read-depth columns found."
    )
    next
  }

  tot_map <- data.table(
    ID = normalize_id(
      sub(
        "_tot$",
        "",
        tot_cols
      )
    ),
    tot_col = tot_cols
  )

  meth_map <- data.table(
    ID = normalize_id(
      sub(
        "_meth$",
        "",
        meth_cols
      )
    ),
    meth_col = meth_cols
  )

  if (anyDuplicated(tot_map$ID)) {
    stop(
      "Duplicated normalized *_tot IDs in chunk ",
      cid
    )
  }

  if (anyDuplicated(meth_map$ID)) {
    stop(
      "Duplicated normalized *_meth IDs in chunk ",
      cid
    )
  }

  pair_map <- merge(
    tot_map,
    meth_map,
    by = "ID",
    all.x = TRUE
  )

  if (analytic_only) {
    pair_map <- pair_map[
      ID %chin% analytic_ids
    ]
  }

  if (!nrow(pair_map)) {
    qc_list[[ii]] <- data.table(
      job_id = job_id,
      chunk_id = cid,
      chunk_file = ff,
      n_cpg_rows = NA_integer_,
      n_sample_pairs = 0L,
      n_noninteger_depth_values = NA_real_,
      status = "NO_ANALYTIC_SAMPLE_PAIRS",
      message = NA_character_
    )
    next
  }

  # Force the analytic subject ordering so all matrices line up.
  pair_map <- pair_map[
    match(
      analytic_ids,
      ID
    )
  ]

  if (anyNA(pair_map$ID)) {
    missing_ids <- analytic_ids[
      is.na(pair_map$ID)
    ]

    stop(
      "Chunk ",
      cid,
      " is missing *_tot columns for ",
      length(missing_ids),
      " analytic subjects; examples: ",
      paste(
        head(
          missing_ids,
          10L
        ),
        collapse = ", "
      )
    )
  }

  select_cols <- unique(
    c(
      "chr",
      "start",
      intersect(
        "end",
        hn
      ),
      pair_map$tot_col,
      pair_map$meth_col[
        !is.na(
          pair_map$meth_col
        )
      ]
    )
  )

  dt <- tryCatch(
    fread(
      ff,
      select = select_cols,
      showProgress = FALSE
    ),
    error = function(e) e
  )

  if (inherits(dt, "error")) {
    qc_list[[ii]] <- data.table(
      job_id = job_id,
      chunk_id = cid,
      chunk_file = ff,
      n_cpg_rows = NA_integer_,
      n_sample_pairs = nrow(pair_map),
      n_noninteger_depth_values = NA_real_,
      status = "READ_ERROR",
      message = conditionMessage(dt)
    )
    next
  }

  dt[
    ,
    chr := normalize_chr(chr)
  ]

  dt[
    ,
    start := suppressWarnings(
      as.integer(start)
    )
  ]

  # CpG x subject total read-depth matrix.
  selected_tot_cols <- pair_map$tot_col

  Tmat <- as.matrix(
    dt[
      ,
      ..selected_tot_cols
    ]
  )

  storage.mode(Tmat) <- "numeric"

  colnames(Tmat) <- pair_map$ID

  # Paired methylation matrix. Missing *_meth columns are represented as NA.
  Mmat <- matrix(
    NA_real_,
    nrow = nrow(dt),
    ncol = n_subjects,
    dimnames = list(
      NULL,
      analytic_ids
    )
  )

  have_meth <- !is.na(
    pair_map$meth_col
  )

  if (any(have_meth)) {
    selected_meth_cols <- pair_map$meth_col[
      have_meth
    ]

    tmp_m <- as.matrix(
      dt[
        ,
        ..selected_meth_cols
      ]
    )

    storage.mode(tmp_m) <- "numeric"

    Mmat[
      ,
      which(have_meth)
    ] <- tmp_m
  }

  all_valid <- is.finite(Tmat) &
    Tmat >= 0

  meth_valid <- is.finite(Mmat) &
    Mmat >= 0 &
    Mmat <= 1

  # Primary GAM-DMR preprocessing states that CpG-subject measurements with
  # read depth < 5 are excluded. Enforce that threshold explicitly here when
  # defining the observations whose depth distribution is relevant to the GAM.
  gam_valid <- all_valid &
    Tmat >= 5 &
    meth_valid

  # QC diagnostic: finite methylation values that nevertheless have depth < 5.
  # Ideally this is zero if the upstream low-depth masking was already applied.
  finite_meth_depth_lt5_n <- sum(
    meth_valid &
      all_valid &
      Tmat < 5,
    na.rm = TRUE
  )

  # Diagnostic only: read depth should be integer-valued.
  noninteger_n <- sum(
    all_valid &
      abs(
        Tmat - round(Tmat)
      ) > 1e-8,
    na.rm = TRUE
  )

  Tall <- Tmat
  Tall[!all_valid] <- NA_real_

  Tgam <- Tmat
  Tgam[!gam_valid] <- NA_real_

  n_cpg <- nrow(dt)

  # ----------------------------- per-CpG stats -------------------------------

  n_all <- rowSums(
    is.finite(Tall)
  )

  n_gam <- rowSums(
    is.finite(Tgam)
  )

  mean_all <- matrixStats::rowMeans2(
    Tall,
    na.rm = TRUE
  )

  mean_gam <- matrixStats::rowMeans2(
    Tgam,
    na.rm = TRUE
  )

  sd_all <- matrixStats::rowSds(
    Tall,
    na.rm = TRUE
  )

  sd_gam <- matrixStats::rowSds(
    Tgam,
    na.rm = TRUE
  )

  q_all <- matrixStats::rowQuantiles(
    Tall,
    probs = c(
      0.25,
      0.50,
      0.75
    ),
    na.rm = TRUE
  )

  q_gam <- matrixStats::rowQuantiles(
    Tgam,
    probs = c(
      0.10,
      0.25,
      0.50,
      0.75,
      0.90,
      0.95
    ),
    na.rm = TRUE
  )

  cpg_dt <- data.table(
    job_id = job_id,
    data_chunk_id = cid,
    chr = dt$chr,
    start = dt$start,
    CpG = paste0(
      dt$chr,
      ":",
      dt$start
    ),
    n_analytic_subjects = n_subjects,
    n_all_depth = n_all,
    all_depth_fraction = n_all /
      n_subjects,
    mean_depth_all = mean_all,
    q25_depth_all = q_all[, 1L],
    median_depth_all = q_all[, 2L],
    q75_depth_all = q_all[, 3L],
    sd_depth_all = sd_all,
    cv_depth_all = safe_cv(
      sd_all,
      mean_all
    ),
    n_gam_observed = n_gam,
    call_rate_gam = n_gam /
      n_subjects,
    mean_depth_gam = mean_gam,
    p10_depth_gam = q_gam[, 1L],
    q25_depth_gam = q_gam[, 2L],
    median_depth_gam = q_gam[, 3L],
    q75_depth_gam = q_gam[, 4L],
    p90_depth_gam = q_gam[, 5L],
    p95_depth_gam = q_gam[, 6L],
    sd_depth_gam = sd_gam,
    cv_depth_gam = safe_cv(
      sd_gam,
      mean_gam
    ),
    iqr_depth_gam = q_gam[, 4L] -
      q_gam[, 2L]
  )

  cpg_dt[
    ,
    frac_depth_lt5_all :=
      safe_fraction(
        rowSums(
          Tall < 5,
          na.rm = TRUE
        ),
        n_all
      )
  ]

  cpg_dt[
    ,
    frac_depth_5_9_gam :=
      safe_fraction(
        rowSums(
          Tgam >= 5 &
            Tgam < 10,
          na.rm = TRUE
        ),
        n_gam
      )
  ]

  cpg_dt[
    ,
    frac_depth_10_19_gam :=
      safe_fraction(
        rowSums(
          Tgam >= 10 &
            Tgam < 20,
          na.rm = TRUE
        ),
        n_gam
      )
  ]

  cpg_dt[
    ,
    frac_depth_ge20_gam :=
      safe_fraction(
        rowSums(
          Tgam >= 20,
          na.rm = TRUE
        ),
        n_gam
      )
  ]

  cpg_list[[ii]] <- cpg_dt

  # --------------------------- per-subject stats -----------------------------

  subject_acc[
    ,
    n_cpg_rows :=
      n_cpg_rows +
      n_cpg
  ]

  subject_acc[
    ,
    n_all_depth :=
      n_all_depth +
      colSums(
        all_valid
      )
  ]

  subject_acc[
    ,
    sum_depth_all :=
      sum_depth_all +
      safe_col_sums(Tall)
  ]

  subject_acc[
    ,
    sumsq_depth_all :=
      sumsq_depth_all +
      safe_col_sums(Tall^2)
  ]

  subject_acc[
    ,
    n_depth_lt5_all :=
      n_depth_lt5_all +
      colSums(
        Tall < 5,
        na.rm = TRUE
      )
  ]

  subject_acc[
    ,
    n_depth_ge5_all :=
      n_depth_ge5_all +
      colSums(
        Tall >= 5,
        na.rm = TRUE
      )
  ]

  subject_acc[
    ,
    n_depth_ge10_all :=
      n_depth_ge10_all +
      colSums(
        Tall >= 10,
        na.rm = TRUE
      )
  ]

  subject_acc[
    ,
    n_depth_ge20_all :=
      n_depth_ge20_all +
      colSums(
        Tall >= 20,
        na.rm = TRUE
      )
  ]

  subject_acc[
    ,
    n_depth_ge30_all :=
      n_depth_ge30_all +
      colSums(
        Tall >= 30,
        na.rm = TRUE
      )
  ]

  subject_acc[
    ,
    n_depth_ge50_all :=
      n_depth_ge50_all +
      colSums(
        Tall >= 50,
        na.rm = TRUE
      )
  ]

  subject_acc[
    ,
    n_gam_observed :=
      n_gam_observed +
      colSums(
        gam_valid
      )
  ]

  subject_acc[
    ,
    sum_depth_gam :=
      sum_depth_gam +
      safe_col_sums(Tgam)
  ]

  subject_acc[
    ,
    sumsq_depth_gam :=
      sumsq_depth_gam +
      safe_col_sums(Tgam^2)
  ]

  subject_acc[
    ,
    n_depth_5_9_gam :=
      n_depth_5_9_gam +
      colSums(
        Tgam >= 5 &
          Tgam < 10,
        na.rm = TRUE
      )
  ]

  subject_acc[
    ,
    n_depth_10_19_gam :=
      n_depth_10_19_gam +
      colSums(
        Tgam >= 10 &
          Tgam < 20,
        na.rm = TRUE
      )
  ]

  subject_acc[
    ,
    n_depth_20_29_gam :=
      n_depth_20_29_gam +
      colSums(
        Tgam >= 20 &
          Tgam < 30,
        na.rm = TRUE
      )
  ]

  subject_acc[
    ,
    n_depth_30_49_gam :=
      n_depth_30_49_gam +
      colSums(
        Tgam >= 30 &
          Tgam < 50,
        na.rm = TRUE
      )
  ]

  subject_acc[
    ,
    n_depth_ge50_gam :=
      n_depth_ge50_gam +
      colSums(
        Tgam >= 50,
        na.rm = TRUE
      )
  ]

  # ----------------------- exact subject depth histograms --------------------

  for (jj in seq_len(n_subjects)) {
    hist_all[jj, ] <- hist_all[jj, ] +
      depth_to_hist(
        Tall[, jj],
        max_hist_depth
      )

    hist_gam[jj, ] <- hist_gam[jj, ] +
      depth_to_hist(
        Tgam[, jj],
        max_hist_depth
      )
  }

  # -------------------------- subject x chromosome ---------------------------

  chr_values <- unique(
    dt$chr
  )

  for (cc in chr_values) {
    ridx <- which(
      dt$chr == cc
    )

    if (!length(ridx)) {
      next
    }

    ta <- Tall[
      ridx,
      ,
      drop = FALSE
    ]

    tg <- Tgam[
      ridx,
      ,
      drop = FALSE
    ]

    chr_list[[length(chr_list) + 1L]] <- data.table(
      ID = analytic_ids,
      chr = cc,
      n_cpg_rows = length(ridx),
      n_all_depth = colSums(
        is.finite(ta)
      ),
      sum_depth_all = safe_col_sums(ta),
      sumsq_depth_all = safe_col_sums(
        ta^2
      ),
      n_gam_observed = colSums(
        is.finite(tg)
      ),
      sum_depth_gam = safe_col_sums(tg),
      sumsq_depth_gam = safe_col_sums(
        tg^2
      )
    )
  }

  qc_list[[ii]] <- data.table(
    job_id = job_id,
    chunk_id = cid,
    chunk_file = ff,
    n_cpg_rows = n_cpg,
    n_sample_pairs = nrow(pair_map),
    n_noninteger_depth_values = noninteger_n,
    n_finite_meth_depth_lt5 = finite_meth_depth_lt5_n,
    status = "COMPLETED",
    message = NA_character_
  )

  rm(
    dt,
    Tmat,
    Mmat,
    Tall,
    Tgam,
    all_valid,
    gam_valid,
    cpg_dt
  )

  invisible(
    gc(FALSE)
  )
}

# --------------------------------- outputs -----------------------------------

cpg_out <- rbindlist(
  cpg_list,
  use.names = TRUE,
  fill = TRUE
)

if (nrow(cpg_out)) {
  setorder(
    cpg_out,
    data_chunk_id,
    chr,
    start
  )
}

saveRDS(
  cpg_out,
  cpg_file,
  compress = "gzip"
)

fwrite(
  subject_acc,
  subject_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

saveRDS(
  list(
    job_id = job_id,
    IDs = analytic_ids,
    max_hist_depth = max_hist_depth,
    hist_all = hist_all,
    hist_gam = hist_gam
  ),
  hist_file,
  compress = "gzip"
)

chr_out <- if (length(chr_list)) {
  rbindlist(
    chr_list,
    use.names = TRUE,
    fill = TRUE
  )[
    ,
    .(
      n_cpg_rows = sum(n_cpg_rows),
      n_all_depth = sum(n_all_depth),
      sum_depth_all = sum(sum_depth_all),
      sumsq_depth_all = sum(sumsq_depth_all),
      n_gam_observed = sum(n_gam_observed),
      sum_depth_gam = sum(sum_depth_gam),
      sumsq_depth_gam = sum(sumsq_depth_gam)
    ),
    by = .(
      ID,
      chr
    )
  ]
} else {
  data.table()
}

fwrite(
  chr_out,
  chr_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

qc_out <- rbindlist(
  qc_list,
  use.names = TRUE,
  fill = TRUE
)

final_qc <- data.table(
  job_id = job_id,
  n_jobs = n_jobs,
  n_chunk_files_total = length(chunk_files),
  n_chunks_assigned = length(my_files),
  n_chunks_completed = sum(
    qc_out$status == "COMPLETED",
    na.rm = TRUE
  ),
  n_cpg_rows_processed = sum(
    qc_out$n_cpg_rows[
      qc_out$status == "COMPLETED"
    ],
    na.rm = TRUE
  ),
  n_analytic_subjects = n_subjects,
  max_hist_depth = max_hist_depth,
  elapsed_minutes = as.numeric(
    difftime(
      Sys.time(),
      run_start,
      units = "mins"
    )
  ),
  completed_at = format(
    Sys.time(),
    tz = "America/Toronto",
    usetz = TRUE
  ),
  job_status = if (
    length(my_files) ==
      sum(
        qc_out$status == "COMPLETED",
        na.rm = TRUE
      )
  ) {
    "COMPLETED"
  } else {
    "COMPLETED_WITH_WARNINGS"
  }
)

fwrite(
  rbindlist(
    list(
      qc_out,
      final_qc
    ),
    use.names = TRUE,
    fill = TRUE
  ),
  qc_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat(
  "\nRead-depth job completed.\n"
)

print(final_qc)

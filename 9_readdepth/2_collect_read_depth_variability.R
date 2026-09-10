#!/usr/bin/env Rscript

# =============================================================================
# 21_collect_read_depth_variability.R
#
# Collect chunk-wise read-depth QC outputs and produce:
#   - subject-level read-depth summaries;
#   - subject x chromosome summaries;
#   - genome-wide CpG read-depth variability summaries;
#   - 1-Mb genomic-bin summaries;
#   - observation-level read-depth histogram;
#   - manuscript/rebuttal-ready text summary;
#   - figures for read-depth variability across subjects and genome.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

data.table::setDTthreads(1L)

# ------------------------------- configuration -------------------------------

PATH_wk <- path.expand(
  Sys.getenv("PATH_WK", "~/scratch/UQAC/meth/")
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

allow_incomplete <- Sys.getenv(
  "ALLOW_INCOMPLETE", "0"
) %in% c("1", "TRUE", "true", "T", "yes", "YES")

partial_root <- file.path(output_root, "partial")
cpg_dir <- file.path(partial_root, "cpg")
subject_dir <- file.path(partial_root, "subject")
hist_dir <- file.path(partial_root, "subject_hist")
chr_dir <- file.path(partial_root, "subject_chr")
qc_dir <- file.path(partial_root, "qc")

table_dir <- file.path(output_root, "tables")
figure_dir <- file.path(output_root, "figures")

for (d in c(
  output_root,
  table_dir,
  figure_dir
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
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

chr_order_value <- function(x) {
  x <- normalize_chr(x)

  out <- suppressWarnings(
    as.integer(x)
  )

  out[
    x == "X"
  ] <- 23L

  out[
    x == "Y"
  ] <- 24L

  out[
    x %in% c("M", "MT")
  ] <- 25L

  out
}

variance_from_sums <- function(n, s, ss) {
  out <- rep(
    NA_real_,
    length(n)
  )

  ok <- is.finite(n) &
    n > 1 &
    is.finite(s) &
    is.finite(ss)

  vv <- (
    ss[ok] -
      (
        s[ok]^2 /
          n[ok]
      )
  ) /
    (
      n[ok] - 1
    )

  vv[
    vv < 0 &
      vv > -1e-8
  ] <- 0

  out[ok] <- sqrt(
    pmax(
      vv,
      0
    )
  )

  out
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

hist_quantile <- function(counts, p, max_depth) {
  counts <- as.numeric(counts)

  n <- sum(
    counts,
    na.rm = TRUE
  )

  if (!is.finite(n) || n <= 0) {
    return(NA_real_)
  }

  target <- max(
    1,
    ceiling(
      p * n
    )
  )

  idx <- which(
    cumsum(counts) >= target
  )[1L]

  if (!length(idx) || is.na(idx)) {
    return(NA_real_)
  }

  # Last bin is > max_depth.
  if (idx == max_depth + 2L) {
    return(NA_real_)
  }

  idx - 1L
}

summarize_variable <- function(x) {
  x <- as.numeric(x)
  x <- x[
    is.finite(x)
  ]

  if (!length(x)) {
    return(
      data.table(
        n = 0L,
        mean = NA_real_,
        sd = NA_real_,
        min = NA_real_,
        q025 = NA_real_,
        q25 = NA_real_,
        median = NA_real_,
        q75 = NA_real_,
        q975 = NA_real_,
        max = NA_real_
      )
    )
  }

  q <- quantile(
    x,
    probs = c(
      0.025,
      0.25,
      0.50,
      0.75,
      0.975
    ),
    names = FALSE,
    na.rm = TRUE
  )

  data.table(
    n = length(x),
    mean = mean(x),
    sd = sd(x),
    min = min(x),
    q025 = q[1L],
    q25 = q[2L],
    median = q[3L],
    q75 = q[4L],
    q975 = q[5L],
    max = max(x)
  )
}

fmt <- function(x, digits = 2L) {
  if (!length(x) || !is.finite(x[1L])) {
    return("NA")
  }

  sprintf(
    paste0(
      "%.",
      digits,
      "f"
    ),
    x[1L]
  )
}

fmt_pct <- function(x, digits = 1L) {
  if (!length(x) || !is.finite(x[1L])) {
    return("NA")
  }

  paste0(
    sprintf(
      paste0(
        "%.",
        digits,
        "f"
      ),
      100 * x[1L]
    ),
    "%"
  )
}

# ----------------------------- validate inputs -------------------------------

expected <- seq_len(
  n_jobs
)

subject_files <- file.path(
  subject_dir,
  sprintf(
    "read_depth_subject_partial_job_%03d.tsv",
    expected
  )
)

hist_files <- file.path(
  hist_dir,
  sprintf(
    "read_depth_subject_hist_job_%03d.rds",
    expected
  )
)

chr_files <- file.path(
  chr_dir,
  sprintf(
    "read_depth_subject_chr_partial_job_%03d.tsv",
    expected
  )
)

cpg_files <- file.path(
  cpg_dir,
  sprintf(
    "read_depth_cpg_job_%03d.rds",
    expected
  )
)

qc_files <- file.path(
  qc_dir,
  sprintf(
    "read_depth_qc_job_%03d.tsv",
    expected
  )
)

all_files <- list(
  subject = subject_files,
  hist = hist_files,
  chr = chr_files,
  cpg = cpg_files,
  qc = qc_files
)

for (nm in names(all_files)) {
  missing <- all_files[[nm]][
    !file.exists(
      all_files[[nm]]
    )
  ]

  if (
    length(missing) &&
      !allow_incomplete
  ) {
    stop(
      "Missing ",
      length(missing),
      " ",
      nm,
      " partial files."
    )
  }

  all_files[[nm]] <- all_files[[nm]][
    file.exists(
      all_files[[nm]]
    )
  ]
}

if (!length(all_files$subject)) {
  stop("No read-depth partial outputs found.")
}

# ------------------------------ phenotype data -------------------------------

pheno_env <- new.env(
  parent = emptyenv()
)

load(
  pheno_file_path,
  envir = pheno_env
)

ph <- as.data.table(
  copy(
    get(
      "pheno_file",
      envir = pheno_env
    )
  )
)

rm(pheno_env)

keep_pheno <- intersect(
  c(
    "ID",
    "FID",
    "AA_only",
    "AgeCalc",
    "Sex",
    "BMI"
  ),
  names(ph)
)

ph <- unique(
  ph[
    ,
    ..keep_pheno
  ],
  by = "ID"
)

ph[
  ,
  ID := normalize_id(ID)
]

if ("FID" %in% names(ph)) {
  ph[
    ,
    FID := as.character(FID)
  ]
}

if ("AA_only" %in% names(ph)) {
  ph[
    ,
    AA_only := suppressWarnings(
      as.integer(
        as.character(AA_only)
      )
    )
  ]
}

# --------------------------- combine subject stats ---------------------------

subject_partial <- rbindlist(
  lapply(
    all_files$subject,
    fread
  ),
  use.names = TRUE,
  fill = TRUE
)

subject <- subject_partial[
  ,
  .(
    n_cpg_rows = sum(n_cpg_rows),
    n_all_depth = sum(n_all_depth),
    sum_depth_all = sum(sum_depth_all),
    sumsq_depth_all = sum(sumsq_depth_all),
    n_depth_lt5_all = sum(n_depth_lt5_all),
    n_depth_ge5_all = sum(n_depth_ge5_all),
    n_depth_ge10_all = sum(n_depth_ge10_all),
    n_depth_ge20_all = sum(n_depth_ge20_all),
    n_depth_ge30_all = sum(n_depth_ge30_all),
    n_depth_ge50_all = sum(n_depth_ge50_all),
    n_gam_observed = sum(n_gam_observed),
    sum_depth_gam = sum(sum_depth_gam),
    sumsq_depth_gam = sum(sumsq_depth_gam),
    n_depth_5_9_gam = sum(n_depth_5_9_gam),
    n_depth_10_19_gam = sum(n_depth_10_19_gam),
    n_depth_20_29_gam = sum(n_depth_20_29_gam),
    n_depth_30_49_gam = sum(n_depth_30_49_gam),
    n_depth_ge50_gam = sum(n_depth_ge50_gam)
  ),
  by = ID
]

subject[
  ,
  mean_depth_all :=
    sum_depth_all /
    n_all_depth
]

subject[
  ,
  sd_depth_all :=
    variance_from_sums(
      n_all_depth,
      sum_depth_all,
      sumsq_depth_all
    )
]

subject[
  ,
  cv_depth_all :=
    safe_cv(
      sd_depth_all,
      mean_depth_all
    )
]

subject[
  ,
  mean_depth_gam :=
    sum_depth_gam /
    n_gam_observed
]

subject[
  ,
  sd_depth_gam :=
    variance_from_sums(
      n_gam_observed,
      sum_depth_gam,
      sumsq_depth_gam
    )
]

subject[
  ,
  cv_depth_gam :=
    safe_cv(
      sd_depth_gam,
      mean_depth_gam
    )
]

subject[
  ,
  gam_call_rate :=
    safe_fraction(
      n_gam_observed,
      n_cpg_rows
    )
]

subject[
  ,
  frac_depth_lt5_all :=
    safe_fraction(
      n_depth_lt5_all,
      n_all_depth
    )
]

subject[
  ,
  frac_depth_5_9_gam :=
    safe_fraction(
      n_depth_5_9_gam,
      n_gam_observed
    )
]

subject[
  ,
  frac_depth_10_19_gam :=
    safe_fraction(
      n_depth_10_19_gam,
      n_gam_observed
    )
]

subject[
  ,
  frac_depth_20_29_gam :=
    safe_fraction(
      n_depth_20_29_gam,
      n_gam_observed
    )
]

subject[
  ,
  frac_depth_30_49_gam :=
    safe_fraction(
      n_depth_30_49_gam,
      n_gam_observed
    )
]

subject[
  ,
  frac_depth_ge50_gam :=
    safe_fraction(
      n_depth_ge50_gam,
      n_gam_observed
    )
]

# -------------------------- aggregate subject histograms ---------------------

hist_objects <- lapply(
  all_files$hist,
  readRDS
)

ref_ids <- hist_objects[[1L]]$IDs
max_hist_depth <- hist_objects[[1L]]$max_hist_depth

hist_all <- matrix(
  0,
  nrow = length(ref_ids),
  ncol = max_hist_depth + 2L,
  dimnames = list(
    ref_ids,
    colnames(
      hist_objects[[1L]]$hist_all
    )
  )
)

hist_gam <- hist_all

for (hh in hist_objects) {
  if (!identical(
    as.character(hh$IDs),
    as.character(ref_ids)
  )) {
    stop("Subject histogram ID ordering differs across jobs.")
  }

  if (hh$max_hist_depth != max_hist_depth) {
    stop("MAX_HIST_DEPTH differs across jobs.")
  }

  hist_all <- hist_all +
    hh$hist_all

  hist_gam <- hist_gam +
    hh$hist_gam
}

hist_id_index <- match(
  subject$ID,
  rownames(hist_gam)
)

if (anyNA(hist_id_index)) {
  stop("Could not match all subjects to histogram matrices.")
}

subject[
  ,
  q25_depth_gam := vapply(
    hist_id_index,
    function(i) {
      hist_quantile(
        hist_gam[i, ],
        0.25,
        max_hist_depth
      )
    },
    numeric(1)
  )
]

subject[
  ,
  median_depth_gam := vapply(
    hist_id_index,
    function(i) {
      hist_quantile(
        hist_gam[i, ],
        0.50,
        max_hist_depth
      )
    },
    numeric(1)
  )
]

subject[
  ,
  q75_depth_gam := vapply(
    hist_id_index,
    function(i) {
      hist_quantile(
        hist_gam[i, ],
        0.75,
        max_hist_depth
      )
    },
    numeric(1)
  )
]

subject[
  ,
  p90_depth_gam := vapply(
    hist_id_index,
    function(i) {
      hist_quantile(
        hist_gam[i, ],
        0.90,
        max_hist_depth
      )
    },
    numeric(1)
  )
]

subject <- merge(
  subject,
  ph,
  by = "ID",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  subject,
  ID
)

fwrite(
  subject,
  file.path(
    table_dir,
    "read_depth_subject_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ---------------------- subject-by-chromosome summaries ----------------------

chr_partial <- rbindlist(
  lapply(
    all_files$chr,
    fread
  ),
  use.names = TRUE,
  fill = TRUE
)

subject_chr <- chr_partial[
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

subject_chr[
  ,
  mean_depth_gam :=
    sum_depth_gam /
    n_gam_observed
]

subject_chr[
  ,
  sd_depth_gam :=
    variance_from_sums(
      n_gam_observed,
      sum_depth_gam,
      sumsq_depth_gam
    )
]

subject_chr[
  ,
  cv_depth_gam :=
    safe_cv(
      sd_depth_gam,
      mean_depth_gam
    )
]

subject_chr[
  ,
  gam_call_rate :=
    safe_fraction(
      n_gam_observed,
      n_cpg_rows
    )
]

fwrite(
  subject_chr,
  file.path(
    table_dir,
    "read_depth_subject_by_chromosome.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ----------------------------- CpG-level summaries ---------------------------

cat(
  "Loading CpG summary partials...\n"
)

cpg <- rbindlist(
  lapply(
    all_files$cpg,
    readRDS
  ),
  use.names = TRUE,
  fill = TRUE
)

if (!nrow(cpg)) {
  stop("No CpG-level read-depth summaries were collected.")
}

cpg[
  ,
  chr := normalize_chr(chr)
]

cpg[
  ,
  chr_order := chr_order_value(chr)
]

cpg[
  ,
  bin_1mb_start :=
    floor(
      start /
        1e6
    ) *
      1e6
]

cpg[
  ,
  bin_1mb_end :=
    bin_1mb_start +
      1e6 -
      1
]

cpg[
  ,
  bin_1mb_mid :=
    (
      bin_1mb_start +
        bin_1mb_end
    ) /
      2
]

cpg_variables <- c(
  "call_rate_gam",
  "mean_depth_gam",
  "median_depth_gam",
  "sd_depth_gam",
  "cv_depth_gam",
  "iqr_depth_gam",
  "frac_depth_5_9_gam",
  "frac_depth_10_19_gam",
  "frac_depth_ge20_gam"
)

cpg_summary <- rbindlist(
  lapply(
    cpg_variables,
    function(v) {
      z <- summarize_variable(
        cpg[[v]]
      )

      z[
        ,
        metric := v
      ]

      z
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

setcolorder(
  cpg_summary,
  c(
    "metric",
    setdiff(
      names(cpg_summary),
      "metric"
    )
  )
)

fwrite(
  cpg_summary,
  file.path(
    table_dir,
    "read_depth_cpg_variability_overall.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

genomic_bins <- cpg[
  ,
  .(
    n_cpgs = .N,
    mean_of_cpg_mean_depth = mean(
      mean_depth_gam,
      na.rm = TRUE
    ),
    median_of_cpg_median_depth = median(
      median_depth_gam,
      na.rm = TRUE
    ),
    mean_cpg_sd_depth = mean(
      sd_depth_gam,
      na.rm = TRUE
    ),
    mean_cpg_cv_depth = mean(
      cv_depth_gam,
      na.rm = TRUE
    ),
    median_cpg_cv_depth = median(
      cv_depth_gam,
      na.rm = TRUE
    ),
    mean_gam_call_rate = mean(
      call_rate_gam,
      na.rm = TRUE
    ),
    median_gam_call_rate = median(
      call_rate_gam,
      na.rm = TRUE
    ),
    fraction_cpg_median_depth_lt10 = mean(
      median_depth_gam < 10,
      na.rm = TRUE
    ),
    fraction_cpg_median_depth_ge20 = mean(
      median_depth_gam >= 20,
      na.rm = TRUE
    )
  ),
  by = .(
    chr,
    chr_order,
    bin_1mb_start,
    bin_1mb_end,
    bin_1mb_mid
  )
]

setorder(
  genomic_bins,
  chr_order,
  bin_1mb_start
)

fwrite(
  genomic_bins,
  file.path(
    table_dir,
    "read_depth_genomic_bins_1Mb.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

chromosome_summary <- cpg[
  ,
  .(
    n_cpgs = .N,
    mean_of_cpg_mean_depth = mean(
      mean_depth_gam,
      na.rm = TRUE
    ),
    median_of_cpg_median_depth = median(
      median_depth_gam,
      na.rm = TRUE
    ),
    mean_cpg_cv_depth = mean(
      cv_depth_gam,
      na.rm = TRUE
    ),
    median_cpg_cv_depth = median(
      cv_depth_gam,
      na.rm = TRUE
    ),
    mean_gam_call_rate = mean(
      call_rate_gam,
      na.rm = TRUE
    ),
    median_gam_call_rate = median(
      call_rate_gam,
      na.rm = TRUE
    )
  ),
  by = .(
    chr,
    chr_order
  )
]

setorder(
  chromosome_summary,
  chr_order
)

fwrite(
  chromosome_summary,
  file.path(
    table_dir,
    "read_depth_by_chromosome.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# Keep a reproducible 100,000-CpG sample for inspection/plotting without
# producing another multi-million-row text file.
set.seed(20260812)

sample_n <- min(
  100000L,
  nrow(cpg)
)

cpg_sample <- cpg[
  sample(
    .N,
    sample_n
  )
]

fwrite(
  cpg_sample,
  file.path(
    table_dir,
    "read_depth_cpg_summary_random100k.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------- observation-level histogram -----------------------

obs_hist_all <- colSums(
  hist_all
)

obs_hist_gam <- colSums(
  hist_gam
)

depth_labels <- colnames(
  hist_gam
)

observation_hist <- data.table(
  depth_bin = depth_labels,
  all_depth_count = as.numeric(
    obs_hist_all
  ),
  gam_observed_count = as.numeric(
    obs_hist_gam
  )
)

observation_hist[
  ,
  all_depth_fraction :=
    all_depth_count /
    sum(all_depth_count)
]

observation_hist[
  ,
  gam_observed_fraction :=
    gam_observed_count /
    sum(gam_observed_count)
]

fwrite(
  observation_hist,
  file.path(
    table_dir,
    "read_depth_observation_histogram.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

overall_gam_mean <- sum(
  subject$sum_depth_gam
) /
  sum(
    subject$n_gam_observed
  )

overall_all_mean <- sum(
  subject$sum_depth_all
) /
  sum(
    subject$n_all_depth
  )

overall_gam_median <- hist_quantile(
  obs_hist_gam,
  0.50,
  max_hist_depth
)

overall_gam_q25 <- hist_quantile(
  obs_hist_gam,
  0.25,
  max_hist_depth
)

overall_gam_q75 <- hist_quantile(
  obs_hist_gam,
  0.75,
  max_hist_depth
)

overall_gam_p90 <- hist_quantile(
  obs_hist_gam,
  0.90,
  max_hist_depth
)

overall_gam_p95 <- hist_quantile(
  obs_hist_gam,
  0.95,
  max_hist_depth
)

depth_numeric <- 0:max_hist_depth

gam_counts_exact <- obs_hist_gam[
  seq_len(
    max_hist_depth + 1L
  )
]

gam_total <- sum(
  obs_hist_gam
)

gam_frac_5_9 <- sum(
  gam_counts_exact[
    depth_numeric >= 5 &
      depth_numeric < 10
  ]
) /
  gam_total

gam_frac_10_19 <- sum(
  gam_counts_exact[
    depth_numeric >= 10 &
      depth_numeric < 20
  ]
) /
  gam_total

gam_frac_20_29 <- sum(
  gam_counts_exact[
    depth_numeric >= 20 &
      depth_numeric < 30
  ]
) /
  gam_total

gam_frac_30_49 <- sum(
  gam_counts_exact[
    depth_numeric >= 30 &
      depth_numeric < 50
  ]
) /
  gam_total

gam_frac_ge50 <- (
  sum(
    gam_counts_exact[
      depth_numeric >= 50
    ]
  ) +
    obs_hist_gam[
      max_hist_depth + 2L
    ]
) /
  gam_total

# ------------------------ subject descriptive summaries ----------------------

subject_metrics <- c(
  "mean_depth_gam",
  "median_depth_gam",
  "sd_depth_gam",
  "cv_depth_gam",
  "gam_call_rate",
  "frac_depth_5_9_gam",
  "frac_depth_10_19_gam",
  "frac_depth_20_29_gam",
  "frac_depth_30_49_gam",
  "frac_depth_ge50_gam"
)

subject_summary <- rbindlist(
  lapply(
    subject_metrics,
    function(v) {
      z <- summarize_variable(
        subject[[v]]
      )

      z[
        ,
        metric := v
      ]

      z
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

setcolorder(
  subject_summary,
  c(
    "metric",
    setdiff(
      names(subject_summary),
      "metric"
    )
  )
)

fwrite(
  subject_summary,
  file.path(
    table_dir,
    "read_depth_subject_variability_overall.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

if ("AA_only" %in% names(subject)) {
  aa_summary <- subject[
    AA_only %in% c(0L, 1L),
    .(
      n_subjects = .N,
      mean_subject_mean_depth = mean(
        mean_depth_gam,
        na.rm = TRUE
      ),
      sd_subject_mean_depth = sd(
        mean_depth_gam,
        na.rm = TRUE
      ),
      median_subject_mean_depth = median(
        mean_depth_gam,
        na.rm = TRUE
      ),
      mean_subject_median_depth = mean(
        median_depth_gam,
        na.rm = TRUE
      ),
      mean_call_rate = mean(
        gam_call_rate,
        na.rm = TRUE
      ),
      mean_subject_cv = mean(
        cv_depth_gam,
        na.rm = TRUE
      )
    ),
    by = AA_only
  ]

  fwrite(
    aa_summary,
    file.path(
      table_dir,
      "read_depth_summary_by_AA_status_DESCRIPTIVE.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

# ---------------------------------- figures ----------------------------------

subject[
  ,
  subject_rank := frank(
    mean_depth_gam,
    ties.method = "first"
  )
]

p_subject <- ggplot(
  subject,
  aes(
    x = subject_rank,
    y = mean_depth_gam
  )
) +
  geom_point(
    size = 1.5,
    alpha = 0.75
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Read-depth variability across analytic subjects",
    x = "Subjects ranked by mean read depth",
    y = "Mean read depth across GAM-observed CpGs"
  )

ggsave(
  file.path(
    figure_dir,
    "Figure_read_depth_A_across_subjects.pdf"
  ),
  p_subject,
  width = 8.5,
  height = 5.5
)

ggsave(
  file.path(
    figure_dir,
    "Figure_read_depth_A_across_subjects.png"
  ),
  p_subject,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

p_cpg_median <- ggplot(
  cpg_sample[
    is.finite(
      median_depth_gam
    )
  ],
  aes(
    x = median_depth_gam
  )
) +
  geom_histogram(
    bins = 80
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Distribution of CpG-specific median read depth",
    x = "Median read depth across analytic subjects",
    y = "CpGs"
  )

ggsave(
  file.path(
    figure_dir,
    "Figure_read_depth_B_CpG_median_distribution.pdf"
  ),
  p_cpg_median,
  width = 8.5,
  height = 5.5
)

ggsave(
  file.path(
    figure_dir,
    "Figure_read_depth_B_CpG_median_distribution.png"
  ),
  p_cpg_median,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

p_cpg_cv <- ggplot(
  cpg_sample[
    is.finite(
      cv_depth_gam
    ) &
      cv_depth_gam >= 0
  ],
  aes(
    x = cv_depth_gam
  )
) +
  geom_histogram(
    bins = 80
  ) +
  coord_cartesian(
    xlim = quantile(
      cpg_sample$cv_depth_gam[
        is.finite(
          cpg_sample$cv_depth_gam
        )
      ],
      c(
        0,
        0.99
      ),
      na.rm = TRUE
    )
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "CpG-specific read-depth variability across subjects",
    subtitle = "x-axis truncated at the 99th percentile for visualization",
    x = "Coefficient of variation of read depth",
    y = "CpGs"
  )

ggsave(
  file.path(
    figure_dir,
    "Figure_read_depth_C_CpG_CV_distribution.pdf"
  ),
  p_cpg_cv,
  width = 8.5,
  height = 5.5
)

ggsave(
  file.path(
    figure_dir,
    "Figure_read_depth_C_CpG_CV_distribution.png"
  ),
  p_cpg_cv,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

genomic_bins_plot <- genomic_bins[
  is.finite(
    chr_order
  ) &
    is.finite(
      mean_of_cpg_mean_depth
    )
]

genomic_bins_plot[
  ,
  chr_factor := factor(
    chr,
    levels = chromosome_summary$chr[
      order(
        chromosome_summary$chr_order
      )
    ]
  )
]

p_genome <- ggplot(
  genomic_bins_plot,
  aes(
    x = bin_1mb_mid /
      1e6,
    y = mean_of_cpg_mean_depth,
    group = 1
  )
) +
  geom_line(
    linewidth = 0.45
  ) +
  facet_wrap(
    ~chr_factor,
    scales = "free_x",
    ncol = 4
  ) +
  theme_bw(
    base_size = 9
  ) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(
      size = 8
    )
  ) +
  labs(
    title = "Read-depth variability across the targeted genome",
    subtitle = "Mean CpG read depth summarized in 1-Mb genomic bins",
    x = "Genomic position (Mb)",
    y = "Mean read depth"
  )

ggsave(
  file.path(
    figure_dir,
    "Figure_read_depth_D_across_genome_1Mb.pdf"
  ),
  p_genome,
  width = 12,
  height = 12
)

ggsave(
  file.path(
    figure_dir,
    "Figure_read_depth_D_across_genome_1Mb.png"
  ),
  p_genome,
  width = 12,
  height = 12,
  dpi = 300
)

if ("AA_only" %in% names(subject)) {
  p_aa <- ggplot(
    subject[
      AA_only %in% c(0L, 1L)
    ],
    aes(
      x = factor(
        AA_only,
        levels = c(
          0,
          1
        ),
        labels = c(
          "Control",
          "AA"
        )
      ),
      y = mean_depth_gam
    )
  ) +
    geom_boxplot(
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.15,
      alpha = 0.5,
      size = 1.2
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Subject-level mean read depth by allergic asthma status",
      x = NULL,
      y = "Mean read depth across GAM-observed CpGs"
    )

  ggsave(
    file.path(
      figure_dir,
      "Figure_read_depth_E_by_AA_status_DESCRIPTIVE.pdf"
    ),
    p_aa,
    width = 6.5,
    height = 5.5
  )

  ggsave(
    file.path(
      figure_dir,
      "Figure_read_depth_E_by_AA_status_DESCRIPTIVE.png"
    ),
    p_aa,
    width = 6.5,
    height = 5.5,
    dpi = 300
  )
}

# ------------------------- manuscript-ready text summary ---------------------

subject_mean_summary <- summarize_variable(
  subject$mean_depth_gam
)

subject_median_summary <- summarize_variable(
  subject$median_depth_gam
)

subject_cv_summary <- summarize_variable(
  subject$cv_depth_gam
)

subject_call_summary <- summarize_variable(
  subject$gam_call_rate
)

cpg_median_summary <- summarize_variable(
  cpg$median_depth_gam
)

cpg_cv_summary <- summarize_variable(
  cpg$cv_depth_gam
)

cpg_call_summary <- summarize_variable(
  cpg$call_rate_gam
)

summary_lines <- c(
  "READ-DEPTH VARIABILITY ANALYSIS",
  "===============================",
  "",
  paste0(
    "Analytic subjects: ",
    nrow(subject)
  ),
  paste0(
    "CpGs evaluated: ",
    format(
      nrow(cpg),
      big.mark = ",",
      scientific = FALSE
    )
  ),
  "",
  "Observation-level read depth among GAM-observed methylation measurements",
  "-----------------------------------------------------------------------",
  paste0(
    "Mean read depth: ",
    fmt(
      overall_gam_mean,
      2
    )
  ),
  paste0(
    "Median read depth (IQR): ",
    fmt(
      overall_gam_median,
      0
    ),
    " (",
    fmt(
      overall_gam_q25,
      0
    ),
    "-",
    fmt(
      overall_gam_q75,
      0
    ),
    ")"
  ),
  paste0(
    "90th percentile: ",
    fmt(
      overall_gam_p90,
      0
    ),
    "; 95th percentile: ",
    fmt(
      overall_gam_p95,
      0
    )
  ),
  paste0(
    "Depth 5-9: ",
    fmt_pct(
      gam_frac_5_9
    )
  ),
  paste0(
    "Depth 10-19: ",
    fmt_pct(
      gam_frac_10_19
    )
  ),
  paste0(
    "Depth 20-29: ",
    fmt_pct(
      gam_frac_20_29
    )
  ),
  paste0(
    "Depth 30-49: ",
    fmt_pct(
      gam_frac_30_49
    )
  ),
  paste0(
    "Depth >=50: ",
    fmt_pct(
      gam_frac_ge50
    )
  ),
  "",
  "Variability across subjects",
  "---------------------------",
  paste0(
    "Subject mean read depth: median=",
    fmt(
      subject_mean_summary$median,
      2
    ),
    ", IQR=",
    fmt(
      subject_mean_summary$q25,
      2
    ),
    "-",
    fmt(
      subject_mean_summary$q75,
      2
    ),
    ", range=",
    fmt(
      subject_mean_summary$min,
      2
    ),
    "-",
    fmt(
      subject_mean_summary$max,
      2
    ),
    "."
  ),
  paste0(
    "Subject median read depth: median=",
    fmt(
      subject_median_summary$median,
      2
    ),
    ", IQR=",
    fmt(
      subject_median_summary$q25,
      2
    ),
    "-",
    fmt(
      subject_median_summary$q75,
      2
    ),
    "."
  ),
  paste0(
    "Subject read-depth CV: median=",
    fmt(
      subject_cv_summary$median,
      3
    ),
    ", IQR=",
    fmt(
      subject_cv_summary$q25,
      3
    ),
    "-",
    fmt(
      subject_cv_summary$q75,
      3
    ),
    "."
  ),
  paste0(
    "Subject GAM-observed CpG call rate: median=",
    fmt_pct(
      subject_call_summary$median
    ),
    ", IQR=",
    fmt_pct(
      subject_call_summary$q25
    ),
    "-",
    fmt_pct(
      subject_call_summary$q75
    ),
    "."
  ),
  "",
  "Variability across CpGs / genome",
  "--------------------------------",
  paste0(
    "CpG median read depth across subjects: median=",
    fmt(
      cpg_median_summary$median,
      2
    ),
    ", IQR=",
    fmt(
      cpg_median_summary$q25,
      2
    ),
    "-",
    fmt(
      cpg_median_summary$q75,
      2
    ),
    "."
  ),
  paste0(
    "CpG read-depth CV across subjects: median=",
    fmt(
      cpg_cv_summary$median,
      3
    ),
    ", IQR=",
    fmt(
      cpg_cv_summary$q25,
      3
    ),
    "-",
    fmt(
      cpg_cv_summary$q75,
      3
    ),
    "."
  ),
  paste0(
    "CpG GAM-observed call rate: median=",
    fmt_pct(
      cpg_call_summary$median
    ),
    ", IQR=",
    fmt_pct(
      cpg_call_summary$q25
    ),
    "-",
    fmt_pct(
      cpg_call_summary$q75
    ),
    "."
  ),
  "",
  paste0(
    "For comparison, the mean read depth across all finite *_tot observations, ",
    "including observations without a valid paired methylation proportion, was ",
    fmt(
      overall_all_mean,
      2
    ),
    "."
  ),
  "",
  "Interpretation note",
  "-------------------",
  paste0(
    "These diagnostics quantify heterogeneity in read-depth measurement precision ",
    "across CpGs/genomic locations and across subjects. They do not by themselves ",
    "replace a read-depth-weighted or count-based sensitivity analysis; rather, ",
    "they show the magnitude and structure of the read-depth variability that the ",
    "reviewer requested be characterized."
  ),
  ""
)

writeLines(
  summary_lines,
  con = file.path(
    output_root,
    "read_depth_results_summary.txt"
  )
)

# ------------------------------------ QC -------------------------------------

qc_all <- rbindlist(
  lapply(
    all_files$qc,
    fread
  ),
  use.names = TRUE,
  fill = TRUE
)

collection_qc <- data.table(
  expected_jobs = n_jobs,
  n_subject_files = length(
    all_files$subject
  ),
  n_hist_files = length(
    all_files$hist
  ),
  n_chr_files = length(
    all_files$chr
  ),
  n_cpg_files = length(
    all_files$cpg
  ),
  n_qc_files = length(
    all_files$qc
  ),
  n_analytic_subjects = nrow(subject),
  n_cpgs = nrow(cpg),
  n_chunk_rows_completed = qc_all[
    status == "COMPLETED",
    .N
  ],
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
    length(
      all_files$cpg
    ) == n_jobs
  ) {
    "COMPLETED"
  } else {
    "COMPLETED_INCOMPLETE_ALLOWED"
  }
)

fwrite(
  collection_qc,
  file.path(
    output_root,
    "read_depth_collection_qc.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat(
  "\nRead-depth variability collection completed.\n\n"
)

print(collection_qc)

cat(
  "\nKey summary:\n"
)

cat(
  paste(
    summary_lines[
      1:min(
        length(summary_lines),
        35L
      )
    ],
    collapse = "\n"
  ),
  "\n"
)

#!/usr/bin/env Rscript

# =============================================================================
# 23_collect_compare_read_depth_weighted_GAM.R
#
# Collect the read-depth-weighted GAM-DMR sensitivity results and compare them
# with the original unweighted results produced by 20_run_region_BMI.R.
#
# Primary comparison:
#   - p-value concordance for s(start):AA_only
#   - BH-FDR concordance
#   - EDF concordance
#   - mean/max absolute effect-size concordance
#   - overlap of final DMR calls
#
# Default final DMR definition:
#   BH-FDR < 0.05
#   EDF for s(start):AA_only > 0.5
#   mean_diff > 0.05
#
# Thresholds are configurable through environment variables.
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

weighted_dir <- path.expand(
  Sys.getenv(
    "WEIGHTED_GAM_OUTPUT",
    file.path(
      PATH_wk,
      "results/15_revision/7_read_depth_weighted_GAM"
    )
  )
)

# Original unweighted GAM-DMR results were produced in two batches.
# BMI_B1 contains the primary forward run; BMI_B2 contains the reverse/recovery
# run for regions that were not completed successfully in BMI_B1.
#
# Both directories are collected below. If the same region is present in both
# batches, BMI_B2 is given priority because it represents the recovery/rerun
# result for that region.
unweighted_dir_B1 <- path.expand(
  Sys.getenv(
    "UNWEIGHTED_GAM_OUTPUT_B1",
    file.path(
      PATH_wk,
      "results/11_mgcv/BMI_B1"
    )
  )
)

unweighted_dir_B2 <- path.expand(
  Sys.getenv(
    "UNWEIGHTED_GAM_OUTPUT_B2",
    file.path(
      PATH_wk,
      "results/11_mgcv/BMI_B2"
    )
  )
)

comparison_dir <- path.expand(
  Sys.getenv(
    "WEIGHTED_GAM_COMPARISON_OUTPUT",
    file.path(
      PATH_wk,
      "results/15_revision/7_read_depth_weighted_GAM/comparison"
    )
  )
)

fdr_cutoff <- as.numeric(
  Sys.getenv(
    "DMR_FDR_CUTOFF",
    "0.05"
  )
)

edf_cutoff <- as.numeric(
  Sys.getenv(
    "DMR_EDF_CUTOFF",
    "0.5"
  )
)

effect_cutoff <- as.numeric(
  Sys.getenv(
    "DMR_MEAN_DIFF_CUTOFF",
    "0.05"
  )
)

if (
  !is.finite(fdr_cutoff) ||
    fdr_cutoff <= 0 ||
    fdr_cutoff > 1
) {
  stop(
    "DMR_FDR_CUTOFF must be in (0,1]."
  )
}

if (
  !is.finite(edf_cutoff) ||
    edf_cutoff < 0
) {
  stop(
    "DMR_EDF_CUTOFF must be >= 0."
  )
}

if (
  !is.finite(effect_cutoff) ||
    effect_cutoff < 0
) {
  stop(
    "DMR_MEAN_DIFF_CUTOFF must be >= 0."
  )
}

table_dir <- file.path(
  comparison_dir,
  "tables"
)

figure_dir <- file.path(
  comparison_dir,
  "figures"
)

for (d in c(
  comparison_dir,
  table_dir,
  figure_dir
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ---------------------------------- helpers ----------------------------------

safe_cor <- function(x, y, method) {
  ok <- is.finite(x) &
    is.finite(y)

  if (sum(ok) < 3L) {
    return(NA_real_)
  }

  suppressWarnings(
    cor(
      x[ok],
      y[ok],
      method = method
    )
  )
}

fmt <- function(x, digits = 3L) {
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

safe_read_tsv <- function(f) {
  if (!file.exists(f)) {
    return(NULL)
  }

  info <- file.info(f)

  if (
    is.na(info$size) ||
      info$size == 0
  ) {
    return(NULL)
  }

  tryCatch(
    fread(
      f,
      showProgress = FALSE
    ),
    error = function(e) NULL
  )
}

safe_read_old <- function(f, old_names) {
  if (!file.exists(f)) {
    return(NULL)
  }

  info <- file.info(f)

  if (
    is.na(info$size) ||
      info$size == 0
  ) {
    return(NULL)
  }

  z <- tryCatch(
    fread(
      f,
      header = FALSE,
      showProgress = FALSE
    ),
    error = function(e) NULL
  )

  if (
    is.null(z) ||
      !nrow(z)
  ) {
    return(NULL)
  }

  if (ncol(z) != length(old_names)) {
    warning(
      "Skipping unweighted file with ",
      ncol(z),
      " columns (expected ",
      length(old_names),
      "): ",
      f
    )
    return(NULL)
  }

  setnames(
    z,
    old_names
  )

  z
}

# --------------------------- collect weighted results ------------------------

weighted_files <- list.files(
  weighted_dir,
  pattern = "^weighted_GAM_results_job_[0-9]+\\.tsv$",
  full.names = TRUE
)

if (!length(weighted_files)) {
  stop(
    "No weighted GAM result files found under ",
    weighted_dir
  )
}

cat(
  "Weighted result files: ",
  length(weighted_files),
  "\n",
  sep = ""
)

weighted <- rbindlist(
  lapply(
    weighted_files,
    safe_read_tsv
  ),
  use.names = TRUE,
  fill = TRUE
)

if (!nrow(weighted)) {
  stop("Weighted GAM result table is empty.")
}

weighted[
  ,
  region_id := as.integer(
    region_id
  )
]

weighted[
  ,
  data_chunk_id := as.integer(
    data_chunk_id
  )
]

if (anyDuplicated(weighted$region_id)) {
  warning(
    "Duplicated weighted region IDs found. ",
    "Keeping the last occurrence."
  )

  weighted <- weighted[
    order(
      region_id
    )
  ][
    ,
    .SD[
      .N
    ],
    by = region_id
  ]
}

weighted[
  ,
  FDR_s_start_AA :=
    p.adjust(
      p_s_start_AA,
      method = "BH"
    )
]

weighted[
  ,
  DMR_weighted :=
    is.finite(
      FDR_s_start_AA
    ) &
    FDR_s_start_AA < fdr_cutoff &
    is.finite(
      edf_s_start_AA
    ) &
    edf_s_start_AA > edf_cutoff &
    is.finite(
      mean_diff
    ) &
    mean_diff > effect_cutoff
]

fwrite(
  weighted,
  file.path(
    table_dir,
    "weighted_GAM_all_regions.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------- collect original unweighted -----------------------

# Exact column order produced by the original unweighted GAM scripts.
old_names <- c(
  "p_Intercept",
  "p_AA_only",
  "p_AgeCalc",
  "p_Sex",
  "p_Non_smoker",
  "p_EOSINOpc",
  "p_LYMPHOpc",
  "p_MONOpc",
  "p_NEUTROpc",
  "p_sv1",
  "p_sv2",
  "p_sv3",
  "p_sv4",
  "p_sv5",
  "p_BMI",
  "p_s_start",
  "p_s_start_AA",
  "p_s_FID",
  "edf_s_start",
  "edf_s_start_AA",
  "edf_s_FID",
  "R2",
  "AIC",
  "Deviance_explained",
  "REML",
  "N_cpgs",
  "N_samples",
  "data_chunk_id",
  "region_id",
  "max_diff",
  "mean_diff"
)

collect_unweighted_batch <- function(dir_path, batch_name, batch_priority, file_pattern) {

  files <- list.files(
    dir_path,
    pattern = file_pattern,
    full.names = TRUE
  )

  if (!length(files)) {
    warning(
      "No matching original result files found under ",
      dir_path,
      " for batch ",
      batch_name,
      "."
    )
    return(data.table())
  }

  cat(
    "Unweighted ",
    batch_name,
    " result files: ",
    length(files),
    "\n",
    sep = ""
  )

  out_list <- lapply(
    files,
    function(f) {
      z <- safe_read_old(
        f,
        old_names = old_names
      )

      if (is.null(z) || !nrow(z)) {
        return(NULL)
      }

      z[
        ,
        `:=`(
          unweighted_batch = batch_name,
          unweighted_batch_priority = batch_priority,
          unweighted_source_file = basename(f)
        )
      ]

      z
    }
  )

  out_list <- Filter(
    Negate(is.null),
    out_list
  )

  if (!length(out_list)) {
    return(data.table())
  }

  rbindlist(
    out_list,
    use.names = TRUE,
    fill = TRUE
  )
}

unweighted_B1 <- collect_unweighted_batch(
  unweighted_dir_B1,
  batch_name = "BMI_B1",
  batch_priority = 1L,
  file_pattern = "^7_results_job_[0-9]+\\.txt$"
)

unweighted_B2 <- collect_unweighted_batch(
  unweighted_dir_B2,
  batch_name = "BMI_B2",
  batch_priority = 2L,
  file_pattern = "^8_results_job_[0-9]+\\.txt$"
)

if (!nrow(unweighted_B1) && !nrow(unweighted_B2)) {
  stop(
    "No compatible unweighted GAM result files could be read from either ",
    unweighted_dir_B1,
    " or ",
    unweighted_dir_B2,
    "."
  )
}

unweighted_all <- rbindlist(
  list(
    unweighted_B1,
    unweighted_B2
  ),
  use.names = TRUE,
  fill = TRUE
)

unweighted_all[
  ,
  region_id := as.integer(region_id)
]

unweighted_all[
  ,
  data_chunk_id := as.integer(data_chunk_id)
]

# Audit duplicated regions across BMI_B1/BMI_B2.
duplicate_unweighted <- unweighted_all[
  duplicated(region_id) |
    duplicated(region_id, fromLast = TRUE)
][
  order(
    region_id,
    unweighted_batch_priority
  )
]

if (nrow(duplicate_unweighted)) {
  warning(
    "Found ",
    uniqueN(duplicate_unweighted$region_id),
    " duplicated unweighted region IDs across BMI_B1/BMI_B2. ",
    "BMI_B2 is given priority for duplicated regions."
  )

  fwrite(
    duplicate_unweighted,
    file.path(
      table_dir,
      "unweighted_GAM_duplicate_regions_B1_B2.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

# Merge by region. BMI_B2 has priority whenever the same region appears in
# both directories. If a region occurs more than once within the same batch,
# the last occurrence is retained.
setorder(
  unweighted_all,
  region_id,
  unweighted_batch_priority
)

unweighted <- unweighted_all[
  ,
  .SD[.N],
  by = region_id
]

if (anyDuplicated(unweighted$region_id)) {
  stop(
    "Internal error: duplicated region_id values remain after merging BMI_B1/BMI_B2."
  )
}

cat(
  "Unique unweighted regions after merging BMI_B1 + BMI_B2: ",
  nrow(unweighted),
  "\n",
  sep = ""
)

cat(
  "  retained from BMI_B1: ",
  unweighted[
    unweighted_batch == "BMI_B1",
    .N
  ],
  "\n",
  sep = ""
)

cat(
  "  retained from BMI_B2: ",
  unweighted[
    unweighted_batch == "BMI_B2",
    .N
  ],
  "\n",
  sep = ""
)

unweighted[
  ,
  FDR_s_start_AA :=
    p.adjust(
      p_s_start_AA,
      method = "BH"
    )
]

unweighted[
  ,
  DMR_unweighted :=
    is.finite(FDR_s_start_AA) &
    FDR_s_start_AA < fdr_cutoff &
    is.finite(edf_s_start_AA) &
    edf_s_start_AA > edf_cutoff &
    is.finite(mean_diff) &
    mean_diff > effect_cutoff
]

fwrite(
  unweighted,
  file.path(
    table_dir,
    "unweighted_GAM_all_regions_reconstructed.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ---------------------------- matched comparison -----------------------------

keep_weighted <- c(
  "region_id",
  "data_chunk_id",
  "p_s_start_AA",
  "FDR_s_start_AA",
  "edf_s_start_AA",
  "mean_diff",
  "max_diff",
  "R2",
  "N_cpgs",
  "N_samples",
  "N_observations",
  "mean_read_depth",
  "median_read_depth",
  "q25_read_depth",
  "q75_read_depth",
  "DMR_weighted"
)

keep_weighted <- intersect(
  keep_weighted,
  names(weighted)
)

w <- weighted[
  ,
  ..keep_weighted
]

setnames(
  w,
  setdiff(
    names(w),
    c(
      "region_id",
      "data_chunk_id"
    )
  ),
  paste0(
    setdiff(
      names(w),
      c(
        "region_id",
        "data_chunk_id"
      )
    ),
    "_weighted"
  )
)

keep_unweighted <- c(
  "region_id",
  "data_chunk_id",
  "p_s_start_AA",
  "FDR_s_start_AA",
  "edf_s_start_AA",
  "mean_diff",
  "max_diff",
  "R2",
  "N_cpgs",
  "N_samples",
  "DMR_unweighted",
  "unweighted_batch",
  "unweighted_source_file"
)

u <- unweighted[
  ,
  ..keep_unweighted
]

u_model_cols <- setdiff(
  names(u),
  c(
    "region_id",
    "data_chunk_id",
    "unweighted_batch",
    "unweighted_source_file"
  )
)

setnames(
  u,
  u_model_cols,
  paste0(
    u_model_cols,
    "_unweighted"
  )
)

setnames(
  u,
  "unweighted_batch",
  "unweighted_batch_source"
)

comparison <- merge(
  u,
  w,
  by = c(
    "region_id",
    "data_chunk_id"
  ),
  all = TRUE,
  sort = TRUE
)

comparison[
  ,
  matched :=
    is.finite(
      p_s_start_AA_unweighted
    ) &
    is.finite(
      p_s_start_AA_weighted
    )
]

comparison[
  ,
  neglog10p_unweighted :=
    -log10(
      pmax(
        p_s_start_AA_unweighted,
        .Machine$double.xmin
      )
    )
]

comparison[
  ,
  neglog10p_weighted :=
    -log10(
      pmax(
        p_s_start_AA_weighted,
        .Machine$double.xmin
      )
    )
]

comparison[
  ,
  delta_neglog10p :=
    neglog10p_weighted -
    neglog10p_unweighted
]

comparison[
  ,
  delta_mean_diff :=
    mean_diff_weighted -
    mean_diff_unweighted
]

comparison[
  ,
  DMR_status :=
    fifelse(
      DMR_unweighted_unweighted %in% TRUE &
        DMR_weighted_weighted %in% TRUE,
      "Both",
      fifelse(
        DMR_unweighted_unweighted %in% TRUE,
        "Unweighted only",
        fifelse(
          DMR_weighted_weighted %in% TRUE,
          "Weighted only",
          "Neither"
        )
      )
    )
]

fwrite(
  comparison,
  file.path(
    table_dir,
    "weighted_vs_unweighted_region_comparison.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------- concordance ---------------------------------

matched <- comparison[
  matched == TRUE
]

n_unweighted_dmr <- matched[
  DMR_unweighted_unweighted %in% TRUE,
  .N
]

n_weighted_dmr <- matched[
  DMR_weighted_weighted %in% TRUE,
  .N
]

n_both_dmr <- matched[
  DMR_unweighted_unweighted %in% TRUE &
    DMR_weighted_weighted %in% TRUE,
  .N
]

n_union_dmr <- matched[
  DMR_unweighted_unweighted %in% TRUE |
    DMR_weighted_weighted %in% TRUE,
  .N
]

jaccard <- if (
  n_union_dmr > 0
) {
  n_both_dmr /
    n_union_dmr
} else {
  NA_real_
}

retention_unweighted <- if (
  n_unweighted_dmr > 0
) {
  n_both_dmr /
    n_unweighted_dmr
} else {
  NA_real_
}

weighted_supported_by_old <- if (
  n_weighted_dmr > 0
) {
  n_both_dmr /
    n_weighted_dmr
} else {
  NA_real_
}

rho_p <- safe_cor(
  matched$neglog10p_unweighted,
  matched$neglog10p_weighted,
  "spearman"
)

rho_effect <- safe_cor(
  matched$mean_diff_unweighted,
  matched$mean_diff_weighted,
  "spearman"
)

rho_edf <- safe_cor(
  matched$edf_s_start_AA_unweighted,
  matched$edf_s_start_AA_weighted,
  "spearman"
)

pearson_p <- safe_cor(
  matched$neglog10p_unweighted,
  matched$neglog10p_weighted,
  "pearson"
)

pearson_effect <- safe_cor(
  matched$mean_diff_unweighted,
  matched$mean_diff_weighted,
  "pearson"
)

concordance <- data.table(
  metric = c(
    "Matched regions",
    "Spearman rho: -log10 p",
    "Pearson r: -log10 p",
    "Spearman rho: mean_diff",
    "Pearson r: mean_diff",
    "Spearman rho: EDF",
    "Unweighted DMRs",
    "Weighted DMRs",
    "DMRs in both",
    "DMR union",
    "Jaccard overlap",
    "Unweighted DMR retention under weighting",
    "Weighted DMRs also found unweighted"
  ),
  value = c(
    nrow(matched),
    rho_p,
    pearson_p,
    rho_effect,
    pearson_effect,
    rho_edf,
    n_unweighted_dmr,
    n_weighted_dmr,
    n_both_dmr,
    n_union_dmr,
    jaccard,
    retention_unweighted,
    weighted_supported_by_old
  )
)

fwrite(
  concordance,
  file.path(
    table_dir,
    "weighted_GAM_concordance_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ---------------------------------- figures ----------------------------------

p_pvalue <- ggplot(
  matched,
  aes(
    x = neglog10p_unweighted,
    y = neglog10p_weighted
  )
) +
  geom_point(
    alpha = 0.20,
    size = 0.7
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Regional association concordance after read-depth weighting",
    subtitle = paste0(
      "Spearman rho = ",
      fmt(
        rho_p,
        3
      )
    ),
    x = expression(
      -log[10](
        p
      )~
        "unweighted"
    ),
    y = expression(
      -log[10](
        p
      )~
        "read-depth weighted"
    )
  )

ggsave(
  file.path(
    figure_dir,
    "Figure_weighted_vs_unweighted_pvalue_concordance.pdf"
  ),
  p_pvalue,
  width = 7,
  height = 6
)

ggsave(
  file.path(
    figure_dir,
    "Figure_weighted_vs_unweighted_pvalue_concordance.png"
  ),
  p_pvalue,
  width = 7,
  height = 6,
  dpi = 300
)

p_effect <- ggplot(
  matched,
  aes(
    x = mean_diff_unweighted,
    y = mean_diff_weighted
  )
) +
  geom_point(
    alpha = 0.20,
    size = 0.7
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Effect-size concordance after read-depth weighting",
    subtitle = paste0(
      "Spearman rho = ",
      fmt(
        rho_effect,
        3
      )
    ),
    x = "Mean absolute AA effect, unweighted",
    y = "Mean absolute AA effect, read-depth weighted"
  )

ggsave(
  file.path(
    figure_dir,
    "Figure_weighted_vs_unweighted_effect_concordance.pdf"
  ),
  p_effect,
  width = 7,
  height = 6
)

ggsave(
  file.path(
    figure_dir,
    "Figure_weighted_vs_unweighted_effect_concordance.png"
  ),
  p_effect,
  width = 7,
  height = 6,
  dpi = 300
)

dmr_counts <- matched[
  ,
  .N,
  by = DMR_status
]

dmr_counts[
  ,
  DMR_status := factor(
    DMR_status,
    levels = c(
      "Both",
      "Unweighted only",
      "Weighted only",
      "Neither"
    )
  )
]

p_dmr <- ggplot(
  dmr_counts[
    DMR_status != "Neither"
  ],
  aes(
    x = DMR_status,
    y = N
  )
) +
  geom_col() +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "DMR overlap after read-depth weighting",
    x = NULL,
    y = "Number of regions"
  )

ggsave(
  file.path(
    figure_dir,
    "Figure_weighted_vs_unweighted_DMR_overlap.pdf"
  ),
  p_dmr,
  width = 7,
  height = 5.5
)

ggsave(
  file.path(
    figure_dir,
    "Figure_weighted_vs_unweighted_DMR_overlap.png"
  ),
  p_dmr,
  width = 7,
  height = 5.5,
  dpi = 300
)

# ------------------------- manuscript-ready summary --------------------------

summary_lines <- c(
  "READ-DEPTH-WEIGHTED GAM-DMR SENSITIVITY ANALYSIS",
  "================================================",
  "",
  paste0(
    "Matched regions evaluated: ",
    format(
      nrow(matched),
      big.mark = ",",
      scientific = FALSE
    ),
    "."
  ),
  paste0(
    "Unique original unweighted regions after merging BMI_B1 and BMI_B2: ",
    format(
      nrow(unweighted),
      big.mark = ",",
      scientific = FALSE
    ),
    " (BMI_B1 retained=",
    unweighted[unweighted_batch == "BMI_B1", .N],
    "; BMI_B2 retained=",
    unweighted[unweighted_batch == "BMI_B2", .N],
    ")."
  ),
  paste0(
    "DMR definition: BH-FDR < ",
    fdr_cutoff,
    ", EDF > ",
    edf_cutoff,
    ", mean absolute effect > ",
    effect_cutoff,
    "."
  ),
  "",
  "Concordance across all matched regions",
  "--------------------------------------",
  paste0(
    "Spearman correlation of -log10 regional AA p-values: ",
    fmt(
      rho_p,
      3
    ),
    "."
  ),
  paste0(
    "Pearson correlation of -log10 regional AA p-values: ",
    fmt(
      pearson_p,
      3
    ),
    "."
  ),
  paste0(
    "Spearman correlation of mean absolute AA effect sizes: ",
    fmt(
      rho_effect,
      3
    ),
    "."
  ),
  paste0(
    "Pearson correlation of mean absolute AA effect sizes: ",
    fmt(
      pearson_effect,
      3
    ),
    "."
  ),
  paste0(
    "Spearman correlation of AA-specific smooth EDF: ",
    fmt(
      rho_edf,
      3
    ),
    "."
  ),
  "",
  "DMR overlap",
  "-----------",
  paste0(
    "Unweighted DMRs among matched regions: ",
    n_unweighted_dmr,
    "."
  ),
  paste0(
    "Read-depth-weighted DMRs among matched regions: ",
    n_weighted_dmr,
    "."
  ),
  paste0(
    "DMRs detected by both analyses: ",
    n_both_dmr,
    "."
  ),
  paste0(
    "Jaccard overlap: ",
    fmt(
      jaccard,
      3
    ),
    "."
  ),
  paste0(
    "Fraction of original unweighted DMRs retained after weighting: ",
    fmt_pct(
      retention_unweighted,
      1
    ),
    "."
  ),
  paste0(
    "Fraction of weighted DMRs also detected by the unweighted analysis: ",
    fmt_pct(
      weighted_supported_by_old,
      1
    ),
    "."
  ),
  "",
  "Interpretation",
  "--------------",
  paste0(
    "Read-depth weighting was implemented as a sensitivity analysis using ",
    "weights proportional to total read depth and normalized to mean 1 within ",
    "each fitted region. The same ASR transformation, covariates, spline basis, ",
    "basis-dimension rule, family random effect, and REML smoothing framework ",
    "were retained. Concordance of regional statistics and DMR calls therefore ",
    "isolates the effect of accounting for heterogeneous read-depth precision."
  ),
  ""
)

writeLines(
  summary_lines,
  con = file.path(
    comparison_dir,
    "weighted_GAM_sensitivity_results_summary.txt"
  )
)

cat(
  paste(
    summary_lines,
    collapse = "\n"
  ),
  "\n"
)

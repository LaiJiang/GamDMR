#!/usr/bin/env Rscript

# =============================================================================
# 0_prepare_sva_sensitivity_input.R
#
# Prepare the methylation matrix used for the SVA vfilter sensitivity analysis.
#
# Input:
#   ~/scratch/UQAC/meth/dat/peek_methylation_sva_0035.csv
#   (85,936 CpGs retained previously using methylation variance > 0.035)
#
# Output:
#   sva_sensitivity_input.rds
#
# Key choices:
#   * restrict to the same analytic subjects represented in pheno_file
#   * keep methylation proportions only (*_meth); *_tot columns are not needed
#   * ASR transform: asin(sqrt(beta))
#   * row-median impute missing ASR values because sva requires a complete matrix
#   * do NOT preselect 10k/25k/50k here; sva::vfilter will do that later
#
# The model specification itself is constructed in the next script.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

data.table::setDTthreads(1L)
options(mc.cores = 1L)

# ------------------------------- Paths ---------------------------------------

PATH_wk <- path.expand(
  Sys.getenv("METH_BASE_DIR", unset = "~/scratch/UQAC/meth")
)

meth_file <- path.expand(
  Sys.getenv(
    "SVA_METH_FILE",
    unset = file.path(PATH_wk, "dat", "peek_methylation_sva_0035.csv")
  )
)

# Prefer the phenotype-only object to avoid loading large methylation objects.
pheno_file_path <- path.expand(
  Sys.getenv(
    "SVA_PHENO_RDATA",
    unset = file.path(
      PATH_wk, "scr", "11_mgcv", "dat", "bsmooth_only_pheno_bmi.RData"
    )
  )
)

out_dir <- path.expand(
  Sys.getenv(
    "SVA_SENS_OUT",
    unset = file.path(
      PATH_wk, "results", "15_revision", "7_sva_threshold_sensitivity"
    )
  )
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_rds <- file.path(out_dir, "sva_sensitivity_input.rds")
qc_csv  <- file.path(out_dir, "sva_input_qc.csv")

cat("Methylation file:", meth_file, "\n")
cat("Phenotype file:  ", pheno_file_path, "\n")
cat("Output directory:", out_dir, "\n")

if (!file.exists(meth_file)) stop("Missing methylation file: ", meth_file)
if (!file.exists(pheno_file_path)) stop("Missing phenotype file: ", pheno_file_path)

# ------------------------------ Phenotype ------------------------------------

pheno_env <- new.env(parent = emptyenv())
load(pheno_file_path, envir = pheno_env)

if (!exists("pheno_file", envir = pheno_env, inherits = FALSE)) {
  stop("The phenotype RData does not contain an object named pheno_file.")
}

pheno <- as.data.frame(get("pheno_file", envir = pheno_env))
rm(pheno_env)
gc(FALSE)

if (!"ID" %in% names(pheno)) {
  stop("pheno_file must contain an ID column.")
}

# Accept the alternative smoking-variable spellings seen in different scripts.
smoke_candidates <- c("Non-smoker", "Non.smoker", "Non_smoker")
smoke_col <- smoke_candidates[smoke_candidates %in% names(pheno)][1]

if (is.na(smoke_col)) {
  stop(
    "Could not find smoking variable. Expected one of: ",
    paste(smoke_candidates, collapse = ", ")
  )
}

# Standardize its name only for this sensitivity workflow.
pheno$Non_smoker <- pheno[[smoke_col]]

required_covariates <- c(
  "ID",
  "AgeCalc",
  "Sex",
  "Non_smoker",
  "EOSINOpc",
  "LYMPHOpc",
  "MONOpc",
  "NEUTROpc",
  "BMI"
)

missing_covariates <- setdiff(required_covariates, names(pheno))
if (length(missing_covariates) > 0L) {
  stop(
    "Phenotype object is missing required SVA covariates: ",
    paste(missing_covariates, collapse = ", "),
    "\nSet SVA_PHENO_RDATA to an RData file containing the full pheno_file."
  )
}

# Keep exactly subjects with complete covariates for the original SVA model.
pheno_sub <- pheno[
  complete.cases(pheno[, required_covariates]),
  ,
  drop = FALSE
]

pheno_sub$ID <- as.character(pheno_sub$ID)

if (anyDuplicated(pheno_sub$ID)) {
  stop("Duplicate IDs detected in pheno_sub.")
}

cat("Phenotype rows before complete-case filter:", nrow(pheno), "\n")
cat("Phenotype rows after complete-case filter: ", nrow(pheno_sub), "\n")

# The primary manuscript analysis used N=349. Stop rather than silently
# changing the cohort.
EXPECTED_N <- as.integer(
  Sys.getenv("SVA_EXPECTED_N", unset = "349")
)

if (!is.na(EXPECTED_N) && nrow(pheno_sub) != EXPECTED_N) {
  stop(
    "Expected ", EXPECTED_N,
    " analytic subjects, but complete-covariate pheno_sub contains ",
    nrow(pheno_sub),
    ". Check the phenotype input before running the sensitivity analysis."
  )
}

# --------------------------- Methylation header ------------------------------

hdr <- fread(meth_file, nrows = 0L, showProgress = FALSE)
all_cols <- names(hdr)

meth_cols_all <- grep("_meth$", all_cols, value = TRUE)
meth_ids_all <- sub("_meth$", "", meth_cols_all)

# Preserve phenotype order.
wanted_meth_cols <- paste0(pheno_sub$ID, "_meth")
found <- wanted_meth_cols %in% all_cols

if (!all(found)) {
  missing_ids <- pheno_sub$ID[!found]
  stop(
    "Missing methylation columns for ", length(missing_ids),
    " analytic subjects. First missing IDs: ",
    paste(head(missing_ids, 20), collapse = ", ")
  )
}

meta_cols <- intersect(c("chr", "start", "end"), all_cols)
select_cols <- c(meta_cols, wanted_meth_cols)

cat("Reading ", length(select_cols), " selected columns from ", meth_file, "\n", sep = "")

meth_dt <- fread(
  meth_file,
  select = select_cols,
  showProgress = TRUE,
  na.strings = c("NA", "NaN", "")
)

cat("Loaded methylation dimensions:", nrow(meth_dt), "x", ncol(meth_dt), "\n")

if (nrow(meth_dt) < 50000L) {
  stop("Input contains fewer than 50,000 CpGs; cannot evaluate vfilter=50,000.")
}

# ----------------------- Build beta / ASR matrix -----------------------------

meth_dt[, (wanted_meth_cols) := lapply(.SD, as.numeric), .SDcols = wanted_meth_cols]

beta_sub <- as.matrix(meth_dt[, ..wanted_meth_cols])
storage.mode(beta_sub) <- "double"

colnames(beta_sub) <- pheno_sub$ID

# Guard against small floating point excursions outside [0,1].
n_below0 <- sum(beta_sub < 0, na.rm = TRUE)
n_above1 <- sum(beta_sub > 1, na.rm = TRUE)

if (n_below0 > 0L || n_above1 > 0L) {
  warning(
    "Found values outside [0,1]: below 0 = ", n_below0,
    "; above 1 = ", n_above1,
    ". Values will be clipped before ASR transformation."
  )
  beta_sub <- pmin(pmax(beta_sub, 0), 1)
}

beta_asr <- asin(sqrt(beta_sub))

# -------------------------- Missing-value QC --------------------------------

n_obs_before <- rowSums(is.finite(beta_asr))
n_missing_before <- sum(!is.finite(beta_asr))

cat("Total missing ASR entries before imputation:", n_missing_before, "\n")
cat(
  "Observed subjects per CpG: median=", median(n_obs_before),
  "; IQR=", paste(quantile(n_obs_before, c(.25, .75)), collapse = "-"),
  "; range=", paste(range(n_obs_before), collapse = "-"),
  "\n",
  sep = ""
)

# SVA does not accept missing values.
# Row-median imputation preserves each CpG's central methylation level and
# avoids introducing between-subject structure.
row_median_impute <- function(mat) {

  if (requireNamespace("matrixStats", quietly = TRUE)) {

    med <- matrixStats::rowMedians(mat, na.rm = TRUE)

  } else {

    med <- apply(mat, 1L, stats::median, na.rm = TRUE)

  }

  # CpGs with no finite value cannot be imputed and are removed.
  keep <- is.finite(med)

  if (!all(keep)) {
    cat("Removing ", sum(!keep), " CpGs with no finite analytic-subject values.\n", sep = "")
    mat <- mat[keep, , drop = FALSE]
    med <- med[keep]
  }

  ij <- which(!is.finite(mat), arr.ind = TRUE)

  if (nrow(ij) > 0L) {
    mat[ij] <- med[ij[, 1L]]
  }

  list(matrix = mat, keep = keep, medians = med)
}

imp <- row_median_impute(beta_asr)
beta_asr_imputed <- imp$matrix
keep_rows <- imp$keep

if (any(!is.finite(beta_asr_imputed))) {
  stop("Non-finite values remain after imputation.")
}

# Corresponding CpG metadata.
cpg_info <- meth_dt[, ..meta_cols]
if (!all(keep_rows)) cpg_info <- cpg_info[keep_rows]

# Remove rows with zero variance after imputation, because they cannot
# contribute to variable-CpG selection.
if (requireNamespace("matrixStats", quietly = TRUE)) {
  row_var_asr <- matrixStats::rowVars(beta_asr_imputed)
} else {
  row_var_asr <- apply(beta_asr_imputed, 1L, stats::var)
}

keep_var <- is.finite(row_var_asr) & row_var_asr > 0

cat("CpGs before zero-variance filter:", nrow(beta_asr_imputed), "\n")
cat("CpGs after zero-variance filter: ", sum(keep_var), "\n")

beta_asr_imputed <- beta_asr_imputed[keep_var, , drop = FALSE]
cpg_info <- cpg_info[keep_var]
row_var_asr <- row_var_asr[keep_var]

if (nrow(beta_asr_imputed) < 50000L) {
  stop(
    "Fewer than 50,000 non-zero-variance CpGs remain after analytic-subject ",
    "subsetting/imputation: ", nrow(beta_asr_imputed)
  )
}

# Save a stable CpG identifier.
if (all(c("chr", "start") %in% names(cpg_info))) {
  cpg_info[, cpg_id := paste0(chr, ":", start)]
} else {
  cpg_info[, cpg_id := paste0("CpG_", seq_len(.N))]
}

cpg_info[, var_asr_analytic := row_var_asr]

# ------------------------------- Save ----------------------------------------

obj <- list(
  beta_sub = beta_asr_imputed,  # rows CpGs x columns subjects
  pheno_sub = pheno_sub,
  cpg_info = cpg_info,
  input_file = meth_file,
  phenotype_file = pheno_file_path,
  transformation = "asin(sqrt(beta))",
  imputation = "CpG-wise median on ASR scale",
  n_cpg_input = nrow(meth_dt),
  n_cpg_final = nrow(beta_asr_imputed),
  n_subjects = ncol(beta_asr_imputed)
)

saveRDS(obj, out_rds, compress = FALSE)

qc <- data.table(
  input_cpgs = nrow(meth_dt),
  final_cpgs = nrow(beta_asr_imputed),
  subjects = ncol(beta_asr_imputed),
  missing_values_before_imputation = n_missing_before,
  median_observed_subjects_per_cpg = median(n_obs_before),
  q1_observed_subjects_per_cpg = unname(quantile(n_obs_before, .25)),
  q3_observed_subjects_per_cpg = unname(quantile(n_obs_before, .75)),
  min_observed_subjects_per_cpg = min(n_obs_before),
  max_observed_subjects_per_cpg = max(n_obs_before)
)

fwrite(qc, qc_csv)

cat("\nPrepared SVA sensitivity input:\n")
cat("  ", out_rds, "\n", sep = "")
cat("QC summary:\n")
print(qc)

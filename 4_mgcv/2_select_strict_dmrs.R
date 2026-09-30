suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

PATH_wk <- path.expand(Sys.getenv("METH_BASE_DIR", "~/scratch/UQAC/meth"))
input_dir <- file.path(PATH_wk, "results", "4_mgcv", "B1")
output_dir <- file.path(PATH_wk, "results", "11_mgcv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(input_dir, pattern = "^results_job_[0-9]+\\.txt$", full.names = TRUE)
if (!length(files)) stop("No GAM result files found in: ", input_dir)

fixed <- c(
  "(Intercept)", "AA_only", "AgeCalc", "Sex", "Non.smoker",
  "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc",
  "sv1", "sv2", "sv3", "sv4", "sv5", "BMI"
)
smooth <- c("s(start)", "s(start):AA_only", "s(FID)")
edf <- paste0("edf_", smooth)
cols <- c(
  fixed, smooth, edf,
  "R2", "AIC", "Deviance_explained", "REML",
  "N_cpgs", "N_samples", "data_chunk_id", "region_id", "max_diff", "mean_diff"
)

all_results <- rbindlist(lapply(files, function(f) {
  x <- fread(f, header = TRUE, sep = "\t")
  if (ncol(x) != length(cols)) stop("Unexpected column count: ", f)
  setnames(x, cols)
  x[, job_id := as.integer(str_extract(basename(f), "(?<=_job_)\\d+"))]
  x
}), use.names = TRUE, fill = TRUE)

if (anyDuplicated(all_results$region_id)) {
  dup <- unique(all_results$region_id[duplicated(all_results$region_id)])
  stop("Duplicate region_id values detected: ", paste(head(dup, 20L), collapse = ", "))
}

all_results[, pvals := as.numeric(`s(start):AA_only`)]
all_results[, edf := as.numeric(`edf_s(start):AA_only`)]
all_results[, pval_adjusted := NA_real_]
valid <- which(is.finite(all_results$pvals))
all_results[valid, pval_adjusted := p.adjust(pvals, method = "BH")]

strict <- all_results[
  is.finite(pvals) & pvals < 1e-5 &
    is.finite(pval_adjusted) & pval_adjusted < 0.05 &
    is.finite(edf) & edf >= 0.5 &
    is.finite(R2) & R2 >= 0.50 &
    is.finite(mean_diff) & mean_diff >= 0.05 &
    N_cpgs >= 10
]
setorder(strict, pvals)

fwrite(all_results, file.path(output_dir, "23_mgcv_all_regions.tsv"), sep = "\t")
fwrite(strict, file.path(output_dir, "24_dmrs_STRICT.tsv"), sep = "\t")

cat("Tested regions:", nrow(all_results), "\n")
cat("Strict GAM-DMRs:", nrow(strict), "\n")

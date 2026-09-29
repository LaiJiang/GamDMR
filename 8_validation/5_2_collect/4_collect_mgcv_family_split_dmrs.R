
suppressPackageStartupMessages({
  library(data.table)
})

PATH_wk <- path.expand(Sys.getenv("METH_BASE_DIR", "~/scratch/UQAC/meth/"))
raw_dir <- file.path(PATH_wk, "results", "15_revision", "5_cv", "1_mgcv")
region_path <- file.path(PATH_wk, "scr", "11_mgcv", "dat", "region_file_1_chunk.csv")
out_dir <- file.path(raw_dir, "collected_by_split", "strict_by_split")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(raw_dir, pattern = "^2_run_mgcv_training_splits_results_[0-9]+\\.txt$", full.names = TRUE)
if (!length(files)) stop("No MGCV validation result files found.")

fixed <- c(
  "(Intercept)", "AA_only", "AgeCalc", "Sex", "Non.smoker",
  "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc",
  "sv1", "sv2", "sv3", "sv4", "sv5", "BMI"
)
smooth <- c("s(start)", "s(start):AA_only", "s(FID)")
edf <- paste0("edf_", smooth)
cols <- c(
  "splitID", fixed, smooth, edf,
  "R2", "AIC", "Deviance_explained", "REML",
  "N_cpgs", "N_samples", "N_train_FIDs",
  "data_chunk_id", "region_id", "max_diff", "mean_diff"
)

all <- rbindlist(lapply(files, function(f) {
  x <- fread(f, header = FALSE)
  if (ncol(x) != length(cols)) stop("Unexpected column count: ", f)
  setnames(x, cols)
  x
}), use.names = TRUE, fill = TRUE)

regions <- fread(region_path)
if (!"region_id" %in% names(regions)) regions[, region_id := .I]
coord_cols <- intersect(c("region_id", "chr", "region_start", "region_end"), names(regions))
all <- merge(all, regions[, ..coord_cols], by = "region_id", all.x = TRUE)

split_summary <- list()
all_strict <- list()

for (sid in sort(unique(all$splitID))) {
  x <- copy(all[splitID == sid])
  x[, FDR := p.adjust(`s(start):AA_only`, method = "BH")]
  strict <- x[
    `s(start):AA_only` < 1e-5 &
      FDR < 0.05 &
      `edf_s(start):AA_only` >= 0.5 &
      R2 >= 0.50 &
      mean_diff >= 0.05 &
      N_cpgs >= 10
  ]
  file <- file.path(out_dir, sprintf("24_dmrs_STRICT_split_%03d.tsv", sid))
  fwrite(strict, file, sep = "\t")
  all_strict[[as.character(sid)]] <- strict
  split_summary[[as.character(sid)]] <- data.table(
    splitID = sid,
    N_tested = nrow(x),
    N_DMR = nrow(strict)
  )
}

fwrite(
  rbindlist(all_strict, fill = TRUE),
  file.path(raw_dir, "24_dmrs_STRICT_all_splits.tsv"),
  sep = "\t"
)
fwrite(
  rbindlist(split_summary),
  file.path(raw_dir, "GAM_DMR_split_summary.tsv"),
  sep = "\t"
)

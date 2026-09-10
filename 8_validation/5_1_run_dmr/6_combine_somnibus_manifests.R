#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

PATH_wk <- path.expand("~/scratch/UQAC/meth/")
output_root <- file.path(
  PATH_wk, "results/15_revision/5_cv/2_somnibus/prep"
)
manifest_dir <- file.path(output_root, "manifests")
N_jobs <- as.integer(Sys.getenv("N_JOBS", "200"))

files <- list.files(
  manifest_dir,
  pattern = "^somnibus_manifest_job_[0-9]{4}\\.tsv$",
  full.names = TRUE
)

if (length(files) == 0L) stop("No manifest files found.")

job_ids <- as.integer(
  sub(
    "^somnibus_manifest_job_([0-9]{4})\\.tsv$",
    "\\1",
    basename(files)
  )
)

missing_jobs <- setdiff(seq_len(N_jobs), job_ids)
if (length(missing_jobs) > 0L) {
  warning(
    "Missing ", length(missing_jobs), " job manifests. First IDs: ",
    paste(head(missing_jobs, 20L), collapse = ", ")
  )
}

manifest <- rbindlist(
  lapply(files, fread),
  use.names = TRUE,
  fill = TRUE
)

setorder(
  manifest, data_chunk_id, parent_region_id, child_index,
  na.last = TRUE
)

ready <- manifest[status == "READY"]

if (anyDuplicated(ready$omnibus_region_id)) {
  stop("Duplicated READY omnibus_region_id values found.")
}

missing_bundles <- unique(
  ready[!file.exists(bundle_file), bundle_file]
)
if (length(missing_bundles) > 0L) {
  stop(
    "READY rows point to missing bundle files: ",
    paste(head(missing_bundles, 10L), collapse = ", ")
  )
}

combined_file <- file.path(output_root, "somnibus_manifest_all.tsv")
summary_file <- file.path(output_root, "somnibus_manifest_summary.tsv")

fwrite(
  manifest, combined_file,
  sep = "\t", quote = FALSE, na = "NA"
)

summary_dt <- manifest[
  ,
  .(
    manifest_rows = .N,
    parent_regions = uniqueN(parent_region_id),
    ready_regions = sum(status == "READY"),
    ready_cpgs = sum(
      fifelse(status == "READY", n_cpgs, 0L),
      na.rm = TRUE
    )
  ),
  by = status
][order(status)]

fwrite(
  summary_dt, summary_file,
  sep = "\t", quote = FALSE, na = "NA"
)

cat("Manifest files found:", length(files), "\n")
cat("Missing job manifests:", length(missing_jobs), "\n")
cat("Total manifest rows:", nrow(manifest), "\n")
cat("READY SOMNiBUS regions:", nrow(ready), "\n")
cat("Combined manifest:", combined_file, "\n")
cat("Summary:", summary_file, "\n")
print(summary_dt)

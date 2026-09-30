
suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- path.expand(Sys.getenv("METH_BASE_DIR", "~/scratch/UQAC/meth"))
paper_dir <- file.path(base_dir, "scr", "14_paper")
input_dir <- file.path(base_dir, "results", "13_dmrcate")
output_file <- file.path(paper_dir, "results", "7_DMRcate_DMRs.csv")
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

files <- list.files(input_dir, pattern = "^dmrcate_job_.*\\.tsv$", full.names = TRUE)
if (length(files) == 0L) stop("No DMRcate DMR files found in: ", input_dir)

x <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)
if (!all(c("region_start", "region_end") %in% names(x))) {
  if (all(c("start", "end") %in% names(x))) {
    x[, `:=`(region_start = as.integer(start), region_end = as.integer(end))]
  } else {
    stop("DMRcate DMR coordinates are missing.")
  }
}
if (!"chr" %in% names(x)) {
  if ("seqnames" %in% names(x)) x[, chr := as.character(seqnames)] else stop("DMRcate chromosome is missing.")
}
x[, chr := as.character(chr)]
x[, width := as.integer(region_end) - as.integer(region_start) + 1L]
x[, call_id := .I]

fwrite(x, output_file)

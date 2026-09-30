
suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- path.expand(Sys.getenv("METH_BASE_DIR", "~/scratch/UQAC/meth"))
paper_dir <- file.path(base_dir, "scr", "14_paper")
input_dir <- file.path(base_dir, "results", "13_bsmooth")
output_file <- file.path(paper_dir, "results", "6_BSmooth_DMRs.csv")
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

files <- list.files(input_dir, pattern = "^dmr_job_.*\\.tsv(\\.gz)?$", full.names = TRUE)
if (length(files) == 0L) stop("No BSmooth DMR files found in: ", input_dir)

x <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)
if (!all(c("chr", "start", "end") %in% names(x))) stop("BSmooth DMR coordinates are missing.")

x[, `:=`(
  chr = as.character(chr),
  region_start = as.integer(start),
  region_end = as.integer(end)
)]
x[, width := region_end - region_start + 1L]
x[, call_id := .I]

fwrite(x, output_file)

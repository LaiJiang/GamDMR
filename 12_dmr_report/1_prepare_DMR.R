
suppressPackageStartupMessages({
  library(data.table)
})

PATH_wk <- path.expand(Sys.getenv("METH_BASE_DIR", "~/scratch/UQAC/meth/"))
paper_results <- file.path(PATH_wk, "scr", "14_paper", "results")
output_dir <- file.path(PATH_wk, "results", "15_revision", "9_dmr_report", "DMR_input")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gam <- fread(file.path(paper_results, "9_table_2.csv"))
somnibus <- fread(file.path(PATH_wk, "results", "9_regional", "16_somnibus_regions.csv"))
bsmooth <- fread(file.path(paper_results, "6_BSmooth_DMRs.csv"))
dmrcate <- fread(file.path(paper_results, "7_DMRcate_DMRs.csv"))

coord_table <- function(x, method) {
  if (!"chr" %in% names(x)) {
    if ("seqnames" %in% names(x)) x[, chr := as.character(seqnames)]
  }
  if (!all(c("region_start", "region_end") %in% names(x))) {
    if (all(c("start", "end") %in% names(x))) {
      x[, `:=`(region_start = as.integer(start), region_end = as.integer(end))]
    }
  }
  required <- c("chr", "region_start", "region_end")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(method, " is missing: ", paste(missing, collapse = ", "))
  x[, .(
    method = method,
    DMR_chr = as.character(chr),
    DMR_region_start = as.integer(region_start),
    DMR_region_end = as.integer(region_end)
  )]
}

all_DMRs_long <- rbindlist(
  list(
    coord_table(gam, "GAM-DMR"),
    coord_table(somnibus, "SOMNiBUS"),
    coord_table(bsmooth, "BSmooth"),
    coord_table(dmrcate, "DMRcate")
  ),
  use.names = TRUE
)
all_DMRs_long[, DMR_width_bp := DMR_region_end - DMR_region_start + 1L]
fwrite(all_DMRs_long, file.path(output_dir, "all_DMRs_long.csv"))
print(all_DMRs_long[, .N, by = method])

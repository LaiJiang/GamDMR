
suppressPackageStartupMessages({
  library(data.table)
})

PATH_wk <- path.expand(Sys.getenv("METH_BASE_DIR", "~/scratch/UQAC/meth/"))
PATH_scr11 <- file.path(PATH_wk, "scr", "11_mgcv")
PATH_save <- file.path(PATH_wk, "results", "15_revision", "5_cv")
dir.create(PATH_save, recursive = TRUE, showWarnings = FALSE)

regions <- fread(file.path(PATH_scr11, "dat", "region_file_1_chunk.csv"))
if (!"region_id" %in% names(regions)) regions[, region_id := .I]
if ("n_cpgs" %in% names(regions)) regions <- regions[n_cpgs >= 10L]

fwrite(regions, file.path(PATH_save, "region_DMRcate_M12.csv"))
fwrite(regions, file.path(PATH_scr11, "dat", "region_DMRcate_M12.csv"))

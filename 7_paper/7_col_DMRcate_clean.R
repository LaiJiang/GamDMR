suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

base_dir <- Sys.getenv("METH_BASE_DIR", unset = getwd())
paper_dir <- file.path(base_dir, "14_paper")
input_dir <- file.path(base_dir, "results", "13_dmrcate")
output_file <- file.path(paper_dir, "results", "7_DMRcate_DMRs.csv")

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

summary_files <- list.files(
  path = input_dir,
  pattern = "^summary_job_.*\\.tsv$",
  full.names = TRUE
)

if (length(summary_files) == 0) {
  stop("No DMRcate summary files found in: ", input_dir)
}

merged_df <- rbindlist(lapply(summary_files, fread), use.names = TRUE, fill = TRUE)
dmrcate_dmrs <- merged_df %>% filter(n_dmrs >= 1)

fwrite(dmrcate_dmrs, file = output_file)

#collect the BSmooth results for the paper

    library(data.table)
    library(dplyr)
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_scr14 <- paste0(PATH_wk, "/scr/14_paper/")


PATH_dat <- paste0(PATH_wk, "results/13_dmrcate/stringent/")


# ----------------------------------------
# Directory containing all summary_job_*.tsv files
# ----------------------------------------
dir_path <- PATH_dat

# ----------------------------------------
# List all TSV files
# ----------------------------------------
files <- list.files(
  path = dir_path,
  pattern = "^summary_job_.*\\.tsv$",
  full.names = TRUE
)

# Print number of files found
cat("Found", length(files), "files.\n")

# ----------------------------------------
# Load all TSV files and merge
# ----------------------------------------
library(data.table)

merged_df <- rbindlist(
  lapply(files, fread),
  use.names = TRUE,
  fill = TRUE
)

# ----------------------------------------
# Output
# ----------------------------------------
cat("Merged rows:", nrow(merged_df), "\n")
cat("Merged columns:", ncol(merged_df), "\n")

# Inspect first few rows
head(merged_df)


table(merged_df$n_dmrs)

DMRcate_DMRs <- merged_df %>%
  filter(n_dmrs >= 1)


write.csv(DMRcate_DMRs, file = paste0(PATH_scr14, "results/7_DMRcate_DMRs.csv"), row.names = FALSE)
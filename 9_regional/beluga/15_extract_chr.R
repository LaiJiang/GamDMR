# Libraries
library(data.table)
library(dplyr)
library(stringr)
library(SOMNiBUS)

PATH_wk <- "~/scratch/UQAC/meth/"
PATH_scr6 <- paste0(PATH_wk,"scr/6_whole/")
PATH_scr8 <- paste0(PATH_wk,"scr/8_rerun/")
PATH_scr9 <- paste0(PATH_wk, "scr/9_regional/")

# Output file
results_file <- paste0(PATH_wk, "results/9_regional/15_extract_chr.csv")

# Load regions
somnibus_regions <- read.csv(paste0(PATH_wk, "results/9_regional/14_eval_pvals.csv"))

# If the file already exists, remove it to start fresh
if (file.exists(results_file)) {
  file.remove(results_file)
}

# Loop over regions
for(i_row in 1:nrow(somnibus_regions)) {

  i_chunk_id <- somnibus_regions$chunk_id[i_row]
  meth_file_loc <- paste0(PATH_wk, "data/meth_split/chunk_", sprintf("%04d", i_chunk_id), ".csv")
  
  # Read the methylation file
  meth_file <- read.csv(file = meth_file_loc, header = TRUE, sep = "\t")
  
  # Extract region boundaries
  region_start <- somnibus_regions[i_row, "region_start"]
  region_end   <- somnibus_regions[i_row, "region_end"]
  
  # Filter methylation CpGs overlapping with the region
  meth_file_overlap <- meth_file %>%
    filter(start >= region_start & start <= region_end)
  
  if (nrow(meth_file_overlap) > 0) {
    chr_string <- paste(meth_file_overlap$chr, collapse = ":")
    pos_string <- paste(meth_file_overlap$start, collapse = ":")
    
    output_data <- data.frame(
      region_id = somnibus_regions[i_row, "region_id"],
      chunk_id = i_chunk_id,
      chr = chr_string,
      start = pos_string
    )
    
    # Append to the results file
    write.table(output_data, file = results_file, sep = ",", 
                row.names = FALSE, col.names = !file.exists(results_file), 
                append = TRUE, quote = FALSE)
  }
}

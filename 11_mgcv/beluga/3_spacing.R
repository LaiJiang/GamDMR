#submit the cpg spacing startegy to whole genome onto beluga.

#SBATCH --array=1-10

###########################################################################
PATH_wk <- "~/scratch/UQAC/meth/"

PATH_scr11 <- paste0(PATH_wk, "scr/11_mgcv/")


# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No command-line arguments supplied. Please pass the SLURM_ARRAY_TASK_ID.")
}

# Convert the first argument to an integer
i_job_id <- as.integer(args[1])
cat("Job ID is:", i_job_id, "\n")

#split array 1:1737 into 10 arrays
vec_list <- split(1:1737, cut(seq_along(1:1737), breaks = 10, labels = FALSE))

i_chunk_id_vector <- vec_list[[i_job_id]]


library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(tidyr)
library(mgcv)

for(i_chunk_id in i_chunk_id_vector) {
  #now load that data chunk
  #i_chunk <- fread(paste0(PATH_wk, "results/10_sanity/chunk_", sprintf("%04d", i_chunk), ".csv"))
  i_chunk <- fread(paste0(PATH_wk, "data/meth_split/chunk_", sprintf("%04d", i_chunk_id), ".csv"))
  
  i_cpg_info <- i_chunk %>%
    select(chr, start) %>%
    distinct() %>%
    rename(position = start)
  
  #######################
  #######################
  #load the function of split data by cpg spacing
  source(paste0(PATH_scr11, "0_split_spacing.R"))
  
  cpg_regions <- splitCpGsBySpacing(i_cpg_info, gap = 250, min_cpgs = 10, max_cpgs = 100) 
  
  cat("Processing chunk: ", i_chunk_id, "\n")

  #save the file, remove column names
  fwrite(cpg_regions, paste0(PATH_wk, "data/region_spacing/spacing_", sprintf("%04d", i_chunk_id), ".csv"), col.names = FALSE)
}

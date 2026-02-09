library(globaltest)

#test run somnibus on a single region
#myabe also bump hunter? https://bioconductor.org/packages/release/bioc/html/bumphunter.html

#following such results: https://www.bioconductor.org/packages/release/bioc/vignettes/SOMNiBUS/inst/doc/SOMNiBUS.html


library(dplyr)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No command-line arguments supplied. Please pass the SLURM_ARRAY_TASK_ID.")
}


if(FALSE){
# Convert the first argument to an integer
task_id <- as.integer(args[1])
cat("SLURM_ARRAY_TASK_ID is:", task_id, "\n")
}

#first data chunk
task_id <- 690
###########################################################################
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_scr6 <- paste0(PATH_wk,"scr/6_beluga/")
PATH_scr8 <- paste0(PATH_wk,"scr/8_rerun/")
PATH_scr9 <- paste0(PATH_wk, "scr/9_regional/")

#first to load one splited data C:\Per\LaiJiang\Project\UQAC\meth\dat\chunk_017.csv
meth_file_loc <- paste0(PATH_wk,"scr/6_beluga/chunk_",sprintf("%04d",task_id),".csv")

#location to save reults
results_loc <-   paste0(PATH_scr9,"results/AA_",sprintf("%04d",task_id))

###########################################################################
meth_file <- read.csv(file=meth_file_loc, header=TRUE, sep="\t")


#load the pheno_file from the created data
load(file=paste0(PATH_scr6,"4_1_data.RData"), verbose=TRUE)
###########################################################################
#

#first split data into regions with 1Mb distance
# Sort by chromosome and start (just in case)
meth_file <- meth_file[order(meth_file$chr, meth_file$start), ]

# Compute distance between consecutive CpGs
gap <- c(0, diff(meth_file$start))

# Identify new regions where gap > 1Mb (1e6)
new_region <- gap > 1e6

# Assign region IDs using cumulative sum of new region indicators
meth_file$region_id <- cumsum(new_region) + 1  # +1 to start region IDs from 1

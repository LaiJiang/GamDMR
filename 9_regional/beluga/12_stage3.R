#this file select the regions of interst 
library(data.table)
library(dplyr)
library(stringr)
library(SOMNiBUS)

PATH_wk <- "~/scratch/UQAC/meth/"


PATH_scr6 <- paste0(PATH_wk,"scr/6_whole/")

PATH_scr8 <- paste0(PATH_wk,"scr/8_rerun/")
PATH_scr9 <- paste0(PATH_wk, "scr/9_regional/")


#sel_covariate <- c("AgeCalc", "meth_AA","Non.smoker", "EOSINOpc", "LYMPHOpc")
#sel_covariate <- c("AgeCalc", "AA_only", "Non.smoker", "EOSINOpc", "LYMPHOpc")
sel_covariate <- c( "AA_only")

######################################################################################
######################################################################################
######################################################################################
#now load a data chunk


# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No command-line arguments supplied. Please pass the SLURM_ARRAY_job_id.")
}


if(TRUE){
# Convert the first argument to an integer
job_id <- as.integer(args[1])
cat("SLURM_ARRAY_job_id is:", job_id, "\n")
}

#first data chunk
#job_id <- 690
###########################################################################


#first to load one splited data C:\Per\LaiJiang\Project\UQAC\meth\dat\chunk_017.csv
region_file_loc <- paste0(PATH_wk,"results/somnibus_stage3/data/region_",sprintf("%d",job_id),".RData")

#location to save reults
results_loc <-   paste0(PATH_wk,"results/somnibus_stage3/results/AA_",sprintf("%d",job_id), ".RData")

###########################################################################

load(region_file_loc,verbose=TRUE)

######################################################################
#if there are cpgs worth looking at, then proceed to the next step



#step: run Somnibus analysis
#we use the default region split method:  which splits methylation data into regions based on the spacing of CpGs.
#estimate n.k by the default method
n_k_dim <- max(3, as.integer(length(unique(region_data$Position)) / 20))



outs <- binomRegMethModel(data = region_data, n.k = rep(n_k_dim,2), p0 = 0.003, 
                         p1 = 0.9, Quasi = FALSE, RanEff = FALSE,
                         verbose = TRUE)


###############################################################################

###########################################

###########################################

#note the task_id is the original data chunk id
data_chunk_id <- task_id

save( data_chunk_id, outs, file=  results_loc   )

    




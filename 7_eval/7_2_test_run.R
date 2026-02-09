#this file test run the rerun of the cpgs that need to get pvalue from M1 and M3


library(lmerTest)
library(glmmLasso)
library(lme4)

if(FALSE){
# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No command-line arguments supplied. Please pass the SLURM_ARRAY_TASK_ID.")
}

# Convert the first argument to an integer
task_id <- as.integer(args[1])
cat("SLURM_ARRAY_TASK_ID is:", task_id, "\n")
}

#first data chunk
task_id <- 1
###########################################################################
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_scr6 <- paste0(PATH_wk,"scr/6_whole/")

#first to load one splited data C:\Per\LaiJiang\Project\UQAC\meth\dat\chunk_017.csv
meth_file_loc <- paste0(PATH_wk,"data/meth_split/chunk_",sprintf("%04d",task_id),".csv")

#location to save reults
results_loc <-   paste0(PATH_wk,"results/whole_genome/AA_",sprintf("%04d",task_id))

###########################################################################
meth_file <- read.csv(file=meth_file_loc, header=TRUE, sep="\t")

#methfile contains the JobID!!!
#the saving lcoation also contains the JobID!!!

#load the pheno_file from the created data
#load(file=paste0(PATH_scr6,"4_1_data.RData"), verbose=TRUE)
load(file=paste0(PATH_wk,"results/4_1_data.RData"), verbose=TRUE  )


#1. match the IDs of AA and control 
# Extract only the methylation proportion columns (those that end with '_meth')
meth_data <- meth_file[, grep("_meth$", colnames(meth_file))]

# Convert the data frame to a matrix
meth_matrix <- as.matrix(meth_data)

#colect the IDs for meth_matrix
meth_ID <- colnames(meth_matrix)
meth_ID <- sub("^X", "", meth_ID)      # Remove the leading "X"
meth_ID <- sub("_meth$", "", meth_ID)   # Remove the trailing "_meth"
meth_ID <- gsub("\\.", "-", meth_ID)     # Replace the dot with a dash

# Check dimensions and a preview of the matrix

meth_matrix <- meth_matrix[,match(pheno_file$ID,meth_ID)]



####################################################################
#then create methylation data for a single cpg

i_cpg <- 1

meth_file$location <- paste(meth_file$chr, paste(meth_file$start, meth_file$end, sep = "-"), sep = ":")


for(i_cpg in 1:nrow(meth_matrix)){
###########################################

###########################################
print(paste0("ignore NA: CpG ",i_cpg))

#the asr transformation
meth_response <- asin(sqrt(meth_matrix[i_cpg,]))

#record cpg location
meth_cpg_info <- meth_file$location[i_cpg]

pheno_file$meth_response <- meth_response

#remove NA rows

pheno_file_feed <- na.omit(pheno_file)
#how many NA in each column


#if there are enough sample sizes: n > 50
if(nrow(pheno_file_feed) > 30){

  if(sd(pheno_file_feed$meth_response)!=0){

source(paste0(PATH_scr6,"/6_0_func_single.R"))

#attach this line to a file
write.table(t(output_complete_results), file=paste0(results_loc,"_complete.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(t(output_simple_results), file=paste0(results_loc,"_simple.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")

        

  }

}


}




library(lmerTest)
library(glmmLasso)
library(lme4)
library(data.table)


  
# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No command-line arguments supplied. Please pass the SLURM_ARRAY_TASK_ID.")
}



# Convert the first argument to an integer
task_id <- as.integer(args[1])
cat("SLURM_ARRAY_TASK_ID is:", task_id, "\n")


###########################################################################
PATH_wk <- "~/scratch/UQAC/meth/"


PATH_scr6 <- paste0(PATH_wk,"scr/6_whole/")
PATH_scr8 <- paste0(PATH_wk,"scr/8_rerun/")
PATH_scr11 <- paste0(PATH_wk,"scr/11_mgcv/")
PATH_scr12 <- paste0(PATH_wk,"scr/12_cpg_DMR/")

#first to load one splited data C:\Per\LaiJiang\Project\UQAC\meth\dat\chunk_017.csv
meth_file_loc <- paste0(PATH_wk,"data/meth_split/chunk_",sprintf("%04d",task_id),".csv")

#location to save reults!!!
results_loc <-   paste0(PATH_wk,"results/12_cpg_DMR/AA_",sprintf("%04d",task_id))

###########################################################################
meth_file <- read.csv(file=meth_file_loc, header=TRUE, sep="\t")

#methfile contains the JobID!!!
#the saving lcoation also contains the JobID!!!

#load the pheno_file from the created data
#load(file=paste0(PATH_scr6,"4_1_data.RData"), verbose=TRUE)
load(file = paste0(PATH_scr11, "dat/18_pheno_BMI.RData"), verbose = TRUE)


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



meth_file$location <- paste(meth_file$chr, paste(meth_file$start, meth_file$end, sep = "-"), sep = ":")

################################################################################################################################
################################################################################################################################################################
################################################################################################################################
#this is what we need to add to the rerun file in addtion to the original run files!!!!!!!!!!!!!

#first load the rerun id file 

rerun_ids <- fread(paste0(PATH_scr8,"/dat/2_col_rerun.txt"), header=FALSE)

rerun_chunk <- rerun_ids[rerun_ids$V1 == task_id,]

rerun_i_cpg <-rerun_chunk$V2


#i_cpg <- 12
#################################!!!!!!!!!!!!!
for(i_cpg in rerun_i_cpg){
###########################################

###########################################
print(paste0("process: CpG ",i_cpg))

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

#!!!!!!!!!!!!!!different
source(paste0(PATH_scr12,"/0_func_single_rerun.R"))

#attach this line to a file
write.table(t(output_complete_results), file=paste0(results_loc,"_complete.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(t(output_simple_results), file=paste0(results_loc,"_simple.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")

        

  }

}


}



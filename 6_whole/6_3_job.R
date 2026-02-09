
library(lmerTest)
library(glmmLasso)
library(lme4)



#first to load one splited data C:\Per\LaiJiang\Project\UQAC\meth\dat\chunk_017.csv
meth_file <- read.csv(file="C:/Per/LaiJiang/Project/UQAC/meth/dat/chunk_017.csv", header=TRUE, sep="\t")


#load the pheno_file from the created data
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
load(file=paste0(PATH_wk,"results/4_1_data.RData"), verbose=TRUE)


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


for(i_cpg in 1:100){
###########################################

###########################################
print(paste0("ignore NA: CpG ",i_cpg))

#the asr transformation
meth_response <- asin(sqrt(meth_matrix[i_cpg,]))

pheno_file$meth_response <- meth_response

#remove NA rows

pheno_file_feed <- na.omit(pheno_file)
#how many NA in each column


#if there are enough sample sizes: n > p, p = 13 fixed effects and 1 FID
if(nrow(pheno_file_feed) > 14){
source(paste0(PATH_wk,"scr/6_whole/6_0_func_single.R"))

#attach this line to a file
write.table(t(output_complete_results), file=paste0(PATH_wk,"results/all_AA_raw_complete.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(t(output_simple_results), file=paste0(PATH_wk,"results/all_AA_raw_simple.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")
}


}



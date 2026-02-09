#this script tests the rerun of a single cpg with M1_pval =0, see what happened
#we found that the reason new M1_pval = 0 is because there are not enough samples for these CpGs (too many missing NA values in methylations).
#therefore the pvalues for these cpgs are not reliable (about 2k of cpgs).

cpg_int <- "17:14807842-14807843"

#load the chunk_0001 data first 
##> sel_cpgs
 #  cpg_index                  CpG   M1_coef M1_pval M1_F_sdv M3_coef M3_pval
#       <int>               <char>     <num>   <num>    <num>   <num>   <num>
#1:       818 17:14807842-14807843  24.58517       0 35.87217       0       0
#2:      2900        2:25505881-NA -30.57163       0 49.06305       0       0
#3:      1191       21:46851263-NA -24.82140       0 56.96279       0       0
#   M3_F_sdv M1_lambda M3_lambda Chunk_index
#      <num>     <int>     <int>      <char>
#1:        0         5         0        0690
#2:        0        10         0        0934
#3:        0         5         0        1116


library(data.table)
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

#chunk_dat <- fread(paste0(PATH_wk,"/scr/6_beluga/chunk_0690.csv"), header=FALSE)


#first data chunk
task_id <- 1116



i_cpg <- 1191

###########################################################################
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_scr6 <- paste0(PATH_wk,"scr/6_beluga/")
PATH_scr8 <- paste0(PATH_wk,"scr/8_rerun/")
#first to load one splited data C:\Per\LaiJiang\Project\UQAC\meth\dat\chunk_017.csv
meth_file_loc <- paste0(PATH_wk,"scr/6_beluga/chunk_",sprintf("%04d",task_id),".csv")

#location to save reults
results_loc <-   paste0(PATH_wk,"results/whole_genome/AA_",sprintf("%04d",task_id))

###########################################################################
meth_file <- read.csv(file=meth_file_loc, header=TRUE, sep="\t")

#methfile contains the JobID!!!
#the saving lcoation also contains the JobID!!!

#load the pheno_file from the created data
load(file=paste0(PATH_scr6,"4_1_data.RData"), verbose=TRUE)


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
library(data.table)


rerun_ids <- fread(paste0(PATH_wk,"/scr/8_rerun/results/2_col_rerun.txt"), header=FALSE)

rerun_chunk <- rerun_ids[rerun_ids$V1 == task_id,]

rerun_i_cpg <-rerun_chunk$V2

#################################!!!!!!!!!!!!!

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

print(dim(pheno_file_feed))




######################################################################
#test the differences between the M1 and M2 results 

load( paste0(PATH_wk,"/scr/8_rerun/results/4_eval.RData"), verbose=TRUE)

f<-M1_manhattan[M1_manhattan$pval<1e-5,]
f[grep("7:",f$CpG),]


g<- AA_simple_results[AA_simple_results$M2_pval<1e-5,]
g[grep("7:",g$CpG),]


#this is to debug which error message come from


#load the results from 3_2_run.R
models_res <- read.table(file=paste0(PATH_wk,"dat/AA_raw_simple.txt"),header=FALSE)

#which number is missing from models_res$V1 
all_cpgs <- c(1:max(models_res$V1))
missing_cpgs <- all_cpgs[!(all_cpgs %in% models_res$V1)]


################
models_res_updates <- read.table(file=paste0(PATH_wk,"results/6_AA_raw_simple.txt"),header=FALSE)

unique(models_res_updates$V1)

which(! missing_cpgs %in% models_res_updates$V1)
#everything is in there. 

library(lmerTest)
library(glmmLasso)
library(lme4)
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"


load(file=paste0(PATH_wk,"results/4_1_data.RData"), verbose=TRUE  )

#Approach 1: AA Status as the Outcome:

#model 1: binary logisitic regression

i_cpg <- missing_cpgs[1]

for(i_cpg in missing_cpgs){
###########################################
meth_response_missing <- beta_asr[missing_cpgs,]
#count the number of NA in each row
meth_response_missing_na <- apply(meth_response_missing, 1, function(x) sum(is.na(x)))



###########################################
print(paste0("ignore NA: CpG ",i_cpg))

#the asr transformation
meth_response <- beta_asr[i_cpg,]
#will try original beta_asr later

pheno_file$meth_response <- meth_response

#remove NA rows

pheno_file_feed <- na.omit(pheno_file)
#how many NA in each column


#if there are enough sample sizes: n > p, p = 13 fixed effects and 1 FID
if(nrow(pheno_file_feed) > 14){
source(paste0(PATH_wk,"scr/6_whole/6_0_func_single.R"))

#attach this line to a file
write.table(t(output_complete_results), file=paste0(PATH_wk,"results/6_AA_raw_complete.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(t(output_simple_results), file=paste0(PATH_wk,"results/6_AA_raw_simple.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")
}


}

##################

#how to load one splited data C:\Per\LaiJiang\Project\UQAC\meth\dat\chunk_017.csv
chunk_017 <- read.csv(file="C:/Per/LaiJiang/Project/UQAC/meth/dat/chunk_017.csv", header=TRUE, sep="\t")

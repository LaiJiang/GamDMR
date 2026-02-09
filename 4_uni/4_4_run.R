

library(lmerTest)
library(glmmLasso)
library(lme4)
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"


load(file=paste0(PATH_wk,"results/4_1_data.RData"), verbose=TRUE  )

#Approach 1: AA Status as the Outcome:

#model 1: binary logisitic regression

i_cpg <- 1



###imputed dataset
for(i_cpg in 1:nrow(meth_matrix)){

print(paste0("Impute NA: CpG ",i_cpg))


#the asr transformation
meth_response <- beta_asr_impute[i_cpg,]

pheno_file_feed <- pheno_file

pheno_file_feed$meth_response <- meth_response
#run models 

source(paste0(PATH_wk,"scr/4_uni/4_3_func_single.R"))



#attach this line to a file
write.table(t(output_complete_results), file=paste0(PATH_wk,"results/AA_impute_complete.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(t(output_simple_results), file=paste0(PATH_wk,"results/AA_impute_simple.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")


}

###raw dataset 

for(i_cpg in 1:nrow(meth_matrix)){

print(paste0("ignore NA: CpG ",i_cpg))

#the asr transformation
meth_response <- beta_asr[i_cpg,]
#will try original beta_asr later

pheno_file$meth_response <- meth_response

#remove NA rows

pheno_file_feed <- na.omit(pheno_file)
#run models 

#if there are enough sample sizes
if(nrow(pheno_file_feed)>100){
source(paste0(PATH_wk,"scr/4_uni/4_3_func_single.R"))



#attach this line to a file
write.table(t(output_complete_results), file=paste0(PATH_wk,"results/AA_raw_complete.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(t(output_simple_results), file=paste0(PATH_wk,"results/AA_raw_simple.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")
}

}
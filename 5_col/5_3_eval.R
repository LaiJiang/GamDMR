#evaluate the selected CpGs




PATH_UQAC <- "C:/Per/LaiJiang/Project/UQAC/"

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"



###########################################################################
#save the first three columns to a different csv file

#only need once

if(FALSE){
meth_file <- read.csv(paste0(PATH_wk, "dat/peek_methylation_sva_short.csv"),sep="\t")

write.csv(meth_file[,1:3], file=paste0(PATH_wk, "dat/peek_methylation_sva_short_3col.csv"), row.names=FALSE)
}

###########################################################################


if(FALSE){
load( file=paste0(PATH_wk,"results/meth_matrix.RData"),verbose=TRUE)

##########now order the rows of meth_matrix_impute according to the variance 
var_cpgs_meth<-apply(meth_matrix_imputed,1,var)
#extract the order of var_cpgs_meth from top to bottom
var_cpgs_meth_order <- order(var_cpgs_meth,decreasing = TRUE)
save(var_cpgs_meth_order,file=paste0(PATH_wk,"results/var_cpgs_meth_order.RData"))

}
###########################################################################
#re-order by the ranking of the variance of the CpGs
if(FALSE){
#read file=paste0(PATH_wk, "dat/peek_methylation_sva_short_3col.csv")
CpG_pos <- read.csv(paste0(PATH_wk, "dat/peek_methylation_sva_short_3col.csv"))

#fill in missing values of end position
CpG_pos$end <- CpG_pos$start + 1

load( file=paste0(PATH_wk,"results/var_cpgs_meth_order.RData"),verbose=TRUE)

CpG_pos <- CpG_pos[var_cpgs_meth_order,]

write.csv(CpG_pos, file=paste0(PATH_wk, "dat/5_3_eval_CpGs.csv"), row.names=FALSE)

}
###########################################################################

#Load CpG information containing the chr, start, end positions of CpGs
CpG_pos <- read.csv(paste0(PATH_wk, "dat/5_3_eval_CpGs.csv"))

#load the selected significnat CpGs from model 3 thaat are associated with AA alergic asthma



#load the results from 3_2_run.R
models_res <- read.table(file=paste0(PATH_wk,"dat/AA_impute_simple.txt"),header=FALSE)

colnames(models_res) <- c("CpG",
                     "model1_beta","model1_pvalue", "model1_FID_sd",
                        "model2_beta","model2_pvalue", "model2_FID_sd",
                     "model3_beta","model3_pvalue", "model3_FID_sd")



# Extract significant CpG indices for each model:
sig_model1 <- models_res$CpG[which(models_res$model1_beta!= 0)]
sig_model2 <- models_res$CpG[which(models_res$model2_pvalue < 0.01)]
sig_model3 <- models_res$CpG[which(models_res$model3_beta != 0)]

cpgs_model1 <- CpG_pos[sig_model1,]
cpgs_model2 <- CpG_pos[sig_model2,]
cpgs_model3 <- CpG_pos[sig_model3,]


##now evaluate the selected CpGs, starting with cpgs_model3

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

###################################
CpGs_selection = cpgs_model1


CpGs_selection$chr <- paste0("chr",as.character(CpGs_selection$chr))
names(CpGs_selection) <- c("chr","pos","end")

annotated_cpgs <- merge(CpGs_selection, ann, by.x = c("chr", "pos"), by.y = c("chr", "pos"), all.x = TRUE)

#print these rows in annotated_cpgs that Name is not NA
unique(annotated_cpgs$UCSC_RefGene_Name)

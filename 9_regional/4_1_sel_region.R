#this file select the regions of interst 
library(data.table)
library(dplyr)
library(stringr)

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

#before rerun, the first stage results
AA_simple_results <- fread(paste0(PATH_wk,"/scr/8_rerun/dat/rerun_simple.txt"), header=FALSE)
colnames(AA_simple_results) <- c("Chunk_index","cpg_index","CpG",
                                  "M1_coef","M1_pval","M1_F_sdv",
                                  "M2_coef","M2_pval","M2_F_sdv",
                                  "M3_coef","M3_pval","M3_F_sdv" ,
                                  "M1_lambda","M3_lambda")                                
                                  
#remove dupliated rows
AA_simple_results <- AA_simple_results[!duplicated(AA_simple_results), ]

#replace the NA in M2_coef as 0 
AA_simple_results$M2_coef[is.na(AA_simple_results$M2_coef)] <- 0
AA_simple_results$M2_pval[is.na(AA_simple_results$M2_pval)] <- 1

min(AA_simple_results$M1_pval)
min(AA_simple_results$M3_pval)

#we extract regions that are relaxed:
#M1, M3, lasso coefficient are not 0
#M2: pval < 0.05

#extract the AA_simple_results$CpG that need to be rerun: AA_simple_results$M1_coef !=0 OR AA_simple_results$M3_coef !=0
#OR AA_simple_results$M2_pval< 0.05

#id_rerun <- AA_simple_results[which(AA_simple_results$M1_coef !=0 | AA_simple_results$M3_coef !=0 | AA_simple_results$M2_pval< 0.05/length(AA_simple_results$M2_pval)),]
#
id_rerun <- AA_simple_results[which(AA_simple_results$M1_coef !=0 | AA_simple_results$M3_coef !=0 | AA_simple_results$M2_pval< 0.05),]



id_rerun <- id_rerun %>%
  mutate(
    Chunk_index = Chunk_index %>%
      str_remove("^AA_") %>%             # drop the prefix
      str_remove("_simple\\.txt$") %>%    # drop the suffix
      as.numeric()                        # convert to numeric
  )

#the cpgs of interest
cpgs_sel <- id_rerun$CpG[!duplicated(id_rerun$CpG)]

#save the id_rerun to a file
write.table(
  cpgs_sel,
  file = paste0(PATH_wk, "/scr/9_regional/results/4_sel_region.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

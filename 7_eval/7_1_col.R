#mahattan plot


library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)


#first load the results 

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
library(data.table)
AA_simple_results <- fread(paste0(PATH_wk,"/scr/6_beluga/all_simple.txt"), header=FALSE)
colnames(AA_simple_results) <- c("Chunk_index","CpG",
                                  "M1_coef","M1_pval","M1_F_sdv",
                                  "M2_coef","M2_pval","M2_F_sdv",
                                  "M3_coef","M3_pval","M3_F_sdv" ,
                                  "M1_lambda","M3_lambda")                                
                                  

#remove dupliated rows
AA_simple_results <- AA_simple_results[!duplicated(AA_simple_results), ]

#replace the NA in M2_coef as 0 
AA_simple_results$M2_coef[is.na(AA_simple_results$M2_coef)] <- 0
AA_simple_results$M2_pval[is.na(AA_simple_results$M2_pval)] <- 1

#plot M1_pval only for those M1_coef!=0
M1_pval <- AA_simple_results$M1_pval[AA_simple_results$M1_coef!=0]
M3_pval <- AA_simple_results$M3_pval[AA_simple_results$M3_coef!=0]
M2_pval <- AA_simple_results$M2_pval


sum(M1_pval < 0.05/length(M1_pval))
sum(M3_pval < 0.05/length(M3_pval))
sum(M2_pval < 0.05/length(M2_pval))


AA_simple_results[which(is.na(M2_pval)),]

#extract the AA_simple_results$CpG that need to be rerun: AA_simple_results$M1_coef !=0 OR AA_simple_results$M3_coef !=0
id_rerun <- AA_simple_results[which(AA_simple_results$M1_coef !=0 | AA_simple_results$M3_coef !=0),]


#save the id_rerun to a file
write.table(id_rerun, file = paste0(PATH_wk,"/scr/7_eval/7_1_id_rerun.txt"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



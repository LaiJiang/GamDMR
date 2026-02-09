#mahattan plot


library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)


#first load the results 

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
library(data.table)
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

#plot M1_pval only for those M1_coef!=0
M1_pval <- AA_simple_results$M1_pval[AA_simple_results$M1_coef!=0]
M3_pval <- AA_simple_results$M3_pval[AA_simple_results$M3_coef!=0]
M2_pval <- AA_simple_results$M2_pval


sum(M1_pval < 0.05/length(M1_pval))
sum(M3_pval < 0.05/length(M3_pval))
sum(M2_pval < 0.05/length(M2_pval))

AA_simple_results[which(is.na(M2_pval)),]

sum(AA_simple_results$M2_pval<0.01)

#extract the AA_simple_results$CpG that need to be rerun: AA_simple_results$M1_coef !=0 OR AA_simple_results$M3_coef !=0
id_rerun <- AA_simple_results[which(AA_simple_results$M1_coef !=0 | AA_simple_results$M3_coef !=0),]


#save the id_rerun to a file
if(FALSE){
write.table(id_rerun, file = paste0(PATH_wk,"/scr/7_eval/7_1_id_rerun.txt"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

}


################

#conver the chunk index to the nuemric value 
id_rerun$Chunk_index <- gsub("AA_", "", id_rerun$Chunk_index)
id_rerun$Chunk_index <- gsub("_simple.txt", "", id_rerun$Chunk_index)

id_rerun$Chunk_index <- as.numeric(id_rerun$Chunk_index)

#add a column M1_rerun = 1 if M1_coef !=0, else 0
id_rerun$M1_rerun <- ifelse(id_rerun$M1_coef !=0, 1, 0)
#add a column M3_rerun = 1 if M3_coef !=0, else 0
id_rerun$M3_rerun <- ifelse(id_rerun$M3_coef !=0, 1, 0)


#select columns Chunk_index,cpg_index, CpG, M1_rerun, M3_rerun
id_rerun <- id_rerun[,c("Chunk_index","cpg_index","CpG","M1_rerun","M3_rerun")]

table(id_rerun$M1_rerun,id_rerun$M3_rerun)


#save the id_rerun to a file
write.table(id_rerun, file = paste0(PATH_wk,"/scr/8_rerun/results/2_col_rerun.txt"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##################now extract the id_rerun to re-run M2 model adding BMI
id_rerun <- AA_simple_results[which(AA_simple_results$M2_pval < 0.05),]



#conver the chunk index to the nuemric value 
id_rerun$Chunk_index <- gsub("AA_", "", id_rerun$Chunk_index)
id_rerun$Chunk_index <- gsub("_simple.txt", "", id_rerun$Chunk_index)

id_rerun$Chunk_index <- as.numeric(id_rerun$Chunk_index)

#add a column M2_rerun = 1 if M2_pval <0.05, else 0
id_rerun$M2_rerun <- ifelse(id_rerun$M2_pval <0.05, 1, 0)


#select columns Chunk_index,cpg_index, CpG, M2_rerun
id_rerun <- id_rerun[,c("Chunk_index","cpg_index","CpG","M2_rerun")]



#save the id_rerun to a file
write.table(id_rerun, file = paste0(PATH_wk,"/scr/8_rerun/results/2_col_rerun_M2.txt"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


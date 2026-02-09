
library(data.table)
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
outdir <- paste0(PATH_wk, "results/11_mgcv/")
PATH_scr11 <- paste0(PATH_wk, "/scr/11_mgcv/")


#first load M1, M3, M2 model results and combine them.



#################################################################################################
#M1 M3 results
rerun_simple_results <- fread(paste0(PATH_wk,"/results/12_cpg_DMR_simple.txt"), header=FALSE)



rerun_simple_results[, chunk_id := sub(".*AA_(.*?)_simple\\.txt.*", "\\1", V1)]

#remove the first column
rerun_simple_results <- rerun_simple_results[, -1, with = FALSE]

rerun_simple_results <- rerun_simple_results[!duplicated(rerun_simple_results), ]
colnames(rerun_simple_results) <- c("cpg_index","CpG",
                                  "M1_coef","M1_pval","M1_F_sdv",
                                  "M3_coef","M3_pval","M3_F_sdv" ,
                                  "M1_lambda","M3_lambda", "Chunk_index")  
#################################################################################################
#now load M2 results
rerun_simple_results_M2 <- fread(paste0(PATH_wk,"/results/12_cpg_DMR_M2_simple.txt"), header=FALSE)


rerun_simple_results_M2[, chunk_id := sub(".*AA_(.*?)_simple\\.txt.*", "\\1", V1)]

#remove the first column
rerun_simple_results_M2 <- rerun_simple_results_M2[, -1, with = FALSE]

rerun_simple_results_M2 <- rerun_simple_results_M2[!duplicated(rerun_simple_results_M2), ]
colnames(rerun_simple_results_M2) <- c("cpg_index",
                                  "M2_coef","M2_pval","M2_F_sdv",
                                   "Chunk_index")  
#################################################################################################
#################################################################################################
#################################################################################################
#now load the previous rerun base results

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


#conver the chunk index to the nuemric value 
AA_simple_results$Chunk_index <- gsub("AA_", "", AA_simple_results$Chunk_index)
AA_simple_results$Chunk_index <- gsub("_simple.txt", "", AA_simple_results$Chunk_index)



#################################################################################################

#assign index for merging
AA_simple_results$chunk_cpg_index <- paste0(AA_simple_results$Chunk_index, "_", AA_simple_results$cpg_index)
rerun_simple_results_M2$chunk_cpg_index <- paste0(rerun_simple_results_M2$Chunk_index, "_", rerun_simple_results_M2$cpg_index)
rerun_simple_results$chunk_cpg_index <- paste0(rerun_simple_results$Chunk_index, "_", rerun_simple_results$cpg_index)

#

# Ensure data.tables
setDT(AA_simple_results)
setDT(rerun_simple_results_M2)

# Keep only the columns we need from the rerun table
r2 <- rerun_simple_results_M2[, .(chunk_cpg_index, M2_coef, M2_pval, M2_F_sdv)]

# Update join: overwrite M2_* in AA_simple_results wherever keys match
AA_simple_results[
  r2,
  `:=`(
    M2_coef = i.M2_coef,
    M2_pval = i.M2_pval,
    M2_F_sdv = i.M2_F_sdv
  ),
  on = "chunk_cpg_index"
]

# (Optional) quick sanity check: how many rows were updated?
n_updated <- AA_simple_results[r2, .N, on = "chunk_cpg_index", nomatch = 0L]
cat("Rows updated:", n_updated, "\n")

##########################################################################
#update M1 and M3 
# Keep only the columns we need from the rerun table
r2 <- rerun_simple_results[, .(chunk_cpg_index, M1_coef, M1_pval, M1_F_sdv, M3_coef, M3_pval, M3_F_sdv, M1_lambda, M3_lambda)]

# Update join: overwrite M2_* in AA_simple_results wherever keys match
AA_simple_results[
  r2,
  `:=`(
    M1_coef = i.M1_coef,
    M1_pval = i.M1_pval,
    M1_F_sdv = i.M1_F_sdv,
    M3_coef = i.M3_coef,
    M3_pval = i.M3_pval,
    M3_F_sdv = i.M3_F_sdv,
    M1_lambda = i.M1_lambda,
    M3_lambda = i.M3_lambda
  ),
  on = "chunk_cpg_index"
]

# (Optional) quick sanity check: how many rows were updated?
n_updated <- AA_simple_results[r2, .N, on = "chunk_cpg_index", nomatch = 0L]
cat("Rows updated:", n_updated, "\n")

##########################################################################

#save AA_simple_results
write.table(AA_simple_results, file = paste0(PATH_wk,"/results/12_cpg_DMR/9_combined_BMI.txt"), 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##########################################################################
#now we compare the results with DMA results

CpG_results <- fread(paste0(PATH_wk,"/results/12_cpg_DMR/9_combined_BMI.txt"), header=TRUE)
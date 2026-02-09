## Pre-ranked GSEA for an Asthma pathway using clusterProfiler

# 1. Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)  # sometimes needed for GSEA


PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(qqman)

library(VennDiagram)

##########################################

#stopped here!!!!!
#curate gene list and corresponding pvalues for each gene.

#load M2 results 

##########################################
##########################################
##########################################
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

#generate CpG_pval data frame and run gene annotations 
M2_manhattan <- AA_simple_results[,c("CpG","M2_pval")]
CpG_pval <- M2_manhattan %>%  dplyr::rename(pval = `M2_pval`)
model_name <- "M2"
#run gene extraction function
source(paste0(PATH_wk,"/scr/8_rerun/0_func_genes_rank.R"))
##########################################
##########################################
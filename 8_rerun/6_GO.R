#mahattan plot


library(ggplot2)
library(dplyr)
library(tidyr)
library(qqman)
library(data.table)
library(VennDiagram)

#first load the results 

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"


#load m2 results

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


#load m1 m3 results

load( paste0(PATH_wk,"/scr/8_rerun/results/4_eval.RData"), verbose=TRUE)

#M1_manhattan$pval_min <- pmin(M1_manhattan$pval, M1_manhattan$pval_old)
#M3_manhattan$pval_min <- pmin(M3_manhattan$pval, M3_manhattan$plva_old)

#select cpgs to run gene annotations: with 4 threshold
cpgs_int <- M1_manhattan$CpG[M1_manhattan$pval < 1e-5]


cpgs_int <- M1_manhattan$CpG[M1_manhattan$pval < 0.05/length(M1_manhattan$pval)]
cpgs_int <- M1_manhattan$CpG[M1_manhattan$pval < 1e-5]
cpgs_int <- M1_manhattan$CpG[M1_manhattan$pval < 5e-8]
cpgs_int <- M1_manhattan$CpG[order(M1_manhattan$pval, decreasing = FALSE)[1:500]]
cpgs_int <- M1_manhattan$CpG[order(M1_manhattan$pval, decreasing = FALSE)[1:1000]]

cpgs_int <- M3_manhattan$CpG[order(M3_manhattan$pval, decreasing = FALSE)[1:1000]]

#then try pval_old combined with pval
##################################################################################
##################################################################################
##################################################################################
#genes annotated with these cpgs 
#cpgs_int <- M3_manhattan$CpG[order(M3_manhattan$pval_min, decreasing = FALSE)[1:1000]]

#M1: gene list 
cpgs_int <- M1_manhattan$CpG[M1_manhattan$pval < 1e-5]
cpgs_int <- cpgs_int[!duplicated(cpgs_int)]
#source(paste0(PATH_wk,"/scr/8_rerun/0_func_enrichment.R"))
#use the updated annotations
source(paste0(PATH_wk,"/scr/8_rerun/0_func_enrichment_updates.R"))
#use the updated annotations and expand the cpg hits by window_size to capture  cis-regulatory effects  that extend beyond the gene body.
#windows_size = 500000
#source(paste0(PATH_wk,"/scr/8_rerun/0_func_enrichment_extend.R"))

M1_genes<- gene_map

M1_manhattan[M1_manhattan$CpG=="10:129706544-129706545",]
df_annot[df_annot$genes_symbol=="VARS1",]
M1_manhattan[M1_manhattan$CpG=="6:31745741-NA",]
sort(M1_manhattan$pval)[1:5]
M1_manhattan[which.min(M1_manhattan$pval),]
df_annot[df_annot$CpG=="4:148539626-148539627",]

cpgs_int <- AA_simple_results$CpG[which(AA_simple_results$M2_pval<1e-5)]
cpgs_int <- cpgs_int[!duplicated(cpgs_int)]
source(paste0(PATH_wk,"/scr/8_rerun/0_func_enrichment_updates.R"))
M2_genes<- gene_map

#M3: gene list 
cpgs_int <- M3_manhattan$CpG[M3_manhattan$pval < 1e-5]
cpgs_int <- cpgs_int[!duplicated(cpgs_int)]
source(paste0(PATH_wk,"/scr/8_rerun/0_func_enrichment_updates.R"))

#use the updated annotations and expand the cpg hits by window_size to capture  cis-regulatory effects  that extend beyond the gene body.
#windows_size = 500000
#source(paste0(PATH_wk,"/scr/8_rerun/0_func_enrichment_extend.R"))

M3_genes<- gene_map


save(M1_genes,M2_genes,M3_genes, file = paste0(PATH_wk,"/scr/8_rerun/results/6_GO_genes.RData"))
##################################################################################
##################################################################################
##################################################################################
# dotplot of top terms
dotplot(ego_bp, showCategory=10) + ggtitle("GO BP Enrichment")

# optional: save results as a data.frame / CSV
go_results_df <- as.data.frame(ego_bp)

#write.csv(go_results_df,file = paste0(PATH_wk,"/scr/7_eval/M1_GO.csv"), row.names=FALSE)
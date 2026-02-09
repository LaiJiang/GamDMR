#redraw manhattan plot for model 1 and 3

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

M2_manhattan <- AA_simple_results[,c("CpG","M2_pval", "M2_coef" , "M2_F_sdv")]
M2_manhattan <- M2_manhattan %>%  dplyr::rename(pval = `M2_pval`) %>%  dplyr::rename(coef = `M2_coef`)

min(M2_manhattan$pval)

sort(M2_manhattan$pval,decreasing = FALSE)[1:10]
####################


load(paste0(PATH_wk,"/scr/8_rerun/results/4_eval.RData"),verbose=TRUE)

M1_manhattan$pval[M1_manhattan$CpG=="1:3094655-3094656"]
M1_manhattan$pval[M1_manhattan$CpG=="6:31745741-NA"]

#add all those cpgs in M1_manhattan that are not in M3_manhattan to M3_manhattan, and put their coef as 0 and pval as 1
M3_manhattan <- M3_manhattan %>%
  bind_rows(M1_manhattan[!(M1_manhattan$CpG %in% M3_manhattan$CpG),] %>%
              mutate(coef = 0, pval = 1))



source(paste0(PATH_wk,"scr/7_eval/func_plot_manhattan.R"))


man_plot1 <- man_plot(M1_manhattan, name ="M1_manhattan")
man_plot3 <- man_plot(M3_manhattan, name ="M3_manhattan")
man_plot2 <- man_plot(M2_manhattan, name ="M2_manhattan")


###################

#generate tables for report in results/

load( file = paste0(PATH_wk,"/scr/8_rerun/results/6_GO_genes.RData"),verbose=TRUE)

write.csv(M1_genes, paste0(PATH_wk,"/scr/8_rerun/results/M1_genes.csv"), row.names = FALSE)
write.csv(M2_genes, paste0(PATH_wk,"/scr/8_rerun/results/M2_genes.csv"), row.names = FALSE)
write.csv(M3_genes, paste0(PATH_wk,"/scr/8_rerun/results/M3_genes.csv"), row.names = FALSE)
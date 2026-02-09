#collect tables and fgiures for the supplementary materials of the paper

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_scr14 <- paste0(PATH_wk, "/scr/14_paper/")

library(data.table)
#collect the M1, M2, M3 results after rerun




# Load the data
M1_infor <- fread(paste0(PATH_wk,"/scr/8_rerun/results/7_gsea_M1_genes.csv"))
M2_infor <- fread(paste0(PATH_wk,"/scr/8_rerun/results/7_gsea_M2_genes.csv"))
M3_infor <- fread(paste0(PATH_wk,"/scr/8_rerun/results/7_gsea_M3_genes.csv"))

load(file = paste0(PATH_wk,"/scr/8_rerun/results/4_eval_M123.RData"),verbose=TRUE)

# Extract significant CpGs
M1_cpgs <- M1_manhattan[M1_manhattan$pval < 1e-5, ]
M3_cpgs <- M3_manhattan[M3_manhattan$pval < 1e-5, ]

# Create table_M1 by merging M1_cpgs with M1_infor
table_M1 <- M1_cpgs[, c("CpG", "pval", "coef")]
table_M1 <- merge(table_M1, 
                  M1_infor[, c("CpG", "genes_entrez", "genes_symbol")], 
                  by = "CpG", 
                  all.x = TRUE)

# Reorder columns if desired
table_M1 <- table_M1[, c("CpG", "pval", "coef", "genes_entrez", "genes_symbol")]

# Similarly for Model 2
table_M2 <- M2_cpgs[, c("CpG", "M2_pval", "M2_coef")]
#now rename columns
colnames(table_M2) <- c("CpG", "pval", "coef")
table_M2 <- merge(table_M2, 
                  M2_infor[, c("CpG", "genes_entrez", "genes_symbol")], 
                  by = "CpG", 
                  all.x = TRUE)
table_M2 <- table_M2[, c("CpG", "pval", "coef", "genes_entrez", "genes_symbol")]

# Similarly for Model 3
table_M3 <- M3_cpgs[, c("CpG", "pval", "coef")]
table_M3 <- merge(table_M3, 
                  M3_infor[, c("CpG", "genes_entrez", "genes_symbol")], 
                  by = "CpG", 
                  all.x = TRUE)
table_M3 <- table_M3[, c("CpG", "pval", "coef", "genes_entrez", "genes_symbol")]


#now combine the three tables into a single table, add a column indicating the model
table_M1$model <- "M1"
table_M2$model <- "M2"
table_M3$model <- "M3"
final_table <- rbind(table_M1, table_M2, table_M3)

#crreate a overlap column, if a cpg or gene appears in multiple models. write the model names in the overlap column for that cpg/gene
final_table$overlap <- NA
for(i in 1:nrow(final_table)){
  cpg <- final_table$CpG[i]
  gene <- final_table$genes_symbol[i]
  models_for_cpg <- unique(final_table$model[final_table$CpG == cpg])
  models_for_gene <- unique(final_table$model[final_table$genes_symbol == gene])
  overlap_models <- unique(c(models_for_cpg, models_for_gene))
  if(length(overlap_models) > 1){
    final_table$overlap[i] <- paste(overlap_models, collapse = ";")
  } else {
    final_table$overlap[i] <- NA
  }
}

#now remove the ";NA" in the overlap column
final_table$overlap <- gsub(";NA", "", final_table$overlap)

#now rename the overlap column to "overlap_gene"
colnames(final_table)[which(colnames(final_table) == "overlap")] <- "overlap_gene"

#write to a csv file
fwrite(final_table, file = paste0(PATH_scr14, "results/9_table_significant_CpGs_M1_M2_M3.csv"))
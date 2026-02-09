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

#note this part only need to run once, it take a long time to run the gene annotations
if(FALSE){
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
##########################################

#noew results of model 1 and model 3 (M1_manhattan,M3_manhattan)
load( paste0(PATH_wk,"/scr/8_rerun/results/4_eval.RData"), verbose=TRUE)

#input cogs_int and model_name
model_name <- "M1"
CpG_pval <- M1_manhattan #data set containg CpG and pval
#run gene extraction function
source(paste0(PATH_wk,"/scr/8_rerun/0_func_genes_rank.R"))
#########
#model 3
model_name <- "M3"
CpG_pval <- M3_manhattan #data set containg CpG and pval
#run gene extraction function
source(paste0(PATH_wk,"/scr/8_rerun/0_func_genes_rank.R"))
#########
}
###########################################
###########################################
###########################################

# 2. Read in your CpG-associated gene stats
#    Assume you have a CSV with two columns: 'genes_symbol' and 'pval' 

model_name <- "M3"
cpg_df <- read.csv( paste0(PATH_wk,"/scr/8_rerun/results/7_gsea_",model_name,"_genes.csv"), stringsAsFactors = FALSE)

cpg_df_sel <- cpg_df[cpg_df$pval < 1e-5,]
#show the row with lowest pval
cpg_df[which.min(cpg_df$pval),] #check the cpgs in the region

#show the rows start with 11:xxx at column CpG
cpg_df_sel[grep("^11:", cpg_df_sel$CpG),] #check the cpgs in the region
cpg_df[cpg_df$CpG == "6:31745741-NA",] #check the cpgs in the region

# 3. Map symbols to Entrez IDs
gene_map <- bitr(
  cpg_df$genes_symbol,
  fromType   = "SYMBOL",
  toType     = "ENTREZID",
  OrgDb      = org.Hs.eg.db
)

# 4. Merge and create a pre-ranked vector
merged_df <- merge(cpg_df, gene_map, by.x = "genes_symbol", by.y = "SYMBOL")
merged_df$score <- -log10(merged_df$pval)
geneList <- setNames(merged_df$score, merged_df$ENTREZID)
geneList <- sort(geneList, decreasing = TRUE)

# 5. Run GSEA against KEGG pathways
gsea_kegg <- gseKEGG(
  geneList       = geneList,
  organism       = "hsa",
  nPerm          = 1000,
  minGSSize      = 10,
  maxGSSize      = 500,
  pvalueCutoff   = 0.05,
  pAdjustMethod  = "BH",
  verbose        = TRUE
)

gsea_kegg@result[1:10,"Description"]

###only significant genes

sig_genes <- merged_df$ENTREZID[merged_df$pval < 0.05]
enrich_terms <- enrichKEGG(gene         = sig_genes,
           organism     = "hsa",
           pvalueCutoff = 0.05)

# 1. Convert to a regular data.frame
kegg_df <- as.data.frame(enrich_terms@result)
# 2. Look at the top hits

print(kegg_df[1:10,c("category","subcategory","Description")] )

# 6. Extract the Asthma pathway (KEGG ID: hsa05310)
asthma_res <- subset(kegg_df, ID == "hsa05310")

# 7. Visualize
library(enrichplot)

#save this plot as a jpeg
ggsave(filename = paste0(PATH_wk,"/scr/8_rerun/results/7_gsea_",model_name,"_asthma.jpeg"),
       width = 10, height = 6, dpi = 300)

gseaplot2(
  x       = enrich_terms,
  geneSetID = "hsa05310",
  title     = "GSEA: KEGG Asthma",
  pvalue_table = TRUE
)

dev.off()


# 8. Save results
#write.csv(asthma_res, file = "GSEA_KEGG_Asthma_results.csv", row.names = FALSE)


#tomorrow:  curate the pvalues for corresponding genes from three models seperately, then do GSEA analysis.
#conclusion: the GSEA analysis turns out the asthma pathway (hsa05310) is not enriched in the significant cpgs from any  model 1, 2 or 3. 



#stopped here!!!!
#tomorrow: why the HLA genes not in the list? whats the pvalue look like in the region?

#check beluga the files to be delted.


#  report include the gene list, venndiagram, and gene set enrichment analysis. 
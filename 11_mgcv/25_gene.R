library(data.table)
library(dplyr)
#now extract the genes
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
outdir <- paste0(PATH_wk, "results/11_mgcv/")

DMRs <- fread(   file.path(outdir, "24_dmrs_STRICT.tsv"))

PATH_scr11 <- "C:/Per/LaiJiang/Project/UQAC/meth/scr/11_mgcv/"
# Load region info
region_file <- fread(paste0(PATH_scr11, "dat/region_file_1_chunk.csv"))
region_file$region_id <- 1:nrow(region_file)
# Merge with results for only intersecting region_ids
DMRs <- merge(DMRs, region_file, by = "region_id", all.x = TRUE)

#added pvals, to chek these negative and 0 pvalues should be removed or not?
DMRs <- DMRs %>%
  dplyr::select(region_id, chr, region_start, region_end,pvals)


  source("C:/Per/LaiJiang/Project/UQAC/meth/scr/9_regional/0_func_genes_16.R")


write.csv(df_annot, file = paste0(PATH_wk, "/results/11_mgcv/25_DMR_strict_genes.csv"), row.names = FALSE)

df_annot <- read.csv(paste0(PATH_wk, "/results/11_mgcv/25_DMR_strict_genes.csv"), stringsAsFactors = FALSE)


#replicate starts here!!! 
#df_annot <- read.csv(paste0(PATH_wk, "/results/11_mgcv/11_2_DMR_genes.csv"), stringsAsFactors = FALSE)
genes_int <- c("IL33","PTPRE", "ORMDL3;GSDMB", "IL1RL1", "TSLP", "IL13", "IL4","IL4R","FCER1A", "ATXN7L1;CDHR3","GATA3", "LRP1;STAT6",
"IRAK3","RORA","IKZF3","ADAM33","RUNX1","DPP10", "TP73", "AGRN", "ORMDL3", "ITPR3")


intersect(genes_int, df_annot$genes_symbol)



################################################################################################################
#first run pathyway enrichment analysis
# ---- Required Packages ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db",  "enrichplot"))

# ---- Load Libraries ----
library(clusterProfiler)
library(org.Hs.eg.db)
#library(ReactomePA)
library(enrichplot)
library(tidyverse)

# ---- Load Gene List ----
gene_file <- paste0(PATH_wk, "/results/11_mgcv/25_DMR_strict_genes.csv")
dmr_data <- read.csv(gene_file, stringsAsFactors = FALSE)

# ---- Extract Unique Gene Symbols ----
gene_symbols <- unique(unlist(strsplit(dmr_data$genes_symbol, ";")))
gene_symbols <- gene_symbols[gene_symbols != ""]  # Remove empty strings

# ---- Convert to Entrez IDs ----
gene_entrez <- bitr(gene_symbols, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# ---- KEGG Pathway Enrichment ----
kegg_enrich <- enrichKEGG(gene         = gene_entrez$ENTREZID,
                          organism     = 'hsa',
                          pvalueCutoff = 0.05)

# ---- GO Biological Process Enrichment ----
go_enrich <- enrichGO(gene          = gene_entrez$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      readable      = TRUE)
#no GO enrichments found

#now run Reactome enrichment
############################################################################################################
# Ensure character Entrez IDs and drop NAs
genes <- as.character(gene_entrez$ENTREZID)
genes <- genes[!is.na(genes)]

universe <- as.character(universe_entrez)
universe <- universe[!is.na(universe)]

# (Optional) keep only Entrez IDs that exist in org.Hs.eg.db
valid_ids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
genes    <- intersect(genes, valid_ids)
universe <- intersect(universe, valid_ids)

# Run Reactome enrichment (NOTE: organism = "human")
reactome_enrich <- enrichPathway(
  gene          = unique(genes),
  organism      = "human",
  universe      = unique(universe),
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE   # maps Entrez back to symbols in result
)


######################################################
############################################################################################################

# ---- Visualize Top Pathways ----
dotplot(go_enrich, showCategory=15, title="GO Enrichment - Biological Process")
dotplot(kegg_enrich, showCategory=15, title="KEGG Pathway Enrichment")
dotplot(reactome_enrich, showCategory=15, title="Reactome Pathway Enrichment")

#save dotplots into jpeg files
ggsave(paste0(PATH_wk, "/results/11_mgcv/25_GO_BP_enrichment_dotplots.jpeg"), plot = dotplot(go_enrich, showCategory=15, title="GO Enrichment - Biological Process"), width = 10, height = 6, dpi = 300)
ggsave(paste0(PATH_wk, "/results/11_mgcv/25_KEGG_enrichment_dotplots.jpeg"), plot = dotplot(kegg_enrich, showCategory=15, title="KEGG Pathway Enrichment"), width = 10, height = 6, dpi = 300)
ggsave(paste0(PATH_wk, "/results/11_mgcv/25_Reactome_enrichment_dotplots.jpeg"), plot = dotplot(reactome_enrich, showCategory=15, title="Reactome Pathway Enrichment"), width = 10, height = 6, dpi = 300)

# ---- Save Results ----
write.csv(as.data.frame(go_enrich), paste0(PATH_wk, "/results/11_mgcv/25_GO_BP_enrichment.csv"), row.names = FALSE)
write.csv(as.data.frame(kegg_enrich), paste0(PATH_wk, "/results/11_mgcv/25_KEGG_enrichment.csv"), row.names = FALSE)
write.csv(as.data.frame(reactome_enrich), paste0(PATH_wk, "/results/11_mgcv/25_Reactome_enrichment.csv"), row.names = FALSE)


p <- cnetplot(reactome_enrich, 
         showCategory = 10, 
         foldChange = NULL,   # optional, if you have gene-level FC values
         circular = FALSE, colorEdge = TRUE) +
     ggtitle("Reactome Network ")

ggsave(paste0(PATH_wk, "/results/11_mgcv/25_reactome_network.jpeg"), plot = p, width = 10, height = 6, dpi = 300)



p <- cnetplot(kegg_enrich, 
         showCategory = 10, 
         foldChange = NULL,   # optional, if you have gene-level FC values
         circular = FALSE, colorEdge = TRUE) +
     ggtitle("KEGG Network ")

ggsave(paste0(PATH_wk, "/results/11_mgcv/25_kegg_network.jpeg"), plot = p, width = 10, height = 6, dpi = 300)



#This KEGG enrichment does not provide conclusive evidence that your genes/DMRs are directly related to allergic asthma.
#However, the enrichment of cAMP signaling is supportive of asthma-relevant biology.
##########################################

#for these genes, reverse back to check their original DMRs.!!! why some of the filterings filtered them out?

library(data.table)
library(dplyr)
#now extract the genes
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
outdir <- paste0(PATH_wk, "results/11_mgcv/")

#DMRs <- fread(   file.path(outdir, "16_dmrs_STRICT.tsv"))
#DMRs <- fread(   file.path(outdir, "16_dmrs_base.tsv"))
#DMRs <- fread(paste0(PATH_wk, "results/11_mgcv/16_results_all_batch_filter.txt"))
DMRs <- fread(paste0(PATH_wk, "results/11_mgcv/16_results_all_batch.txt"))




PATH_scr11 <- "C:/Per/LaiJiang/Project/UQAC/meth/scr/11_mgcv/"
# Load region info
region_file <- fread(paste0(PATH_scr11, "dat/region_file_1_chunk.csv"))
region_file$region_id <- 1:nrow(region_file)
# Merge with results for only intersecting region_ids
DMRs <- merge(DMRs, region_file, by = "region_id", all.x = TRUE)


DMRs <- DMRs %>%
  dplyr::select(region_id, chr, region_start, region_end)


  source("C:/Per/LaiJiang/Project/UQAC/meth/scr/9_regional/0_func_genes_16.R")


#replicate starts here!!! 
#df_annot <- read.csv(paste0(PATH_wk, "/results/11_mgcv/11_2_DMR_genes.csv"), stringsAsFactors = FALSE)
genes_int <- c("IL33","PTPRE", "ORMDL3;GSDMB", "IL1RL1", "TSLP", "IL13", "IL4","IL4R","FCER1A", "ATXN7L1;CDHR3","GATA3", "LRP1;STAT6",
"IRAK3","RORA","IKZF3","ADAM33","RUNX1","DPP10", "TP73", "AGRN", "ORMDL3", "ITPR3")


intersect(genes_int, df_annot$genes_symbol)

#merge DMRs with df_annot 
DMRs_genes <- merge(DMRs, df_annot, by = "region_id", all.x = TRUE)

all_results  <- fread(paste0(PATH_wk, "results/11_mgcv/24_results_all_batch.txt"))

DMRs_genes <- merge(DMRs_genes, all_results, by = "region_id", all.x = TRUE)


DMRs_genes$interesting <- ifelse(DMRs_genes$genes_symbol %in% genes_int, 1, 0)


pvals <- DMRs_genes$`s(start):AA_only`
# ---- Basic ggplot2 Scatter Plot ----
p <- ggplot(DMRs_genes, aes(x = `s(start):AA_only`, y = `edf_s(start):AA_only`)) +
  geom_point(aes(color = factor(interesting)), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("0" = "grey", "1" = "red"),
                     labels = c("Not Interesting", "Interesting"),
                     name = "Group") +
  theme_minimal(base_size = 14) +
  labs(title = "DMR Gene Associations",
       x = "pvalue",
       y = "edf")

# ---- Save the Plot ----
ggsave(paste0(PATH_wk, "results/11_mgcv/25_temporary1.png"), plot = p, width = 8, height = 6, dpi = 300)

plot(DMRs_genes$`s(start):AA_only`[DMRs_genes$interesting == 1] )
#pvalue < 0.01?

plot(DMRs_genes$`edf_s(start):AA_only`[DMRs_genes$interesting == 1] )
#edf > 0.5?

#mean_diff > 0.005 Or max_diff > 0.01.

plot(DMRs_genes$mean_diff[DMRs_genes$interesting == 1] ,DMRs_genes$max_diff[DMRs_genes$interesting == 1] )
plot(DMRs_genes$AIC[DMRs_genes$interesting == 1] ,DMRs_genes$R2[DMRs_genes$interesting == 1] )

#R2 > 0.2?

#now filter DMR_genes based on these criteria
DMRs_genes_filtered <- DMRs_genes %>%
    filter(
        `s(start):AA_only` >=0, 
             `s(start):AA_only` < 5e-8, # Adjusted p-value threshold
             `edf_s(start):AA_only` > 1,
             max_diff > 0.02,
             R2 > 0.6)

DMRs_genes_filtered$genes_symbol[DMRs_genes_filtered$interesting == 1]

length(unique(DMRs_genes_filtered$genes_symbol))

#the smallest pvalues? the highest effect sizes?
#highest edf? who are these genes?
library(dplyr)

DMRs_genes_filtered %>%
  tibble::as_tibble() %>%                      # or as.data.frame()
  dplyr::arrange(`s(start):AA_only`,
                 dplyr::desc(max_diff),
                 dplyr::desc(`edf_s(start):AA_only`)) %>%
  dplyr::select(region_id, genes_symbol, `s(start):AA_only`,
                max_diff, `edf_s(start):AA_only`) %>%
  dplyr::slice_head(n = 10)

# ---- Load Gene List ----
dmr_data <- DMRs_genes_filtered
#remove rows with NA in genes_symbol
dmr_data <- dmr_data[!is.na(dmr_data$genes_symbol), ]

min(dmr_data$"s(start):AA_only")

# ---- Extract Unique Gene Symbols ----
gene_symbols <- unique(unlist(strsplit(dmr_data$genes_symbol, ";")))
gene_symbols <- gene_symbols[gene_symbols != ""]  # Remove empty strings

# ---- Convert to Entrez IDs ----
gene_entrez <- bitr(gene_symbols, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# ---- KEGG Pathway Enrichment ----
kegg_enrich <- enrichKEGG(gene         = gene_entrez$ENTREZID,
                          organism     = 'hsa',
                          pvalueCutoff = 0.05)

# ---- GO Biological Process Enrichment ----
go_enrich <- enrichGO(gene          = gene_entrez$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      readable      = TRUE)

library(ReactomePA)
reactome_enrich <- enrichPathway(
  gene      = intersect(gene_entrez$ENTREZID, universe_entrez), # Entrez IDs
  organism  = "human",
  universe  = universe_entrez, 
  pvalueCutoff = 0.05
)

cnetplot(reactome_enrich, 
         showCategory = 10, 
         foldChange = NULL,   # optional, if you have gene-level FC values
         circular = FALSE, colorEdge = TRUE)

#print the top 10 genes with lowest p-values
library(dplyr)

top_genes <- dmr_data %>%
  tibble::as_tibble() %>%
  dplyr::slice_min(`s(start):AA_only`, n = 10, with_ties = TRUE) %>%
  dplyr::select(region_id, genes_symbol, `s(start):AA_only`,
                max_diff, `edf_s(start):AA_only`)

dmr_data$genes_symbol[which(dmr_data$`s(start):AA_only` ==0)]

#conclusion: the top genes with lowest pvalues are most related to 
#these potential confounding factors: weight, height, BMI, cholesterol levels, blood pressure.
#food allergy
#bone tissue density
#protein measurement
#!!!!!!!!!!!!!!!!!!!


#what are these top genes? what phenotypes are they related to? can we adjust them in the model?!!!!

#print the top 10 genes with highest edf
library(dplyr)

top_genes <- dmr_data %>%
  tibble::as_tibble() %>% 
  dplyr::arrange(dplyr::desc(`edf_s(start):AA_only`)) %>%
  dplyr::select(region_id, genes_symbol, `s(start):AA_only`, 
                max_diff, `edf_s(start):AA_only`, R2) %>%
  dplyr::slice_head(n = 10)

top_genes


#heightx3, protein measurement, weight
intersect(genes_int, dmr_data$genes_symbol)


##print the top 10 genes with highest max_diff
top_genes <- dmr_data %>%
  tibble::as_tibble() %>% 
  dplyr::arrange(dplyr::desc(max_diff)) %>%
  dplyr::select(region_id, genes_symbol, `s(start):AA_only`, 
                max_diff, `edf_s(start):AA_only`, R2) %>%
  dplyr::slice_head(n = 10)

top_genes



##print the top 10 genes with highest mean_diff
top_genes <- dmr_data %>%
  tibble::as_tibble() %>% 
  dplyr::arrange(dplyr::desc(mean_diff)) %>%
  dplyr::select(region_id, genes_symbol, `s(start):AA_only`, 
                max_diff, `edf_s(start):AA_only`, R2) %>%
  dplyr::slice_head(n = 10)

top_genes

##filter by highest edf may help, and also. consider filter edf harder, then include the genes with non-zero pvalues. 
#run analysis again.
#tomorrow: also run pathway anayslis for strict and base first. see if the results are similar.
#then try different fitlering strategies to see if that helps.


#mahattan plot


library(ggplot2)
library(dplyr)
library(tidyr)
library(qqman)

library(VennDiagram)

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

M1_cpgs <- AA_simple_results$CpG[AA_simple_results$M2_pval < 1e-5]

##################################################################################
#genes annotated with these cpgs 

#top 1000 cpgs 
cpgs_int <- AA_simple_results$CpG[order(AA_simple_results$M2_pval, decreasing = FALSE)[1:1000]]
M1_cpgs <- cpgs_int[!duplicated(cpgs_int)]



# load libraries
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)



# 1) Parse into a data.frame
df_coords <- do.call(rbind, strsplit(M1_cpgs, "[:-]")) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  setNames(c("chr","start","end")) %>%
  mutate(
    start = as.integer(start),
    end   = as.integer(end),
    CpG   = M1_cpgs,
    # <<< add the "chr" prefix so it matches the TxDb style >>>
    chr   = paste0("chr", chr)
  )

#if end is NA, set it to start + 1
df_coords$end[is.na(df_coords$end)] <- df_coords$start[is.na(df_coords$end)] + 1


# 2) Build GRanges for your CpGs
gr_cpg <- makeGRangesFromDataFrame(df_coords,
  seqnames.field   = "chr",
  start.field      = "start",
  end.field        = "end",
  keep.extra.columns = TRUE
)

# 3) Load the UCSC hg19 knownGene models
txdb     <- TxDb.Hsapiens.UCSC.hg19.knownGene
gr_genes <- genes(txdb)   # has seqnames "chr1","chr2",…

# 4) Find overlaps
hits <- findOverlaps(gr_cpg, gr_genes)

df_hits <- data.frame(
  CpG     = gr_cpg$CpG[queryHits(hits)],
  gene_id = gr_genes$gene_id[subjectHits(hits)],
  stringsAsFactors = FALSE
)

# 5) Map Entrez IDs → gene symbols
df_annot <- df_hits %>%
  left_join(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = unique(df_hits$gene_id),
      keytype = "ENTREZID",
      columns = "SYMBOL"
    ),
    by = c("gene_id" = "ENTREZID")
  ) %>%
  group_by(CpG) %>%
  summarize(
    genes_entrez = paste(unique(gene_id), collapse=";"),
    genes_symbol = paste(unique(SYMBOL), collapse=";")
  )

print(df_annot)

write.csv(df_annot, file = paste0(PATH_wk,"/scr/7_eval/M2_cpgs_genes.csv"), row.names = FALSE)

#now exract the gene names from the cpgs

# — load libraries —
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# — 1) extract unique gene symbols from your df_annot —
#    assume df_annot has columns CpG, genes_entrez, genes_symbol
gene_symbols <- df_annot$genes_symbol %>%
  strsplit(";") %>%        # split the semicolon lists
  unlist() %>%             # flatten
  unique() %>%             # unique symbols
  na.omit()                # drop any NAs

# — 2) map SYMBOL → ENTREZID —
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = gene_symbols,
  keytype = "SYMBOL",
  columns = c("ENTREZID")
)
entrez_ids <- unique(gene_map$ENTREZID)

# — 3) run GO enrichment (Biological Process) —
ego_bp <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",           # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE            # convert IDs → gene symbols in output
)

# — 4) inspect & visualize —
# top 10 terms
head(ego_bp, n = 10)

# dotplot of top terms
dotplot(ego_bp, showCategory=10) + ggtitle("GO BP Enrichment")

# optional: save results as a data.frame / CSV
go_results_df <- as.data.frame(ego_bp)


# optional: save results as a data.frame / CSV
go_results_df <- as.data.frame(ego_bp)
write.csv(go_results_df,file = paste0(PATH_wk,"/scr/7_eval/M2_GO.csv"), row.names=FALSE)
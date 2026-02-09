#this script take input of cpgs results and output the ranked genes and corresponding pvalues 


# load libraries
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)


#input CpG_pval: the dataframe containg CpG and pval 
#and model_name for the output file name

cpgs_int <- CpG_pval$CpG
# 1) Parse into a data.frame
df_coords <- do.call(rbind, strsplit(cpgs_int, "[:-]")) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  setNames(c("chr","start","end")) %>%
  mutate(
    start = as.integer(start),
    end   = as.integer(end),
    CpG   = cpgs_int,
    # <<< add the "chr" prefix so it matches the TxDb style >>>
    chr   = paste0("chr", chr)
  )

#if end is NA, set it to start + 1
df_coords$end[is.na(df_coords$end)] <- df_coords$start[is.na(df_coords$end)] + 1

df_coords[df_coords$start == "31745741", ]

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

print(" 5) Map Entrez IDs → gene symbols")
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

#now attach the pvalues to the genes
# 1. Join your p‐values
library(dplyr)

print(" df_annot")

df_annot <- df_annot %>%
  left_join(
    CpG_pval %>% dplyr::select(CpG, pval),   # force dplyr::select
    by = "CpG"
  ) %>%
  mutate(
    pval = replace_na(pval, 1)                   # require tidyr
  ) %>%
  arrange(pval) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(CpG, genes_entrez, genes_symbol, pval, rank)


library(dplyr)
print(" df_top_by_gene")

# Assume df_annot already has CpG, genes_symbol, pval, rank
df_top_by_gene <- df_annot %>%
  group_by(genes_symbol) %>%
  slice_min(order_by = pval, n = 1, with_ties = FALSE) %>% 
  ungroup()%>%
  arrange(pval)

# Inspect
#df_top_by_gene[df_top_by_gene$genes_symbol == "LOXL4", ]
#df_annot[df_annot$genes_symbol == "LOXL4", ]

write.csv(df_top_by_gene,  file = paste0(PATH_wk,"/scr/8_rerun/results/7_gsea_",model_name,"_genes.csv"), row.names = FALSE)
# load libraries from Biomart annotation 
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(EnsDb.Hsapiens.v75)



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

df_coords[df_coords$CpG=="6:31745741-NA", ]



# 1) Build gr_cpg exactly as before (with "chr6")
gr_cpg <- makeGRangesFromDataFrame(df_coords,
  seqnames.field   = "chr",
  start.field      = "start",
  end.field        = "end",
  keep.extra.columns = TRUE
)

# 2) Build GRanges for your CpGs


# 2) Load GENCODE genes and convert to UCSC style
edb      <- EnsDb.Hsapiens.v75
gr_genes <- genes(edb)
seqlevelsStyle(gr_genes) <- "UCSC"


# 3) Overlap
hits <- findOverlaps(gr_cpg, gr_genes)

# 4) Check that your CpG is now caught
any(gr_cpg$CpG[queryHits(hits)] == "6:31745741-NA")




df_hits <- data.frame(
  CpG     = gr_cpg$CpG[queryHits(hits)],
  gene_id = gr_genes$gene_id[subjectHits(hits)],
  stringsAsFactors = FALSE
)

df_hits[df_hits$CpG=="6:31745741-NA", ] #check the cpgs in the region


# 5) Map ENSG IDs → gene symbols
library(dplyr)
library(AnnotationDbi)

# 5) Map ENSEMBL gene IDs → gene symbols
df_annot <- df_hits %>%
  left_join(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = unique(df_hits$gene_id),
      keytype = "ENSEMBL",
      columns = "SYMBOL"
    ),
    by = c("gene_id" = "ENSEMBL")
  ) %>%
  group_by(CpG) %>%
  summarize(
    genes_ensembl = paste(unique(gene_id), collapse = ";"),
    genes_symbol  = paste(unique(SYMBOL), collapse = ";"),
    .groups       = "drop"
  )

# Inspect
df_annot[df_annot$CpG == "6:31745741-NA", ]


print("gene annotations:")
print(df_annot)

#write.csv(df_annot, file = paste0(PATH_wk,"/scr/7_eval/cpgs_int_genes.csv"), row.names = FALSE)

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
print("top 10 terms:")
print(head(ego_bp, n = 10))


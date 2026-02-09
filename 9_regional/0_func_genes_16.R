# Load libraries
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)

# Your DMRs input
# Assume DMRs is your data.frame with: region_id, chr, region_start, region_end
# Make sure chr column does NOT have "chr" prefix yet
DMRs$chr <- paste0("chr", DMRs$chr)  # Add "chr" prefix

# Build GRanges for DMRs
gr_dmr <- makeGRangesFromDataFrame(DMRs,
  seqnames.field = "chr",
  start.field = "region_start",
  end.field = "region_end",
  keep.extra.columns = TRUE
)

# Load the UCSC hg19 knownGene models
txdb     <- TxDb.Hsapiens.UCSC.hg19.knownGene
gr_genes <- genes(txdb)

# Find overlaps between DMRs and genes
hits <- findOverlaps(gr_dmr, gr_genes)

# Create dataframe of overlaps
df_hits <- data.frame(
  region_id = gr_dmr$region_id[queryHits(hits)],
  gene_id   = gr_genes$gene_id[subjectHits(hits)],
  stringsAsFactors = FALSE
)

# Map Entrez IDs to gene symbols
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
  group_by(region_id) %>%
  summarize(
    genes_entrez = paste(unique(gene_id), collapse = ";"),
    genes_symbol = paste(unique(SYMBOL), collapse = ";"),
    .groups = 'drop'
  )

# Optional: If you have p-values per region (e.g., DMRs$pval), you can join them here:
# DMRs$pval <- YOUR_PVALS   # Add p-values to DMRs if available

# df_annot <- df_annot %>%
#   left_join(DMRs %>% select(region_id, pval), by = "region_id") %>%
#   arrange(pval) %>%
#   mutate(rank = row_number())

# Save the results
#write.csv(df_annot, file = paste0(PATH_wk, "/scr/9_regional/results/DMR_genes.csv"), row.names = FALSE)

# View
head(df_annot)

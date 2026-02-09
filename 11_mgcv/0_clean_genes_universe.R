library(data.table)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGGREST)

## 1) Clean panel universe: split semicolons, keep numeric Entrez IDs
# df_annot$genes_entrez is assumed loaded already
universe_raw <- df_annot$genes_entrez

universe_entrez <- universe_raw |>
  str_split(";") |>
  unlist() |>
  str_trim() |>
  # keep pure numeric Entrez IDs
  (\(x) x[grepl("^[0-9]+$", x)])() |>
  unique()

cat("Universe (raw) length:", length(universe_raw), "\n")
cat("Universe (clean Entrez) length:", length(universe_entrez), "\n")

## 2) Build KEGG-recognized background IDs (strip 'hsa:' prefix)
kegg_bg_ids <- names(KEGGREST::keggList("hsa"))           # 'hsa:xxxx'
kegg_bg_ids <- sub("^hsa:", "", kegg_bg_ids)              # 'xxxx' as character

cat("KEGG background size:", length(kegg_bg_ids), "\n")

## 3) Intersections / diagnostics
cat("Universe ∩ KEGG size:", length(intersect(universe_entrez, kegg_bg_ids)), "\n")
cat("DMR gene set size (mapped):", length(unique(gene_entrez$ENTREZID)), "\n")
cat("DMR ∩ KEGG size:", length(intersect(gene_entrez$ENTREZID, kegg_bg_ids)), "\n")

## 4) Use KEGG-filtered universe
universe_entrez_kegg <- intersect(universe_entrez, kegg_bg_ids)
cat("universe genes entrez:", length(universe_entrez_kegg), "\n")

## 5) Rerun KEGG with custom universe
kegg_enrich_universe <- enrichKEGG(
  gene         = intersect(gene_entrez$ENTREZID, kegg_bg_ids),  # ensure KEGG-mappable
  organism     = "hsa",
  universe     = universe_entrez_kegg,
  pvalueCutoff = 0.05
)

## 6) Quick check
if (nrow(as.data.frame(kegg_enrich_universe)) == 0) {
  message("No KEGG terms after cleaning. Consider also trying Reactome/GO or checking mapping rate.")
} else {
  print(head(as.data.frame(kegg_enrich_universe)[, c("ID","Description","p.adjust","Count")]))
}

## 7) Dotplot (save as JPEG)
if (nrow(as.data.frame(kegg_enrich_universe)) > 0) {
  library(enrichplot)
  library(ggplot2)
  p <- dotplot(kegg_enrich_universe, showCategory = 15) +
       ggtitle("KEGG Enrichment (custom universe, cleaned)")
  ggsave(filename = file.path(PATH_wk, "results/11_mgcv/27_KEGG_enrichment_custom_universe.jpeg"),
         plot = p, width = 10, height = 6, dpi = 300)
}


##conclusion: no KEGG terms after cleaning. 

library(ReactomePA)
reactome_enrich <- enrichPathway(
  gene      = intersect(gene_entrez$ENTREZID, universe_entrez), # Entrez IDs
  organism  = "human",
  universe  = universe_entrez, 
  pvalueCutoff = 0.05
)

library(enrichplot)
library(ggplot2)

dotplot(reactome_enrich, showCategory = 15) +
  ggtitle("Reactome Pathway Enrichment")

cnetplot(reactome_enrich, 
         showCategory = 10, 
         foldChange = NULL,   # optional, if you have gene-level FC values
         circular = FALSE, colorEdge = TRUE)


p <- dotplot(reactome_enrich, showCategory = 20) +
     ggtitle("")  +
         xlab("Gene Ratio") 

ggsave(paste0(PATH_wk, "/results/11_mgcv/0_reactome_enrichment_dot.jpeg"), plot = p, width = 10, height = 6, dpi = 300)


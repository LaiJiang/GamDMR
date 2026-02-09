#this file to highlight the Focal adhesion / Rap1 signaling pathway for the AIRS “transport” theme
library(data.table)
# assumed from your pipeline
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

dmr_df   <- fread(file.path(PATH_wk, "results/11_mgcv/25_DMR_strict_genes.csv"))
panel_df <- fread(file.path(PATH_wk, "results/11_mgcv/27_panel_genes.csv"))


#1) Freeze inputs (reproducibility)

#Lock the exact DMR gene list you’ll use (the same one behind your KEGG/Reactome results).

#Lock the panel/universe Entrez IDs you used for enrichment.

#Keep both “with BMI” as primary and “without BMI” as a sensitivity set (optional: “with vs without IgE”).

# Universe (clean Entrez)
library(stringr); library(dplyr); library(org.Hs.eg.db); library(AnnotationDbi)
panel_symbols <- unique(str_trim(unlist(strsplit(panel_df$genes_symbol, ";"))))
map_universe <- AnnotationDbi::select(org.Hs.eg.db, keys=panel_symbols,
                                      keytype="SYMBOL", columns=c("ENTREZID","SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>% distinct()
universe_entrez <- unique(map_universe$ENTREZID)

# DMR genes (clean Entrez)
dmr_symbols <- unique(str_trim(unlist(strsplit(dmr_df$genes_symbol, ";"))))
map_dmr <- AnnotationDbi::select(org.Hs.eg.db, keys=dmr_symbols,
                                 keytype="SYMBOL", columns=c("ENTREZID","SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>% distinct()
dmr_entrez <- unique(intersect(map_dmr$ENTREZID, universe_entrez))


#2) Pull KEGG gene sets for the two pathways and intersect with your DMRs

#KEGG IDs: hsa04510 (Focal adhesion) and hsa04015 (Rap1 signaling).

#Map to Entrez → Symbols; intersect with dmr_entrez.

library(KEGGREST)

kegg_to_df <- function(pid) {
  k <- keggGet(pid)[[1]]
  genes <- k$GENE                               # vector: ENTREZ, desc alternating
  entrez <- genes[seq(1, length(genes), 2)]     # keep Entrez only
  data.frame(PathwayID = pid, ENTREZID = entrez, stringsAsFactors = FALSE)
}

kegg_df <- bind_rows(kegg_to_df("hsa04510"), kegg_to_df("hsa04015"))

# Map to symbols
kegg_df <- AnnotationDbi::select(org.Hs.eg.db, keys=kegg_df$ENTREZID,
                                 keytype="ENTREZID", columns="SYMBOL") %>%
  left_join(kegg_df, by="ENTREZID") %>% distinct()

# Intersect with your DMR gene set
kegg_hits <- kegg_df %>%
  mutate(In_DMR = ENTREZID %in% dmr_entrez) %>%
  filter(In_DMR)

# Quick counts (for abstract)
table(kegg_hits$PathwayID)
length(unique(kegg_hits$ENTREZID))

#3) Build a transport-focused table to show in the talk

#Include pathway, gene, and (if available) region-level stats: region_id, AA effect direction/size, p/FDR.

# OPTIONAL: join region-level stats if you have them (edit column names to your result schema)
# Example assumes you have region_id, genes_symbol, AA_only_coef, AA_only_p
# all_results <- fread(file.path(PATH_wk, "results/11_mgcv/24_results_all_batch.txt"))
# Map regions -> genes -> Entrez, then filter to kegg_hits$ENTREZID


transport_tbl <- unique(kegg_hits[, c("PathwayID", "ENTREZID", "SYMBOL")])

fwrite(transport_tbl, file.path(PATH_wk, "results/11_mgcv/31_airs_transport_table.csv"))


#4) Do a focused enrichment test (hypergeometric) just for these pathways

#Report simple over-representation numbers that match the theme.

hyper_wrap <- function(path_id) {
  set_genes <- kegg_df %>% filter(PathwayID == path_id) %>% pull(ENTREZID) %>% unique()
  K <- length(intersect(set_genes, dmr_entrez))            # overlap
  M <- length(intersect(set_genes, universe_entrez))       # set size in universe
  N <- length(universe_entrez)                             # universe size
  n <- length(dmr_entrez)                                  # DMR set size
  # one-sided hypergeometric (enrichment)
  p <- phyper(q = K-1, m = M, n = N-M, k = n, lower.tail = FALSE)
  data.frame(PathwayID=path_id, Overlap=K, SetSize=M, DMRsize=n, Universe=N, Pvalue=p)
}

enrich_focus <- bind_rows(hyper_wrap("hsa04510"), hyper_wrap("hsa04015")) %>%
  arrange(Pvalue)
enrich_focus


#5) Make figures that won’t get you scooped

#1 slide Manhattan plot (whole-genome EWAS) with no gene labels.

#1 slide transport figure:

#Dotplot highlighting only “Focal adhesion / Rap1 / Leukocyte migration”.

#CNET (gene–pathway network) limited to those two KEGG pathways.

#Per-gene lollipop/bar: AA effect size (or Δmethylation) for top 5 integrin/Rap1 genes.

#1 backup slide: sensitivity (with vs without BMI/IgE) showing pathway persists.

library(clusterProfiler); library(enrichplot); library(ggplot2)

# Rebuild a tiny enrichResult for just these two pathways (optional, for dotplot aesthetics)
# Or just filter your existing KEGG/GO results to the transport terms and plot.

# Example cnetplot limited to your two KEGG pathways:
# Suppose `kegg_enrich_full` is your previous KEGG enrichResult (custom universe applied)
transport_ids <- c("hsa04510","hsa04015")


# Load KEGG results from CSV
kegg_enrich_df <- read.csv(paste0(PATH_wk, "/results/11_mgcv/25_KEGG_enrichment.csv"))


# Convert back to enrichResult-like object
kegg_enrich_full <- new("enrichResult",
                   result = kegg_enrich_df,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.2,
                   organism = "hsa",
                   keytype = "ENTREZID")

kegg_transport <- kegg_enrich_full
kegg_transport@result <- subset(as.data.frame(kegg_enrich_full), ID %in% transport_ids)

p_dot <- dotplot(kegg_transport, showCategory=2) + ggtitle("Transport-focused KEGG pathways")
ggsave(file.path(PATH_wk, "results/11_mgcv/31_fig_transport_kegg_dot.jpeg"), p_dot, width=7, height=5, dpi=300)

# Network (genes↔pathways) – limit to 10–15 genes to avoid clutter
#p_cnet <- cnetplot(kegg_transport, showCategory=2, circular=FALSE)
#ggsave(file.path(PATH_wk, "results/11_mgcv/31_fig_transport_kegg_cnet.jpeg"), p_cnet, width=7, height=5, dpi=300)


#Extra: one-liner checks you can report
# How many DMR genes hit transport pathways?
length(unique(kegg_hits$SYMBOL))

# Which transport genes recur across KEGG+GO?
go_bp <- fread(file.path(PATH_wk, "results/11_mgcv/topGO/topGO_BP_results.csv"))
transport_go <- go_bp %>% filter(grepl("leukocyte migration|cell adhesion", Term, ignore.case=TRUE))
transport_go[1:5, c("GO.ID","Term","Significant","weight01Fisher")]

# Sensitivity: does Focal adhesion / Rap1 stay enriched with/without BMI (repeat your enrichments and compare IDs/FDR)

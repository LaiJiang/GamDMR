suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(data.table)
  library(stringr)
  library(scales)
  library(org.Hs.eg.db)
  library(ReactomePA)
})

base_dir <- Sys.getenv("METH_BASE_DIR", unset = getwd())
results14 <- file.path(base_dir, "results", "14_paper")
results11 <- file.path(base_dir, "results", "11_mgcv")
results10 <- file.path(base_dir, "results", "10_sanity")
results9 <- file.path(base_dir, "results", "9_regional")

dir.create(results14, recursive = TRUE, showWarnings = FALSE)

kegg_enrich_df <- read.csv(file.path(results11, "25_KEGG_enrichment.csv"))
kegg_enrich <- new(
  "enrichResult",
  result = kegg_enrich_df,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  organism = "hsa",
  keytype = "ENTREZID"
)
kegg_enrich@result$Description <- str_wrap(kegg_enrich@result$Description, width = 28)

p1 <- dotplot(kegg_enrich, showCategory = 15) +
  xlab("Gene Ratio") +
  ylab(NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(file.path(results14, "12_1_KEGG_enrichment.pdf"), plot = p1, width = 107, height = 125, units = "mm")

gene_file <- file.path(results11, "25_DMR_strict_genes.csv")
dmr_data <- read.csv(gene_file, stringsAsFactors = FALSE)
gene_symbols <- unique(unlist(strsplit(dmr_data$genes_symbol, ";")))
gene_symbols <- gene_symbols[gene_symbols != "" & !is.na(gene_symbols)]
gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

panel_df <- read.csv(file.path(results11, "27_panel_genes.csv"), stringsAsFactors = FALSE)
universe_raw <- panel_df$genes_entrez
universe_entrez <- unique(trimws(unlist(strsplit(universe_raw, ";"))))
universe_entrez <- universe_entrez[grepl("^[0-9]+$", universe_entrez)]

reactome_enrich <- enrichPathway(
  gene = intersect(gene_entrez$ENTREZID, universe_entrez),
  organism = "human",
  universe = universe_entrez,
  pvalueCutoff = 0.05
)
reactome_enrich@result$Description <- str_wrap(reactome_enrich@result$Description, width = 28)

p2 <- dotplot(reactome_enrich, showCategory = 20) +
  xlab("Gene Ratio") +
  ylab(NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(file.path(results14, "12_2_reactome_enrichment.pdf"), plot = p2, width = 107, height = 145, units = "mm")

i_gene <- 1
load(file.path(results10, paste0("meth_dat_", i_gene, ".RData")))
all_df <- read.csv(file.path(results9, "16_eval_dat.csv"))

i_df <- all_df |>
  dplyr::filter(chr == i_position$Chromosome) |>
  dplyr::filter(
    (region_start < i_position$gene_start & region_end > i_position$gene_start) |
    (region_start < i_position$gene_end & region_end > i_position$gene_end) |
    (region_start > i_position$gene_start & region_end < i_position$gene_end)
  )

i_meth_long <- meth_long |>
  dplyr::filter(start >= i_df$region_start & start <= i_df$region_end)

p3 <- ggplot(i_meth_long, aes(x = start, y = meth, color = factor(AA_only))) +
  geom_point(alpha = 0.3, size = 0.3) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = 0.8) +
  xlab("Genomic Position") +
  ylab("Methylation Level") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(file.path(results14, paste0("12_3_region_", i_gene, ".pdf")), plot = p3, width = 107, height = 85, units = "mm")

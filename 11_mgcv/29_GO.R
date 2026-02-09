###############################
# GO Enrichment on DMR Genes
# (paper-style, topGO + Bioconductor)
###############################
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

## Packages (install if needed)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(topGO)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

## ====== INPUTS ======
# Folder prefix (edit if needed)
# PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth"

# 1) DMR gene list (from your pipeline; must contain 'genes_symbol' column)
dmr_gene_file  <- file.path(PATH_wk, "results/11_mgcv/25_DMR_strict_genes.csv")

# 2) Panel gene annotation (universe): a table with 'genes_symbol' (and ideally 'genes_entrez')
# If you already have df_annot from earlier code, you can skip reading:
panel_gene_file <- file.path(PATH_wk, "results/11_mgcv/27_panel_genes.csv")

# Output folder
outdir <- file.path(PATH_wk, "results/11_mgcv/topGO")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ====== LOAD DATA ======
dmr_df <- fread(dmr_gene_file)    # columns include: genes_symbol (semicolon-separated)
panel_df <- fread(panel_gene_file) # columns include: genes_symbol, genes_entrez (possibly semicolon-separated)

## ====== CLEAN + BUILD UNIVERSE ======
# Expand panel gene symbols (some rows may contain "A;B")
panel_symbols <- panel_df$genes_symbol %>%
  str_split(";") %>% unlist() %>% str_trim()
panel_symbols <- panel_symbols[nzchar(panel_symbols)]
panel_symbols <- unique(panel_symbols)

# Map SYMBOL -> ENTREZ (universe)
map_universe <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = unique(panel_symbols),
  columns  = c("ENTREZID","SYMBOL"),
  keytype  = "SYMBOL"
) %>% filter(!is.na(ENTREZID)) %>% distinct()

# Universe as ENTREZ IDs (named by SYMBOL for readability)
universe_entrez <- unique(map_universe$ENTREZID)

## ====== EXTRACT DMR-SELECTED GENES ======
dmr_symbols <- dmr_df$genes_symbol %>%
  str_split(";") %>% unlist() %>% str_trim()
dmr_symbols <- unique(dmr_symbols[nzchar(dmr_symbols)])

# Map DMR symbols to ENTREZ; keep only those present in universe
map_dmr <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = dmr_symbols,
  columns  = c("ENTREZID","SYMBOL"),
  keytype  = "SYMBOL"
) %>% filter(!is.na(ENTREZID)) %>% distinct()

dmr_entrez <- intersect(unique(map_dmr$ENTREZID), universe_entrez)

message("Universe size (Entrez): ", length(universe_entrez))
message("DMR gene set size (Entrez in universe): ", length(dmr_entrez))

## ====== MAKE topGO INPUTS ======
# topGO wants: a named numeric/binary vector over the *universe* (names = ENTREZ),
# values = 1 if gene is "selected", else 0
allGenes <- as.integer(universe_entrez %in% dmr_entrez)
names(allGenes) <- universe_entrez

# Helper function for one ontology
run_topgo <- function(ontology = c("BP","MF","CC"), nodeSize = 10) {
  ontology <- match.arg(ontology)

  # Map ENTIRE universe ENTREZ -> GO terms for the ontology
  gene2GO_df <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys     = universe_entrez,
    keytype  = "ENTREZID",
    columns  = c("GO","ONTOLOGY")
  ) %>% filter(!is.na(GO), ONTOLOGY == ontology) %>% distinct()

  # Build gene2GO list (Entrez -> vector of GO IDs)
  gene2GO_list <- split(gene2GO_df$GO, gene2GO_df$ENTREZID)

  # topGO object
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = allGenes,
                annot    = annFUN.gene2GO,
                gene2GO  = gene2GO_list,
                nodeSize = nodeSize)

  # Two commonly reported tests: Fisher classic & weight01
  fisher_classic  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  fisher_weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

  # Produce results table
  res_tab <- GenTable(GOdata,
                      classicFisher  = fisher_classic,
                      weight01Fisher = fisher_weight01,
                      orderBy        = "weight01Fisher",
                      topNodes       = min(200, length(score(fisher_weight01)))) %>%
    as.data.frame()

  # Add term names & p-values as numeric
  res_tab <- res_tab %>%
    mutate(classicFisher_num  = suppressWarnings(as.numeric(classicFisher)),
           weight01Fisher_num = suppressWarnings(as.numeric(weight01Fisher)))

  list(GOdata = GOdata,
       table  = res_tab)
}
###################fix

## --- PREAMBLE: make sure these exist from your earlier steps ---
## universe_entrez : character Entrez IDs of panel genes
## dmr_entrez      : character Entrez IDs of DMR genes (subset of universe)
## org.Hs.eg.db    : loaded

suppressPackageStartupMessages({
  library(dplyr)
  library(topGO)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(data.table)
})

# 0) Coerce to clean character IDs
universe_entrez <- unique(as.character(universe_entrez))
dmr_entrez      <- unique(as.character(dmr_entrez))

# 1) Build topGO gene list: 1 for selected (DMR), 0 otherwise
geneList <- as.integer(universe_entrez %in% dmr_entrez)
names(geneList) <- universe_entrez

# Selection function: “selected” if value == 1
selFun <- function(allScore) allScore == 1

message("Universe size: ", length(universe_entrez))
message("DMR in universe: ", sum(geneList))

# 2) Helper to run one ontology safely
run_topgo <- function(ontology = c("BP","MF","CC"), nodeSize = 10) {
  ontology <- match.arg(ontology)

  # Map ENTIRE universe to GO terms for the chosen ontology
  m <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = universe_entrez,
    keytype = "ENTREZID",
    columns = c("GO","ONTOLOGY")
  ) %>%
    filter(!is.na(GO), ONTOLOGY == ontology) %>%
    distinct()

  # Build gene2GO (Entrez -> vector of GO IDs)
  gene2GO <- split(m$GO, m$ENTREZID)

  # Sanity checks
  overlap_ids <- intersect(names(gene2GO), names(geneList))
  if (length(overlap_ids) < 10L) {
    stop(sprintf("Too few genes map to GO (%s ontology). Overlap with universe: %d",
                 ontology, length(overlap_ids)))
  }
  message(sprintf("[%s] genes with GO terms in universe: %d",
                  ontology, length(overlap_ids)))

  # Construct topGOdata (now with geneSelectionFun)
  GOdata <- new("topGOdata",
                ontology         = ontology,
                allGenes         = geneList,
                geneSelectionFun = selFun,
                annot            = annFUN.gene2GO,
                gene2GO          = gene2GO,
                nodeSize         = nodeSize)

  # Run tests
  fisher_classic  <- runTest(GOdata, algorithm = "classic",  statistic = "fisher")
  fisher_weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

  # Results table
  res <- GenTable(GOdata,
                  classicFisher  = fisher_classic,
                  weight01Fisher = fisher_weight01,
                  orderBy        = "weight01Fisher",
                  topNodes       = min(200, length(score(fisher_weight01)))) %>%
    as.data.frame()

  # Numeric p-values + FDR
  res$classicFisher_num  <- suppressWarnings(as.numeric(res$classicFisher))
  res$weight01Fisher_num <- suppressWarnings(as.numeric(res$weight01Fisher))
  res$weight01_FDR       <- p.adjust(res$weight01Fisher_num, method = "BH")

  list(GOdata = GOdata, table = res)
}

# 3) Run all three
set.seed(1)
res_BP <- run_topgo("BP", nodeSize = 10)
res_MF <- run_topgo("MF", nodeSize = 10)
res_CC <- run_topgo("CC", nodeSize = 10)

# 4) Save
outdir <- file.path(PATH_wk, "results/11_mgcv/topGO")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
fwrite(res_BP$table, file.path(outdir, "topGO_BP_results.csv"))
fwrite(res_MF$table, file.path(outdir, "topGO_MF_results.csv"))
fwrite(res_CC$table, file.path(outdir, "topGO_CC_results.csv"))

# 5) Quick peek: top hits by FDR
print(head(res_BP$table[order(res_BP$table$weight01_FDR), c("GO.ID","Term","Annotated","Significant","Expected","weight01Fisher","weight01_FDR")], 10))
print(head(res_MF$table[order(res_MF$table$weight01_FDR), c("GO.ID","Term","Annotated","Significant","Expected","weight01Fisher","weight01_FDR")], 10))
print(head(res_CC$table[order(res_CC$table$weight01_FDR), c("GO.ID","Term","Annotated","Significant","Expected","weight01Fisher","weight01_FDR")], 10))

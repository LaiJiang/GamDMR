#!/usr/bin/env Rscript

# =============================================================================
# 5_collect_compare_read_depth_GAM_gene.R
#
# Purpose
# -------
# Starting from the STRICT DMR results produced by
# 4_collect_compare_read_depth_weighted_GAM*.R, further filter BOTH:
#
#   1) original unweighted strict GAM-DMRs
#   2) read-depth-weighted strict GAM-DMRs
#
# by requiring overlap with at least one annotated gene, using the same rule as
# the historical annotation code:
#
#   gr_genes <- genes(txdb)
#   hits <- findOverlaps(gr_dmr, gr_genes)
#   map Entrez gene IDs -> SYMBOL with org.Hs.eg.db
#   group by region_id
#   keep only DMRs represented in the overlap table
#
# Genome build:
#   hg19, matching the methylation/region coordinates used in this project.
#
# IMPORTANT
# ---------
# This script DOES NOT redefine DMR significance. It reads DMR_weighted and
# DMR_unweighted from the strict comparison output and then applies ONLY the
# gene-overlap filter. It also independently verifies that the DMR flags are
# consistent with the exact historical strict rule:
#
#   p_s_start_AA < 1e-5
#   FDR_s_start_AA < 0.01
#   edf_s_start_AA >= 1
#   N_cpgs >= 10
#   R2 >= 0.30
#   mean_diff >= 0.01
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
})

data.table::setDTthreads(1L)

# ------------------------------- configuration -------------------------------

PATH_wk <- path.expand(
  Sys.getenv(
    "PATH_WK",
    "~/scratch/UQAC/meth/"
  )
)

comparison_dir <- path.expand(
  Sys.getenv(
    "WEIGHTED_GAM_COMPARISON_OUTPUT",
    file.path(
      PATH_wk,
      "results/15_revision/7_read_depth_weighted_GAM/comparison"
    )
  )
)

table_dir_4 <- file.path(
  comparison_dir,
  "tables"
)

weighted_file <- path.expand(
  Sys.getenv(
    "WEIGHTED_STRICT_RESULTS",
    file.path(
      table_dir_4,
      "weighted_GAM_all_regions.tsv"
    )
  )
)

unweighted_file <- path.expand(
  Sys.getenv(
    "UNWEIGHTED_STRICT_RESULTS",
    file.path(
      table_dir_4,
      "unweighted_GAM_all_regions_reconstructed.tsv"
    )
  )
)

region_file_path <- path.expand(
  Sys.getenv(
    "REGION_FILE",
    file.path(
      PATH_wk,
      "scr/11_mgcv/dat/region_file_1_chunk.csv"
    )
  )
)

output_dir <- path.expand(
  Sys.getenv(
    "GENE_FILTER_OUTPUT",
    file.path(
      comparison_dir,
      "gene_filter"
    )
  )
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ----------------------------- strict-rule constants --------------------------

# These are used ONLY to verify the DMR flags inherited from step 4.
STRICT_P <- 1e-5
STRICT_FDR <- 0.01
STRICT_EDF <- 1
STRICT_N_CPG <- 10L
STRICT_R2 <- 0.30
STRICT_MEAN_DIFF <- 0.01

cat(
  paste0(
    "\nGene-overlap filtering of STRICT GAM-DMRs\n",
    "Strict DMR definition expected from step 4:\n",
    "  p_s_start_AA < 1e-5\n",
    "  FDR_s_start_AA < 0.01\n",
    "  edf_s_start_AA >= 1\n",
    "  N_cpgs >= 10\n",
    "  R2 >= 0.30\n",
    "  mean_diff >= 0.01\n\n"
  )
)

# ---------------------------------- helpers ----------------------------------

normalize_chr_ucsc <- function(x) {

  x <- as.character(x)
  x <- sub(
    "^chr",
    "",
    x,
    ignore.case = TRUE
  )

  # Handle occasional numeric coding for sex chromosomes if present.
  x[
    x == "23"
  ] <- "X"

  x[
    x == "24"
  ] <- "Y"

  x[
    x %in% c(
      "25",
      "M",
      "MT"
    )
  ] <- "M"

  paste0(
    "chr",
    x
  )
}

strict_rule <- function(dt) {

  required <- c(
    "p_s_start_AA",
    "FDR_s_start_AA",
    "edf_s_start_AA",
    "N_cpgs",
    "R2",
    "mean_diff"
  )

  missing <- setdiff(
    required,
    names(dt)
  )

  if (length(missing)) {
    stop(
      "Strict-rule verification failed because columns are missing: ",
      paste(
        missing,
        collapse = ", "
      )
    )
  }

  !is.na(dt$p_s_start_AA) &
    dt$p_s_start_AA < STRICT_P &
    !is.na(dt$FDR_s_start_AA) &
    dt$FDR_s_start_AA < STRICT_FDR &
    !is.na(dt$edf_s_start_AA) &
    dt$edf_s_start_AA >= STRICT_EDF &
    !is.na(dt$N_cpgs) &
    dt$N_cpgs >= STRICT_N_CPG &
    !is.na(dt$R2) &
    dt$R2 >= STRICT_R2 &
    !is.na(dt$mean_diff) &
    dt$mean_diff >= STRICT_MEAN_DIFF
}

check_dmr_flag <- function(dt, flag_col, label) {

  if (!flag_col %in% names(dt)) {
    stop(
      label,
      " does not contain required DMR flag column: ",
      flag_col
    )
  }

  expected <- strict_rule(dt)

  observed <- as.logical(
    dt[[flag_col]]
  )

  observed[
    is.na(observed)
  ] <- FALSE

  mismatch <- which(
    expected != observed
  )

  cat(
    label,
    ": ",
    sum(observed),
    " strict DMRs inherited from step 4.\n",
    sep = ""
  )

  if (length(mismatch)) {
    stop(
      label,
      " contains ",
      length(mismatch),
      " rows where ",
      flag_col,
      " does not match the exact historical strict-DMR definition. ",
      "Please regenerate step 4 with the exact strict collector before ",
      "running the gene-overlap step."
    )
  }

  invisible(TRUE)
}

prepare_dmr_coordinates <- function(
  result_dt,
  flag_col,
  region_dt,
  label
) {

  dmr <- copy(
    result_dt[
      as.logical(
        get(flag_col)
      ) %in% TRUE
    ]
  )

  if (!nrow(dmr)) {
    warning(
      "No strict DMRs found for ",
      label,
      "."
    )

    return(dmr)
  }

  if (anyDuplicated(dmr$region_id)) {
    stop(
      "Duplicated region_id values detected among ",
      label,
      " strict DMRs."
    )
  }

  # Keep region coordinates from the original region definition file.
  coord_cols <- c(
    "region_id",
    "chr",
    "region_start",
    "region_end",
    "data_chunk_id"
  )

  coord_cols <- intersect(
    coord_cols,
    names(region_dt)
  )

  region_coord <- unique(
    region_dt[
      ,
      ..coord_cols
    ],
    by = "region_id"
  )

  # Avoid duplicate data_chunk_id column after joining.
  if (
    "data_chunk_id" %in% names(dmr) &&
      "data_chunk_id" %in% names(region_coord)
  ) {
    region_coord[
      ,
      data_chunk_id := NULL
    ]
  }

  dmr <- merge(
    dmr,
    region_coord,
    by = "region_id",
    all.x = TRUE,
    sort = FALSE
  )

  missing_coord <- dmr[
    is.na(chr) |
      is.na(region_start) |
      is.na(region_end),
    region_id
  ]

  if (length(missing_coord)) {
    stop(
      "Could not recover genomic coordinates for ",
      length(missing_coord),
      " ",
      label,
      " DMRs. Example region_id(s): ",
      paste(
        head(
          missing_coord,
          10L
        ),
        collapse = ", "
      )
    )
  }

  dmr[
    ,
    chr_ucsc := normalize_chr_ucsc(chr)
  ]

  dmr
}

annotate_gene_overlaps <- function(
  dmr_dt,
  gr_genes,
  label
) {

  if (!nrow(dmr_dt)) {
    return(
      list(
        annotation = data.table(),
        filtered = data.table()
      )
    )
  }

  gr_dmr <- GRanges(
    seqnames = dmr_dt$chr_ucsc,
    ranges = IRanges(
      start = as.integer(
        dmr_dt$region_start
      ),
      end = as.integer(
        dmr_dt$region_end
      )
    ),
    region_id = as.character(
      dmr_dt$region_id
    )
  )

  # Keep only sequence names shared with the hg19 TxDb.
  common_seq <- intersect(
    seqlevels(gr_dmr),
    seqlevels(gr_genes)
  )

  if (!length(common_seq)) {
    stop(
      "No common chromosome names between ",
      label,
      " DMRs and hg19 TxDb."
    )
  }

  gr_dmr <- keepSeqlevels(
    gr_dmr,
    common_seq,
    pruning.mode = "coarse"
  )

  hits <- findOverlaps(
    gr_dmr,
    gr_genes,
    ignore.strand = TRUE
  )

  if (!length(hits)) {
    warning(
      "No gene overlaps found for ",
      label,
      " DMRs."
    )

    return(
      list(
        annotation = data.table(),
        filtered = dmr_dt[0]
      )
    )
  }

  # This reproduces the old code logic:
  #   region_id = gr_dmr$region_id[queryHits(hits)]
  #   gene_id   = gr_genes$gene_id[subjectHits(hits)]
  gene_ids <- if (
    "gene_id" %in% names(
      mcols(gr_genes)
    )
  ) {
    as.character(
      gr_genes$gene_id
    )
  } else {
    as.character(
      names(gr_genes)
    )
  }

  df_hits <- data.table(
    region_id = as.integer(
      mcols(gr_dmr)$region_id[
        queryHits(hits)
      ]
    ),
    gene_id = as.character(
      gene_ids[
        subjectHits(hits)
      ]
    )
  )

  df_hits <- unique(
    df_hits[
      !is.na(region_id) &
        !is.na(gene_id)
    ]
  )

  # Map Entrez IDs to HGNC symbols exactly as in the historical code.
  gene_map <- as.data.table(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(
        df_hits$gene_id
      ),
      keytype = "ENTREZID",
      columns = "SYMBOL"
    )
  )

  setnames(
    gene_map,
    "ENTREZID",
    "gene_id"
  )

  df_hits <- merge(
    df_hits,
    gene_map,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE,
    allow.cartesian = TRUE
  )

  # Same historical grouping rule:
  #   genes_entrez = paste(unique(gene_id), collapse=";")
  #   genes_symbol = paste(unique(SYMBOL), collapse=";")
  df_annot <- df_hits[
    ,
    .(
      genes_entrez = paste(
        unique(gene_id),
        collapse = ";"
      ),
      genes_symbol = paste(
        unique(SYMBOL),
        collapse = ";"
      )
    ),
    by = region_id
  ]

  # The gene filter is simply an INNER join to the overlap annotation table:
  # strict DMRs without a gene overlap are dropped.
  filtered <- merge(
    dmr_dt,
    df_annot,
    by = "region_id",
    all = FALSE,
    sort = FALSE
  )

  list(
    annotation = df_annot,
    filtered = filtered
  )
}

split_gene_symbols <- function(x) {

  x <- x[
    !is.na(x) &
      nzchar(x)
  ]

  if (!length(x)) {
    return(character())
  }

  z <- unlist(
    strsplit(
      x,
      ";",
      fixed = TRUE
    ),
    use.names = FALSE
  )

  unique(
    z[
      !is.na(z) &
        nzchar(z) &
        z != "NA"
    ]
  )
}

# -------------------------------- load inputs --------------------------------

for (f in c(
  weighted_file,
  unweighted_file,
  region_file_path
)) {
  if (!file.exists(f)) {
    stop(
      "Required input file does not exist: ",
      f
    )
  }
}

weighted <- fread(
  weighted_file,
  showProgress = FALSE
)

unweighted <- fread(
  unweighted_file,
  showProgress = FALSE
)

region_file <- fread(
  region_file_path,
  showProgress = FALSE
)

# The original region file did not explicitly store region_id; the old GAM
# runner used the ROW NUMBER of region_file as region_id. Reconstruct that
# identifier if needed.
if (!"region_id" %in% names(region_file)) {
  region_file[
    ,
    region_id := .I
  ]
}

region_file[
  ,
  region_id := as.integer(region_id)
]

weighted[
  ,
  region_id := as.integer(region_id)
]

unweighted[
  ,
  region_id := as.integer(region_id)
]

required_region_cols <- c(
  "region_id",
  "chr",
  "region_start",
  "region_end"
)

missing_region_cols <- setdiff(
  required_region_cols,
  names(region_file)
)

if (length(missing_region_cols)) {
  stop(
    "Region file is missing required columns: ",
    paste(
      missing_region_cols,
      collapse = ", "
    )
  )
}

# ----------------------- verify inherited strict DMRs ------------------------

check_dmr_flag(
  weighted,
  "DMR_weighted",
  "Read-depth-weighted GAM"
)

check_dmr_flag(
  unweighted,
  "DMR_unweighted",
  "Original unweighted GAM"
)

# ------------------------ prepare DMR genomic ranges -------------------------

weighted_dmr <- prepare_dmr_coordinates(
  weighted,
  "DMR_weighted",
  region_file,
  "weighted"
)

unweighted_dmr <- prepare_dmr_coordinates(
  unweighted,
  "DMR_unweighted",
  region_file,
  "unweighted"
)

cat(
  "\nBefore gene-overlap filtering:\n",
  "  unweighted strict DMRs = ",
  nrow(unweighted_dmr),
  "\n",
  "  weighted strict DMRs   = ",
  nrow(weighted_dmr),
  "\n",
  sep = ""
)

# ------------------------------- hg19 genes ----------------------------------

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Same gene definition used in the old code:
gr_genes <- GenomicFeatures::genes(
  txdb
)

# Make the Entrez IDs explicit and robust to package-version differences.
if (
  !"gene_id" %in% names(
    mcols(gr_genes)
  )
) {
  mcols(gr_genes)$gene_id <- names(
    gr_genes
  )
}

# ------------------------------ annotate/filter ------------------------------

weighted_gene <- annotate_gene_overlaps(
  weighted_dmr,
  gr_genes,
  "weighted"
)

unweighted_gene <- annotate_gene_overlaps(
  unweighted_dmr,
  gr_genes,
  "unweighted"
)

weighted_gene_dmr <- weighted_gene$filtered
unweighted_gene_dmr <- unweighted_gene$filtered

cat(
  "\nAfter requiring overlap with >=1 annotated gene:\n",
  "  unweighted gene-overlapping strict DMRs = ",
  nrow(unweighted_gene_dmr),
  "\n",
  "  weighted gene-overlapping strict DMRs   = ",
  nrow(weighted_gene_dmr),
  "\n",
  sep = ""
)

# ------------------------------- save outputs --------------------------------

# Historical-style region -> gene annotation tables.
fwrite(
  unweighted_gene$annotation,
  file.path(
    output_dir,
    "unweighted_DMR_strict_genes.csv"
  ),
  sep = ",",
  quote = TRUE,
  na = "NA"
)

fwrite(
  weighted_gene$annotation,
  file.path(
    output_dir,
    "weighted_DMR_strict_genes.csv"
  ),
  sep = ",",
  quote = TRUE,
  na = "NA"
)

# Full strict-DMR result rows, after the gene-overlap filter.
fwrite(
  unweighted_gene_dmr,
  file.path(
    output_dir,
    "unweighted_DMR_strict_gene_filtered.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  weighted_gene_dmr,
  file.path(
    output_dir,
    "weighted_DMR_strict_gene_filtered.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ---------------------- compare gene-filtered DMR regions --------------------

u_regions <- unique(
  unweighted_gene_dmr$region_id
)

w_regions <- unique(
  weighted_gene_dmr$region_id
)

all_regions <- sort(
  union(
    u_regions,
    w_regions
  )
)

region_comparison <- data.table(
  region_id = all_regions
)

region_comparison[
  ,
  unweighted_gene_DMR :=
    region_id %in% u_regions
]

region_comparison[
  ,
  weighted_gene_DMR :=
    region_id %in% w_regions
]

region_comparison[
  ,
  DMR_status :=
    fifelse(
      unweighted_gene_DMR &
        weighted_gene_DMR,
      "Both",
      fifelse(
        unweighted_gene_DMR,
        "Unweighted only",
        "Weighted only"
      )
    )
]

u_annot <- copy(
  unweighted_gene$annotation
)

w_annot <- copy(
  weighted_gene$annotation
)

if (nrow(u_annot)) {
  setnames(
    u_annot,
    c(
      "genes_entrez",
      "genes_symbol"
    ),
    c(
      "genes_entrez_unweighted",
      "genes_symbol_unweighted"
    )
  )

  region_comparison <- merge(
    region_comparison,
    u_annot,
    by = "region_id",
    all.x = TRUE,
    sort = FALSE
  )
}

if (nrow(w_annot)) {
  setnames(
    w_annot,
    c(
      "genes_entrez",
      "genes_symbol"
    ),
    c(
      "genes_entrez_weighted",
      "genes_symbol_weighted"
    )
  )

  region_comparison <- merge(
    region_comparison,
    w_annot,
    by = "region_id",
    all.x = TRUE,
    sort = FALSE
  )
}

setorder(
  region_comparison,
  region_id
)

fwrite(
  region_comparison,
  file.path(
    output_dir,
    "weighted_vs_unweighted_gene_filtered_DMR_regions.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ----------------------------- compare gene sets ------------------------------

u_genes <- split_gene_symbols(
  unweighted_gene$annotation$genes_symbol
)

w_genes <- split_gene_symbols(
  weighted_gene$annotation$genes_symbol
)

shared_genes <- intersect(
  u_genes,
  w_genes
)

union_genes <- union(
  u_genes,
  w_genes
)

shared_regions <- intersect(
  u_regions,
  w_regions
)

union_regions <- union(
  u_regions,
  w_regions
)

summary_dt <- data.table(
  metric = c(
    "Unweighted strict DMRs before gene filter",
    "Weighted strict DMRs before gene filter",
    "Unweighted strict DMRs overlapping >=1 gene",
    "Weighted strict DMRs overlapping >=1 gene",
    "Unweighted gene-overlap retention fraction",
    "Weighted gene-overlap retention fraction",
    "Shared gene-filtered DMR regions",
    "Union gene-filtered DMR regions",
    "Jaccard overlap of gene-filtered DMR regions",
    "Unique genes in unweighted gene-filtered DMRs",
    "Unique genes in weighted gene-filtered DMRs",
    "Shared gene symbols",
    "Union gene symbols",
    "Jaccard overlap of gene symbols"
  ),
  value = c(
    nrow(unweighted_dmr),
    nrow(weighted_dmr),
    length(u_regions),
    length(w_regions),
    if (nrow(unweighted_dmr) > 0) {
      length(u_regions) /
        nrow(unweighted_dmr)
    } else {
      NA_real_
    },
    if (nrow(weighted_dmr) > 0) {
      length(w_regions) /
        nrow(weighted_dmr)
    } else {
      NA_real_
    },
    length(shared_regions),
    length(union_regions),
    if (length(union_regions) > 0) {
      length(shared_regions) /
        length(union_regions)
    } else {
      NA_real_
    },
    length(u_genes),
    length(w_genes),
    length(shared_genes),
    length(union_genes),
    if (length(union_genes) > 0) {
      length(shared_genes) /
        length(union_genes)
    } else {
      NA_real_
    }
  )
)

fwrite(
  summary_dt,
  file.path(
    output_dir,
    "gene_filtered_DMR_comparison_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

gene_comparison <- data.table(
  gene_symbol = sort(
    union_genes
  )
)

gene_comparison[
  ,
  in_unweighted :=
    gene_symbol %in% u_genes
]

gene_comparison[
  ,
  in_weighted :=
    gene_symbol %in% w_genes
]

gene_comparison[
  ,
  status :=
    fifelse(
      in_unweighted &
        in_weighted,
      "Both",
      fifelse(
        in_unweighted,
        "Unweighted only",
        "Weighted only"
      )
    )
]

fwrite(
  gene_comparison,
  file.path(
    output_dir,
    "weighted_vs_unweighted_gene_symbol_comparison.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ---------------------------- text result summary -----------------------------

region_jaccard <- if (
  length(union_regions) > 0
) {
  length(shared_regions) /
    length(union_regions)
} else {
  NA_real_
}

gene_jaccard <- if (
  length(union_genes) > 0
) {
  length(shared_genes) /
    length(union_genes)
} else {
  NA_real_
}

summary_lines <- c(
  "GENE-FILTERED READ-DEPTH-WEIGHTED GAM-DMR SENSITIVITY ANALYSIS",
  "================================================================",
  "",
  "Gene-overlap rule:",
  "A strict GAM-DMR is retained only if its hg19 genomic interval overlaps at least one gene returned by genes(TxDb.Hsapiens.UCSC.hg19.knownGene). Entrez IDs are mapped to symbols using org.Hs.eg.db, reproducing the historical gene-annotation rule.",
  "",
  paste0(
    "Unweighted strict DMRs before gene filtering: ",
    nrow(unweighted_dmr),
    "."
  ),
  paste0(
    "Weighted strict DMRs before gene filtering: ",
    nrow(weighted_dmr),
    "."
  ),
  paste0(
    "Unweighted strict DMRs retained after gene filtering: ",
    length(u_regions),
    " (",
    sprintf(
      "%.1f",
      100 *
        length(u_regions) /
        max(
          1,
          nrow(unweighted_dmr)
        )
    ),
    "%)."
  ),
  paste0(
    "Weighted strict DMRs retained after gene filtering: ",
    length(w_regions),
    " (",
    sprintf(
      "%.1f",
      100 *
        length(w_regions) /
        max(
          1,
          nrow(weighted_dmr)
        )
    ),
    "%)."
  ),
  "",
  paste0(
    "Shared gene-filtered DMR regions: ",
    length(shared_regions),
    "."
  ),
  paste0(
    "Jaccard overlap of gene-filtered DMR regions: ",
    ifelse(
      is.finite(region_jaccard),
      sprintf(
        "%.3f",
        region_jaccard
      ),
      "NA"
    ),
    "."
  ),
  paste0(
    "Unique gene symbols in unweighted gene-filtered DMRs: ",
    length(u_genes),
    "."
  ),
  paste0(
    "Unique gene symbols in weighted gene-filtered DMRs: ",
    length(w_genes),
    "."
  ),
  paste0(
    "Shared gene symbols: ",
    length(shared_genes),
    "."
  ),
  paste0(
    "Jaccard overlap of gene symbols: ",
    ifelse(
      is.finite(gene_jaccard),
      sprintf(
        "%.3f",
        gene_jaccard
      ),
      "NA"
    ),
    "."
  ),
  ""
)

writeLines(
  summary_lines,
  con = file.path(
    output_dir,
    "gene_filtered_DMR_results_summary.txt"
  )
)

cat(
  "\n",
  paste(
    summary_lines,
    collapse = "\n"
  ),
  "\n",
  sep = ""
)

cat(
  "\nOutputs written to:\n  ",
  output_dir,
  "\n",
  sep = ""
)


suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

base_dir <- path.expand(Sys.getenv("METH_BASE_DIR", "~/scratch/UQAC/meth"))
outdir <- file.path(base_dir, "results", "11_mgcv")
script11_dir <- file.path(base_dir, "scr", "11_mgcv")
results12 <- file.path(base_dir, "scr", "12_cpg_DMR", "results")
paper_dir <- file.path(base_dir, "scr", "14_paper", "results")
dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(results12, "5_eval.RData"))
m2_cpg_signals <- M2_manhattan[M2_manhattan$pval < 1e-5, ]

dmrs <- fread(file.path(outdir, "24_dmrs_STRICT.tsv"))
region_file <- fread(file.path(script11_dir, "dat", "region_file_1_chunk.csv"))
region_file[, region_id := .I]
dmrs <- merge(dmrs, region_file, by = "region_id", all.x = TRUE) %>%
  select(region_id, chr, region_start, region_end, pvals)

dmr_gene_file <- fread(file.path(outdir, "25_DMR_strict_genes.csv"))
dmr_dt <- as.data.table(dmrs)[, .(
  chr = as.integer(chr),
  start = as.integer(region_start),
  end = as.integer(region_end),
  region_id
)]
setkey(dmr_dt, chr, start, end)

cpg_dt <- as.data.table(m2_cpg_signals)[, .(CpG, pval, coef, pval_old)]
chr_part <- tstrsplit(cpg_dt$CpG, ":", fixed = TRUE)
pos_part <- tstrsplit(chr_part[[2]], "-", fixed = TRUE)
chr_num <- suppressWarnings(as.integer(chr_part[[1]]))
pos1 <- suppressWarnings(as.integer(pos_part[[1]]))
pos2 <- suppressWarnings(as.integer(ifelse(pos_part[[2]] %in% c("NA", ""), NA, pos_part[[2]])))

cpg_iv <- data.table(
  idx = seq_len(nrow(cpg_dt)),
  chr = chr_num,
  start = ifelse(is.na(pos2), pos1, pmin(pos1, pos2)),
  end = ifelse(is.na(pos2), pos1, pmax(pos1, pos2))
)
cpg_iv <- cpg_iv[!is.na(chr) & !is.na(start) & !is.na(end)]
setkey(cpg_iv, chr, start, end)

hits <- foverlaps(cpg_iv, dmr_dt, nomatch = 0L)
flag <- integer(nrow(cpg_dt))
flag[unique(hits$idx)] <- 1L
m2_cpg_signals$overlap_DMR <- flag

hits2 <- foverlaps(dmr_dt, cpg_iv, nomatch = 0L)
dmrs$overlap_CpG <- as.integer(dmrs$region_id %in% unique(hits2$region_id))

dmr_gene_file <- merge(
  dmr_gene_file,
  dmrs[, .(region_id, chr, region_start, region_end, pvals, overlap_CpG)],
  by = "region_id",
  all.x = TRUE
)

dmr_m2 <- dmrs
m1_env <- new.env(parent = emptyenv())
load(file.path(results12, "6_overlap_results.RData"), envir = m1_env)
if (!exists("dmrs", envir = m1_env, inherits = FALSE)) stop("dmrs is missing from 6_overlap_results.RData")
dmr_m1 <- get("dmrs", envir = m1_env)

save(m2_cpg_signals, dmrs, cpg_iv, dmr_gene_file, dmr_dt,
     file = file.path(paper_dir, "4_overlap_M2.RData"))
save(dmr_m1, dmr_m2, dmr_gene_file,
     file = file.path(paper_dir, "4_DMR_M1_M2_overlap.RData"))

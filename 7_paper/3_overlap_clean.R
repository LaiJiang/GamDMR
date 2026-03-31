suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

base_dir <- Sys.getenv("METH_BASE_DIR", unset = getwd())
outdir <- file.path(base_dir, "results", "11_mgcv")
script11_dir <- file.path(base_dir, "scr", "11_mgcv")
results12 <- file.path(base_dir, "scr", "12_cpg_DMR", "results")

load(file.path(results12, "5_eval.RData"))

m1_cpg_signals <- M1_manhattan[M1_manhattan$pval < 1e-5, ]
cpgs_int <- unique(M1_manhattan$CpG[M1_manhattan$pval < 1e-5])

source(file.path(base_dir, "scr", "8_rerun", "0_func_enrichment_updates.R"))
m1_genes <- gene_map

dmrs <- fread(file.path(outdir, "24_dmrs_STRICT.tsv"))
region_file <- fread(file.path(script11_dir, "dat", "region_file_1_chunk.csv"))
region_file$region_id <- seq_len(nrow(region_file))
dmrs <- merge(dmrs, region_file, by = "region_id", all.x = TRUE) %>%
  select(region_id, chr, region_start, region_end, pvals)

dmr_gene_file <- fread(file.path(base_dir, "results", "11_mgcv", "25_DMR_strict_genes.csv"))

dmr_dt <- as.data.table(dmrs)[, .(
  chr = as.integer(chr),
  start = as.integer(region_start),
  end = as.integer(region_end),
  region_id
)]
setkey(dmr_dt, chr, start, end)

cpg_dt <- as.data.table(m1_cpg_signals)[, .(CpG, pval, coef, pval_old)]
chr_part <- tstrsplit(cpg_dt$CpG, ":", fixed = TRUE)
pos_part <- tstrsplit(chr_part[[2]], "-", fixed = TRUE)

chr_num <- suppressWarnings(as.integer(chr_part[[1]]))
pos1 <- suppressWarnings(as.integer(pos_part[[1]]))
pos2 <- suppressWarnings(as.integer(ifelse(pos_part[[2]] %in% c("NA", ""), NA, pos_part[[2]])))

cpg_iv <- data.table(
  idx = seq_len(nrow(cpg_dt)),
  chr = chr_num,
  start = ifelse(is.na(pos2), pos1, pmin(pos1, pos2)),
  end = ifelse(is.na(pos2), pos1, pmin(pos1, pos2))
)
cpg_iv <- cpg_iv[!is.na(chr) & !is.na(start) & !is.na(end)]
setkey(cpg_iv, chr, start, end)

hits <- foverlaps(cpg_iv, dmr_dt, nomatch = 0L)
overlap_idx <- unique(hits$idx)
overlap_flag <- integer(nrow(cpg_dt))
overlap_flag[overlap_idx] <- 1L
m1_cpg_signals$overlap_DMR <- overlap_flag

hits2 <- foverlaps(dmr_dt, cpg_iv, nomatch = 0L)
overlap_dmr_ids <- unique(hits2$region_id)
dmrs$overlap_CpG <- as.integer(dmrs$region_id %in% overlap_dmr_ids)

dmr_gene_file <- merge(
  dmr_gene_file,
  dmrs[, .(region_id, chr, region_start, region_end, pvals, overlap_CpG)],
  by = "region_id",
  all.x = TRUE
)

save(m1_cpg_signals, dmrs, cpg_iv, dmr_gene_file, dmr_dt, file = file.path(results12, "6_overlap_results.RData"))

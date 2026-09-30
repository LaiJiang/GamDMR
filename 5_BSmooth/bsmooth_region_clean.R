#!/usr/bin/env Rscript

mem_rss_mb <- function() {
  p <- readLines("/proc/self/status")
  as.numeric(sub(".*:\\s+([0-9]+).*", "\\1", grep("^VmRSS:", p, value = TRUE))) / 1024
}

suppressPackageStartupMessages({
  library(bsseq)
  library(BiocParallel)
  library(data.table)
  library(dplyr)
  library(stringr)
})

data.table::setDTthreads(1L)
BiocParallel::register(BiocParallel::SerialParam())
options(mc.cores = 1L)

min_cpgs <- 10L
n_jobs <- 300L
PATH_wk <- path.expand(Sys.getenv("METH_BASE_DIR", "~/scratch/UQAC/meth"))
PATH_scr11 <- file.path(PATH_wk, "scr", "11_mgcv")
output_dir <- file.path(PATH_wk, "results", "13_bsmooth")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) stop("Usage: Rscript bsmooth_region_clean.R <job_id>")
job_id <- as.integer(args[1])
if (is.na(job_id)) stop("job_id must be an integer.")

region_file <- fread(file.path(PATH_scr11, "dat", "region_file_1_chunk.csv"))
region_file[, region_id := .I]
region_index_list <- split(
  seq_len(nrow(region_file)),
  cut(seq_len(nrow(region_file)), breaks = n_jobs, labels = FALSE)
)
if (job_id < 1L || job_id > length(region_index_list)) quit(save = "no", status = 0L)
region_ids <- region_index_list[[job_id]]
if (!length(region_ids)) quit(save = "no", status = 0L)

env <- new.env(parent = emptyenv())
load(file.path(PATH_scr11, "dat", "bsmooth_only_pheno_bmi.RData"), envir = env)
if (!exists("pheno_file", envir = env, inherits = FALSE)) stop("pheno_file is missing.")
pheno_file <- as_tibble(get("pheno_file", envir = env)) %>%
  transmute(ID = as.character(ID), AA_only = AA_only)
rm(env)

build_bsseq_matrices <- function(dt_region) {
  meth_cols <- grep("_meth$", names(dt_region), value = TRUE)
  if (!length(meth_cols)) stop("No methylation columns found.")

  sample_ids <- str_replace(meth_cols, "_meth$", "")
  cov_cols <- paste0(sample_ids, "_tot")
  have_cov <- cov_cols %in% names(dt_region)

  dt_region[, (meth_cols) := lapply(.SD, as.numeric), .SDcols = meth_cols]
  if (any(have_cov)) {
    dt_region[, (cov_cols[have_cov]) := lapply(.SD, as.numeric), .SDcols = cov_cols[have_cov]]
  }
  setorder(dt_region, start)

  meth <- as.matrix(dt_region[, ..meth_cols])
  if (all(have_cov)) {
    cov <- as.matrix(dt_region[, ..cov_cols])
  } else {
    cov <- matrix(20, nrow = nrow(meth), ncol = ncol(meth))
  }

  meth[is.na(meth)] <- 0
  cov[is.na(cov)] <- 0
  meth <- pmin(pmax(meth, 0), 1)
  m_counts <- round(meth * cov)

  colnames(m_counts) <- sample_ids
  colnames(cov) <- sample_ids
  list(M = m_counts, Cov = cov, sample_names = sample_ids, pos = dt_region$start)
}

tstat_out_path <- file.path(output_dir, sprintf("bsmooth_tstat_job_%d.rds", job_id))
dmr_out_path <- file.path(output_dir, sprintf("dmr_job_%d.tsv.gz", job_id))
summary_path <- file.path(output_dir, sprintf("summary_job_%d.tsv", job_id))

regions_this_job <- copy(region_file[region_ids])
chunks <- split(regions_this_job, regions_this_job$data_chunk_id)
summary_rows <- list()
dmr_rows <- list()
all_tstats <- list()

for (chunk_name in names(chunks)) {
  chunk_regions <- as.data.table(chunks[[chunk_name]])
  chunk_id <- as.integer(chunk_name)
  chunk_file <- file.path(PATH_wk, "data", "meth_split", sprintf("chunk_%04d.csv", chunk_id))
  if (!file.exists(chunk_file)) next

  header <- fread(chunk_file, nrows = 0L, showProgress = FALSE)
  all_cols <- names(header)
  file_sample_ids <- gsub("_meth$|_tot$", "", setdiff(all_cols, c("chr", "start", "end")))
  keep_ids <- intersect(file_sample_ids, pheno_file$ID)
  if (!length(keep_ids)) next

  selected <- c("chr", "start", "end", paste0(keep_ids, "_meth"), paste0(keep_ids, "_tot"))
  selected <- intersect(selected, all_cols)
  chunk_dt <- fread(chunk_file, select = selected, showProgress = FALSE)
  chunk_dt[, chr := as.character(chr)]

  for (j in seq_len(nrow(chunk_regions))) {
    info <- chunk_regions[j]
    region_dt <- chunk_dt[
      chr == as.character(info$chr) &
        start >= info$region_start &
        start <= info$region_end
    ]
    if (nrow(region_dt) < min_cpgs) next

    mats <- tryCatch(build_bsseq_matrices(region_dt), error = function(e) NULL)
    if (is.null(mats)) next

    ph <- pheno_file %>% filter(ID %in% mats$sample_names, !is.na(AA_only))
    if (nrow(ph) < 2L) next
    idx <- match(ph$ID, mats$sample_names)
    idx <- idx[!is.na(idx)]
    if (length(idx) < 2L) next

    m2 <- mats$M[, idx, drop = FALSE]
    cov2 <- mats$Cov[, idx, drop = FALSE]
    sample_ids2 <- mats$sample_names[idx]
    grp <- as.integer(ph$AA_only)
    n_case <- sum(grp == 1L, na.rm = TRUE)
    n_ctrl <- sum(grp == 0L, na.rm = TRUE)
    if (n_case < 1L || n_ctrl < 1L) next

    bs <- BSseq(
      M = m2,
      Cov = cov2,
      chr = rep(as.character(info$chr), nrow(m2)),
      pos = mats$pos,
      sampleNames = sample_ids2
    )

    bs_sm <- tryCatch(
      BSmooth(bs, ns = 25, h = 400, maxGap = 5e7, verbose = FALSE, BPPARAM = SerialParam()),
      error = function(e) NULL
    )
    if (is.null(bs_sm)) next

    tstat <- tryCatch(
      BSmooth.tstat(
        bs_sm,
        group1 = which(grp == 1L),
        group2 = which(grp == 0L),
        estimate.var = "group2",
        local.correct = TRUE,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(tstat)) next

    all_tstats[[as.character(info$region_id)]] <- tstat
    st <- bsseq::getStats(tstat)
    stat_vec <- if (is.matrix(st)) {
      preferred <- intersect(c("tstat.corrected", "tstat", "coef", "stat"), colnames(st))
      if (!length(preferred)) NULL else st[, preferred[1L], drop = TRUE]
    } else {
      st
    }
    if (is.null(stat_vec) || !is.numeric(stat_vec) || !any(is.finite(stat_vec))) next

    dmrs <- tryCatch(
      dmrFinder(tstat, cutoff = c(-4.417, 4.417)),
      error = function(e) NULL
    )
    n_dmrs <- if (is.null(dmrs) || !nrow(dmrs)) 0L else nrow(dmrs)

    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      region_id = as.integer(info$region_id),
      data_chunk_id = as.integer(info$data_chunk_id),
      chr = as.character(info$chr),
      region_start = as.integer(info$region_start),
      region_end = as.integer(info$region_end),
      n_cpg = nrow(m2),
      n_samples = ncol(m2),
      n_case = n_case,
      n_ctrl = n_ctrl,
      n_dmrs = n_dmrs
    )

    if (n_dmrs > 0L) {
      dmrs$region_id <- as.integer(info$region_id)
      dmrs$data_chunk_id <- as.integer(info$data_chunk_id)
      dmrs$chr_region <- as.character(info$chr)
      dmrs$reg_start <- as.integer(info$region_start)
      dmrs$reg_end <- as.integer(info$region_end)
      dmrs$n_cpgs <- nrow(m2)
      dmrs$n_samples <- ncol(m2)
      if ("areaStat" %in% names(dmrs)) {
        dmrs$score <- abs(dmrs$areaStat)
        dmrs$direction <- sign(dmrs$areaStat)
      } else if ("maxStat" %in% names(dmrs)) {
        dmrs$score <- abs(dmrs$maxStat)
        dmrs$direction <- sign(dmrs$maxStat)
      } else {
        dmrs$score <- NA_real_
        dmrs$direction <- NA_integer_
      }
      dmrs$midpoint <- floor((dmrs$start + dmrs$end) / 2)
      if (!("width" %in% names(dmrs)) && all(c("start", "end") %in% names(dmrs))) {
        dmrs$width <- dmrs$end - dmrs$start + 1L
      }
      dmr_rows[[length(dmr_rows) + 1L]] <- as.data.table(dmrs)
    }
  }
}

if (length(summary_rows)) fwrite(rbindlist(summary_rows, fill = TRUE), summary_path, sep = "\t")
if (length(dmr_rows)) fwrite(rbindlist(dmr_rows, fill = TRUE), dmr_out_path, sep = "\t")
saveRDS(all_tstats, tstat_out_path)

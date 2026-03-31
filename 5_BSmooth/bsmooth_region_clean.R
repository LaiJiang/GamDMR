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

data.table::setDTthreads(1)
BiocParallel::register(BiocParallel::SerialParam())
options(mc.cores = 1)

min_cpgs <- 10
n_jobs <- 300

base_dir <- Sys.getenv("METH_BASE_DIR", unset = getwd())
script_dir <- file.path(base_dir, "11_mgcv")
output_dir <- file.path(base_dir, "results", "13_bsmooth")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript bsmooth_region_clean.R <job_id>")
job_id <- as.integer(args[1])
if (is.na(job_id)) stop("job_id must be an integer.")

region_file <- fread(file.path(script_dir, "dat", "region_file_1_chunk.csv"))
region_index_list <- split(
  seq_len(nrow(region_file)),
  cut(seq_len(nrow(region_file)), breaks = n_jobs, labels = FALSE)
)
region_ids <- region_index_list[[job_id]]
if (length(region_ids) == 0 || all(is.na(region_ids))) {
  cat("No regions assigned to this job.\n")
  quit(save = "no", status = 0)
}

load(file.path(script_dir, "dat", "bsmooth_only_pheno_bmi.RData"))
pheno_file <- as_tibble(pheno_file) %>% select(ID, AA_only)

build_bsseq_matrices <- function(dt_region) {
  meth_cols <- grep("_meth$", names(dt_region), value = TRUE)
  if (length(meth_cols) == 0) stop("No *_meth columns found.")
  sample_ids <- str_replace(meth_cols, "_meth$", "")
  cov_cols <- paste0(sample_ids, "_tot")
  have_cov <- cov_cols %in% names(dt_region)

  dt_region[, (meth_cols) := lapply(.SD, as.numeric), .SDcols = meth_cols]
  if (any(have_cov)) {
    dt_region[, (cov_cols[have_cov]) := lapply(.SD, as.numeric), .SDcols = cov_cols[have_cov]]
  }

  setorder(dt_region, start)
  meth <- as.matrix(dt_region[, ..meth_cols])
  cov <- if (all(have_cov)) {
    as.matrix(dt_region[, ..cov_cols])
  } else {
    warning("Some coverage columns are missing; using constant coverage = 20 where needed.")
    matrix(20, nrow = nrow(meth), ncol = ncol(meth))
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

if (!file.exists(summary_path)) {
  cat(
    "region_id", "data_chunk_id", "chr", "region_start", "region_end",
    "n_cpg", "n_samples", "n_case", "n_ctrl", "n_dmrs",
    sep = "\t", file = summary_path, append = FALSE
  )
  cat("\n", file = summary_path, append = TRUE)
}

regions_this_job <- copy(region_file[region_ids])
regions_this_job$region_id <- region_ids
chunks <- split(regions_this_job, regions_this_job$data_chunk_id)

all_tstats <- list()

for (chunk_name in names(chunks)) {
  chunk_regions <- chunks[[chunk_name]]
  chunk_id_num <- as.integer(chunk_name)

  chunk_fp <- file.path(base_dir, "data", "meth_split", sprintf("chunk_%04d.csv", chunk_id_num))
  if (!file.exists(chunk_fp)) {
    cat("Missing chunk:", chunk_fp, "\n")
    next
  }

  hdr <- fread(chunk_fp, nrows = 0, showProgress = FALSE)
  all_cols <- names(hdr)
  file_sample_ids <- gsub("_meth$|_tot$", "", setdiff(all_cols, c("chr", "start", "end")))
  keep_ids <- intersect(file_sample_ids, pheno_file$ID)
  if (length(keep_ids) == 0L) {
    cat("No overlapping samples in chunk", chunk_id_num, "\n")
    next
  }

  sel_cols <- c("chr", "start", "end", paste0(keep_ids, "_meth"), paste0(keep_ids, "_tot"))
  sel_cols <- intersect(sel_cols, all_cols)

  i_chunk <- fread(chunk_fp, select = sel_cols, showProgress = FALSE)
  i_chunk[, chr := as.character(chr)]

  for (row_i in seq_len(nrow(chunk_regions))) {
    i_info <- chunk_regions[row_i, ]
    i_info$chr <- as.character(i_info$chr)

    i_region_dt <- i_chunk[
      chr == i_info$chr & start >= i_info$region_start & start <= i_info$region_end
    ]

    if (nrow(i_region_dt) < min_cpgs) next

    mats <- tryCatch(build_bsseq_matrices(i_region_dt), error = function(e) NULL)
    if (is.null(mats)) next

    ph_sub <- pheno_file %>% filter(ID %in% mats$sample_names, !is.na(AA_only))
    if (nrow(ph_sub) < 2) next

    idx <- match(ph_sub$ID, mats$sample_names)
    idx <- idx[!is.na(idx)]
    if (length(idx) < 2) next

    m2 <- mats$M[, idx, drop = FALSE]
    cov2 <- mats$Cov[, idx, drop = FALSE]
    sample_ids2 <- mats$sample_names[idx]
    grp_vec <- as.integer(ph_sub$AA_only)

    n_case <- sum(grp_vec == 1, na.rm = TRUE)
    n_ctrl <- sum(grp_vec == 0, na.rm = TRUE)
    if (n_case < 1 || n_ctrl < 1) next

    bs <- BSseq(
      M = m2,
      Cov = cov2,
      chr = rep(i_info$chr, nrow(m2)),
      pos = mats$pos,
      sampleNames = sample_ids2
    )

    bs_sm <- tryCatch(
      BSmooth(bs, ns = 25, h = 400, maxGap = 5e7, verbose = FALSE, BPPARAM = SerialParam()),
      error = function(e) NULL
    )
    if (is.null(bs_sm)) {
      rm(bs)
      gc(FALSE)
      next
    }

    g1 <- which(grp_vec == 1)
    g2 <- which(grp_vec == 0)

    tstat <- tryCatch(
      BSmooth.tstat(
        bs_sm,
        group1 = g1,
        group2 = g2,
        estimate.var = "group2",
        local.correct = TRUE,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(tstat)) {
      rm(bs, bs_sm)
      gc(FALSE)
      next
    }

    all_tstats[[as.character(i_info$region_id)]] <- tstat

    st <- bsseq::getStats(tstat)
    stat_vec <- if (is.matrix(st)) {
      preferred <- intersect(c("tstat.corrected", "tstat", "coef", "stat"), colnames(st))
      if (length(preferred) == 0L) NULL else st[, preferred[1], drop = TRUE]
    } else {
      st
    }

    if (is.null(stat_vec) || !is.numeric(stat_vec) || !any(is.finite(stat_vec))) {
      rm(bs, bs_sm, tstat)
      gc(FALSE)
      next
    }

    dmrs <- tryCatch(
      dmrFinder(tstat, cutoff = c(-4.417, 4.417)),
      error = function(e) NULL
    )

    n_dmrs <- if (is.null(dmrs) || nrow(dmrs) == 0) 0L else nrow(dmrs)

    cat(
      i_info$region_id, i_info$data_chunk_id, i_info$chr,
      i_info$region_start, i_info$region_end,
      nrow(m2), ncol(m2), n_case, n_ctrl, n_dmrs,
      sep = "\t", file = summary_path, append = TRUE
    )
    cat("\n", file = summary_path, append = TRUE)

    if (n_dmrs > 0) {
      dmrs$region_id <- i_info$region_id
      dmrs$data_chunk_id <- i_info$data_chunk_id
      dmrs$chr_region <- i_info$chr
      dmrs$reg_start <- i_info$region_start
      dmrs$reg_end <- i_info$region_end
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

      fwrite(
        dmrs,
        file = dmr_out_path,
        append = file.exists(dmr_out_path),
        sep = "\t",
        col.names = !file.exists(dmr_out_path)
      )
    }

    rm(i_region_dt, mats, m2, cov2, bs, bs_sm, tstat, dmrs, stat_vec)
    gc(FALSE)
  }

  rm(i_chunk, hdr, all_cols, file_sample_ids, keep_ids, sel_cols, chunk_regions)
  gc(FALSE)
}

saveRDS(all_tstats, tstat_out_path)

cat("Done.\n")
cat("T-statistics:", tstat_out_path, "\n")
cat("DMRs:", dmr_out_path, "\n")
cat("Summary:", summary_path, "\n")

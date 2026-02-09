#!/usr/bin/env Rscript
# BSmooth on a per-chunk, per-region workflow aligned to your mgcv pipeline
# Changes applied:
#  (3) Read only needed columns from each chunk via fread(select=...)
#  (4) Load each chunk ONCE, then iterate its regions
# Also keeps prior hardening (tstat validation, NA phenotype guard, memory cleanup).

#memory tracker
mem_rss_mb <- function() {
  p <- readLines("/proc/self/status")
  as.numeric(sub(".*:\\s+([0-9]+).*", "\\1", grep("^VmRSS:", p, value = TRUE)))/1024
}


suppressPackageStartupMessages({
  library(bsseq)
  
  library(BiocParallel)

  library(data.table)  
  data.table::setDTthreads(1)

  library(dplyr)
  library(stringr)
  # library(BiocParallel)  # Uncomment on Linux/SLURM if you want parallel smoothing with BPPARAM
})


# Hard-stop any parallel backend
BiocParallel::register(BiocParallel::SerialParam())
options(mc.cores = 1)

#
## ---- Config (mirror your mgcv script) ---------------------------------
min_cpgs   <- 10
N_jobs     <- 300
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
#PATH_wk <- "~/scratch/UQAC/meth/"
PATH_scr11 <- file.path(PATH_wk, "scr/11_mgcv/")
PATH_output <- file.path(PATH_wk, "results/13_bsmooth/test/")
if (!dir.exists(PATH_output)) dir.create(PATH_output, recursive = TRUE)

## ---- Args --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("No SLURM_ARRAY_TASK_ID provided.")
i_job_id <- as.integer(args[1])
cat("BSmooth job id:", i_job_id, "\n")

i_job_id = 1
## ---- Load region index -------------------------------------------------
region_file <- fread(file.path(PATH_scr11, "dat/region_file_1_chunk.csv"))
vec_list <- split(seq_len(nrow(region_file)),
                  cut(seq_len(nrow(region_file)), breaks = N_jobs, labels = FALSE))
i_region_id_vector <- vec_list[[i_job_id]]

## ---- Load phenotype (AA_only used for groups) --------------------------
load(file = file.path(PATH_scr11, "dat/bsmooth_only_pheno_bmi.RData"), verbose = TRUE)
#load(file = file.path(PATH_scr11, "dat/18_pheno_BMI.RData"), verbose = TRUE)
# expect 'pheno_file' with columns: ID, AA_only, plus covariates
pheno_file <- as_tibble(pheno_file) %>% select(ID, AA_only)

## ---- Helpers -----------------------------------------------------------
# Build M and Cov matrices from a region data.table with *_meth and *_tot columns
.build_bsseq_matrices <- function(dt_region) {
  meth_cols <- grep("_meth$", names(dt_region), value = TRUE)
  if (length(meth_cols) == 0) stop("No *_meth columns found in region chunk.")
  samp_ids <- str_replace(meth_cols, "_meth$", "")

  cov_cols <- paste0(samp_ids, "_tot")
  have_cov <- cov_cols %in% names(dt_region)

  # Coerce numeric
  dt_region[, (meth_cols) := lapply(.SD, as.numeric), .SDcols = meth_cols]
  if (any(have_cov)) {
    dt_region[, (cov_cols[have_cov]) := lapply(.SD, as.numeric), .SDcols = cov_cols[have_cov]]
  }

  # Order CpGs by genomic start
  setorder(dt_region, start)

  # Assemble matrices [rows = CpGs, cols = samples]
  Meth <- as.matrix(dt_region[, ..meth_cols])
  if (all(have_cov)) {
    Cov <- as.matrix(dt_region[, ..cov_cols])
  } else {
    warning("Coverage columns not found for all samples; using constant coverage = 20.")
    Cov <- matrix(20, nrow = nrow(Meth), ncol = ncol(Meth))
  }

  # Ensure bounds & integer counts
  Meth[is.na(Meth)] <- 0
  Cov[is.na(Cov)]   <- 0
  Meth <- pmin(pmax(Meth, 0), 1)
  M_counts <- round(Meth * Cov)

  # IMPORTANT: set assay colnames to the bare sample IDs
  colnames(M_counts) <- samp_ids
  colnames(Cov)      <- samp_ids

  # Return sorted positions as well (belt-and-suspenders)
  pos_sorted <- dt_region$start

  list(M = M_counts, Cov = Cov, sample_names = samp_ids, pos = pos_sorted)
}

## ---- Outputs -----------------------------------------------------------
tstat_out_path <- file.path(PATH_output, sprintf("bsmooth_tstat_job_%d.rds", i_job_id))
dmr_out_path   <- file.path(PATH_output, sprintf("dmr_job_%d.tsv.gz", i_job_id))

# BEFORE the loops — define path and ensure dir exists
log_summary <- file.path(PATH_output, sprintf("summary_job_%d.tsv", i_job_id))
dir.create(dirname(log_summary), recursive = TRUE, showWarnings = FALSE)
cat("Summary path:", log_summary, "\n")

# Write header exactly once if file doesn't exist
if (!file.exists(log_summary)) {
  cat(
    "region_id", "data_chunk_id", "chr", "region_start", "region_end",
    "n_cpg", "n_samples", "n_case", "n_ctrl", "n_dmrs",
    sep = "\t", file = log_summary, append = FALSE
  )
  cat("\n", file = log_summary, append = TRUE)
}

cat("** per-chunk mode enabled (select columns) **\n")

## ---- Group regions by chunk (data_chunk_id) ----------------------------
# Work only on the subset of regions for this task, grouped by chunk id
regions_this_job <- region_file[i_region_id_vector]
#!!!
regions_this_job$region_id <- i_region_id_vector


chunks <- split(regions_this_job, regions_this_job$data_chunk_id)  # named by chunk id

## ---- Main loop: per-chunk, then per-region -----------------------------
#!!!
print("names of chunks:")
print(names(chunks))

#!!!
regions_this_job[regions_this_job$region_start == 10469,]

for (chunk_name in  names(chunks) ) {
  chunk_regions <- chunks[[chunk_name]]
  chunk_id_num <- as.integer(chunk_name)

#!!!
  chunk_regions <- chunks[["1"]]
  chunk_id_num <- 1

  # Build chunk file path
  #chunk_fp <- file.path(PATH_wk, "data/meth_split", sprintf("chunk_%04d.csv", chunk_id_num))
    chunk_fp <- file.path(PATH_wk, "dat/", sprintf("chunk_%04d.csv", chunk_id_num))


  if (!file.exists(chunk_fp)) {
    cat("Missing chunk:", chunk_fp, " — skipping all regions in this chunk\n")
    next
  }

  # ---- (3) Read only needed columns from this chunk --------------------
  hdr <- data.table::fread(chunk_fp, nrows = 0, showProgress = FALSE)
  all_cols <- names(hdr)
  file_sample_ids <- gsub("_meth$|_tot$", "", setdiff(all_cols, c("chr","start","end")))
  keep_ids <- intersect(file_sample_ids, pheno_file$ID)
  if (length(keep_ids) == 0L) {
    cat("No overlapping samples in chunk", chunk_id_num, "with phenotype — skipping\n")
    next
  }
  sel_cols <- c("chr","start","end", paste0(keep_ids, "_meth"), paste0(keep_ids, "_tot"))
  sel_cols <- intersect(sel_cols, all_cols)

  cat(sprintf("Reading chunk %s with %d selected columns (%d samples)\n",
              chunk_id_num, length(sel_cols), length(keep_ids)))
  i_chunk <- data.table::fread(chunk_fp, select = sel_cols, showProgress = FALSE)

  i_chunk[, chr := as.character(chr)]
  print ("i_chunk chrs:")
  print(unique(i_chunk$chr))

  # ---- Iterate regions within this chunk --------------------------------
  for (row_i in seq_len(nrow(chunk_regions))) {

    row_i <- which(chunk_regions$region_start == 10469)

    i_info <- chunk_regions[row_i, ]
    
    i_info$chr <- as.character(i_info$chr)
    print ("i_info chrs:")
    print(unique(i_info$chr))

    print ("i_info:")
    print(dim(i_info))
    # Subset to region
    i_region_dt <- i_chunk %>%
      filter(chr == i_info$chr,
             start >= i_info$region_start,
             start <= i_info$region_end) %>%
      as.data.table()

    print ("i_region_dt:")
    print(dim(i_region_dt))

    if (nrow(i_region_dt) < min_cpgs) {
      cat("Region", i_info$region_id, "insufficient CpGs (", nrow(i_region_dt), ") — skip\n")
      next
    }

    # Build bsseq matrices
    mats <- .build_bsseq_matrices(i_region_dt)
    M   <- mats$M
    Cov <- mats$Cov
    samp_ids <- mats$sample_names

    # Intersect with phenotype and build groups (drop NA phenotypes)
    ph_sub <- pheno_file %>% filter(ID %in% samp_ids, !is.na(AA_only))
    if (nrow(ph_sub) < 2) {
      cat("Region", i_info$region_id, "— too few samples with non-NA phenotype — skip\n")
      next
    }
    if (!any(samp_ids %in% ph_sub$ID)) {
      cat("Region", i_info$region_id, "— no overlap between samples and phenotype — skip\n")
      next
    }

    # Reorder to phenotype order
    idx <- match(ph_sub$ID, samp_ids)
    idx <- idx[!is.na(idx)]

    M2         <- M[, idx, drop = FALSE]
    Cov2       <- Cov[, idx, drop = FALSE]
    samp_ids2  <- samp_ids[idx]
    grp_vec    <- as.integer(ph_sub$AA_only)

    # Sanity checks
    stopifnot(identical(colnames(M2), samp_ids2))
    stopifnot(identical(colnames(Cov2), samp_ids2))

    n_case <- sum(grp_vec == 1, na.rm = TRUE)
    n_ctrl <- sum(grp_vec == 0, na.rm = TRUE)
    if (n_case < 1 || n_ctrl < 1) {
      cat("Region", i_info$region_id, "— need at least 1 case and 1 control — skip\n")
      next
    }

    cat(sprintf("Region %d: %s:%d-%d | CpGs=%d | samples=%d\n",
                i_info$region_id, i_info$chr, i_info$region_start, i_info$region_end,
                nrow(M2), ncol(M2)))

    # Build BSseq object for this region


    # Build BSseq object for this region
    cat("Building BSseq object...\n")
    print(dim(M2))
    print(dim(Cov2))
    cat("RSS before BSseq build (MB):", mem_rss_mb(), "\n")

    BS <- BSseq(
      M = M2,
      Cov = Cov2,
      chr = rep(i_info$chr, nrow(M2)),
      pos = mats$pos,
      sampleNames = samp_ids2
    )

    cat("RSS after BSseq build (MB):", mem_rss_mb(), "\n")

    # Smoothing (tunable)
    cat("RSS before BSmooth (MB):", mem_rss_mb(), "\n")
    BS_sm <- BSmooth(BS, ns = 25, h = 400, maxGap = 5e7, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
    cat("RSS after BSmooth (MB):", mem_rss_mb(), "\n")

    # Group indices
    g1 <- which(grp_vec == 1)
    g2 <- which(grp_vec == 0)

    #sanity checks on group variances
    table(grp_vec)                      # check sizes
    anyNA(grp_vec); length(unique(grp_vec))
    var_g1 <- apply(M2[, grp_vec==1, drop=FALSE], 1, var, na.rm=TRUE)
    var_g2 <- apply(M2[, grp_vec==0, drop=FALSE], 1, var, na.rm=TRUE)
    summary(var_g1); summary(var_g2)



    # t-statistics
    cat("RSS before BSmooth.tstat (MB):", mem_rss_mb(), "\n")
    tstat <- BSmooth.tstat(
      BS_sm,
      group1 = g1,
      group2 = g2,
      estimate.var = "group2",
      local.correct = TRUE,
      verbose = FALSE
    )
    cat("RSS after BSmooth.tstat (MB):", mem_rss_mb(), "\n")

    print("tstat computed.")
    print(tstat)

    # --- validate t-stat object before dmrFinder ---
    ## --- validate t-stat object before dmrFinder ---
    # Extract the statistic we want to threshold on
    st <- bsseq::getStats(tstat)

    # Inspect once (optional)
    if (is.matrix(st)) {
      cat("getStats returned matrix with columns:", paste(colnames(st), collapse=", "), "\n")
    }

    # Prefer a corrected t-stat if present; otherwise fall back sanely
    pick_col <- NULL
    if (is.matrix(st)) {
      # ordered preference; adjust to match your bsseq version's colnames
      prefs <- c("tstat.corrected", "tstat", "coef", "stat")
      pick_col <- intersect(prefs, colnames(st))
      if (length(pick_col) == 0L) {
        stop("Unknown getStats columns: ", paste(colnames(st), collapse=", "))
      }
      stat_vec <- st[, pick_col[1], drop = TRUE]
    } else {
      stat_vec <- st
    }

    # Basic sanity checks
    if (!is.numeric(stat_vec) || length(stat_vec) == 0L) {
      cat("Region", i_info$region_id, "— invalid t-stat; skipping\n"); next
    }
    if (!any(is.finite(stat_vec))) {
      cat("Region", i_info$region_id, "— all t-stat are NA/Inf; skipping\n"); next
    }

    # Length must match number of loci
    gr_len <- length(GenomicRanges::granges(tstat))
    if (length(stat_vec) != gr_len) {
      cat("Region", i_info$region_id, "— length(stat) != length(granges); skipping\n"); next
    }

    # Diagnostics (concise)
    qs <- quantile(abs(stat_vec), c(.90,.95,.99,.995), na.rm=TRUE)
    cat("t-stat column used:", if (is.null(pick_col)) "vector" else pick_col[1],
        " | length:", length(stat_vec),
        " | finite:", sum(is.finite(stat_vec)),
        " | |t| quantiles:", paste(signif(qs,4), collapse=", "), "\n")


    # DMRs (cutoff = c(-4.417, 4.417) ≈ two-sided p ~ 1e-5)
    dmrs <- try({
      dmrFinder(tstat, cutoff = c(-3, 3))
    }, silent = TRUE)



    if (inherits(dmrs, "try-error")) {
      cat("Region", i_info$region_id, "— dmrFinder error:", as.character(dmrs), "\n")
      n_dmrs <- 0L
    } else {
      n_dmrs <- if (is.null(dmrs) || nrow(dmrs) == 0) 0L else nrow(dmrs)
    }

    # Per-region summary line
    cat(
    i_info$region_id,
    i_info$data_chunk_id,
    i_info$chr,
    i_info$region_start,
    i_info$region_end,
    nrow(M2),                 # n_cpg
    ncol(M2),                 # n_samples
    n_case,
    n_ctrl,
    n_dmrs,
    sep = "\t", file = log_summary, append = TRUE
    )
    cat("\n", file = log_summary, append = TRUE)


    # Annotate and write out DMRs (append)
    if (n_dmrs > 0) {
      dmrs$region_id     <- i_info$region_id
      dmrs$data_chunk_id <- i_info$data_chunk_id
      dmrs$chr_region    <- i_info$chr
      dmrs$reg_start     <- i_info$region_start
      dmrs$reg_end       <- i_info$region_end
      dmrs$n_cpgs        <- nrow(M2)
      dmrs$n_samples     <- ncol(M2)

      # Ranking/plotting helpers
      if ("areaStat" %in% names(dmrs)) {
        dmrs$score <- abs(dmrs$areaStat)
      } else if ("maxStat" %in% names(dmrs)) {
        dmrs$score <- abs(dmrs$maxStat)
      } else {
        dmrs$score <- NA_real_
      }

      dmrs$direction <- if ("areaStat" %in% names(dmrs)) {
        sign(dmrs$areaStat)
      } else if ("meanDiff" %in% names(dmrs)) {
        sign(dmrs$meanDiff)
      } else if ("maxStat" %in% names(dmrs)) {
        sign(dmrs$maxStat)
      } else {
        NA_integer_
      }

      dmrs$midpoint <- floor((dmrs$start + dmrs$end) / 2)
      if (!("width" %in% names(dmrs)) && all(c("start","end") %in% names(dmrs))) {
        dmrs$width <- (dmrs$end - dmrs$start + 1L)
      }

      data.table::fwrite(
        dmrs,
        file = dmr_out_path,      # e.g. "dmr_job_1.tsv"
        append = file.exists(dmr_out_path),
        sep = "\t",
        col.names = !file.exists(dmr_out_path)
      )

    }

    # --- free per-region memory ---
    rm(i_region_dt, mats, M, Cov, M2, Cov2, BS, BS_sm, tstat, dmrs, stat_vec)
    gc(FALSE)
  } # end for each region in chunk

  # --- free per-chunk memory ---
  rm(i_chunk, hdr, all_cols, file_sample_ids, keep_ids, sel_cols, chunk_regions)
  gc(FALSE)
} # end for each chunk

# Close and finish
cat("Done. Outputs:\n",
    "  DMR table: ", dmr_out_path, "\n",
    "  summary:   ", log_summary, "\n", sep = "")

# References:
# bsseq Bioconductor docs: https://bioconductor.org/packages/bsseq/
# data.table fread(select): https://rdatatable.gitlab.io/data.table/reference/fread.html
# BSmooth method (Hansen et al., 2012): https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-r83

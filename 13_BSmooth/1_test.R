#!/usr/bin/env Rscript
# BSmooth on a single-region workflow aligned to your mgcv pipeline
# Inputs/paths mirror your mgcv script; grouping by AA_only.
# If *_cov columns are absent, uses a default coverage (20) as a fallback (logged).

suppressPackageStartupMessages({
  library(bsseq)
  library(data.table)
  library(dplyr)
  library(stringr)
})

## ---- Config (mirror your mgcv script) ---------------------------------
min_cpgs   <- 10
N_jobs     <- 300
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
#PATH_wk    <- "~/scratch/UQAC/meth/"
PATH_scr11 <- file.path(PATH_wk, "scr/11_mgcv/")
PATH_output <- paste0(PATH_wk, "results/13_bsmooth/test/")

#PATH_output <- file.path(PATH_wk, "results/12_bsmooth/BMI_B1/")
if (!dir.exists(PATH_output)) dir.create(PATH_output, recursive = TRUE)

## ---- Args --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("No SLURM_ARRAY_TASK_ID provided.")
i_job_id <- as.integer(args[1])

#!!!!
i_job_id <- 1
cat("BSmooth job id:", i_job_id, "\n")

## ---- Load region index -------------------------------------------------
region_file <- fread(file.path(PATH_scr11, "dat/region_file_1_chunk.csv"))
vec_list <- split(seq_len(nrow(region_file)),
                  cut(seq_len(nrow(region_file)), breaks = N_jobs, labels = FALSE))
i_region_id_vector <- vec_list[[i_job_id]]

## ---- Load phenotype (AA_only used for groups) --------------------------
load(file = file.path(PATH_scr11, "dat/18_pheno_BMI.RData"), verbose = TRUE)
# expect 'pheno_file' with columns: ID, AA_only, plus covariates
pheno_file <- as_tibble(pheno_file) %>% select(ID, AA_only)

## ---- Helpers -----------------------------------------------------------
# build M and Cov matrices from a region data.table with *_meth and *_cov columns
.build_bsseq_matrices <- function(dt_region) {
  meth_cols <- grep("_meth$", names(dt_region), value = TRUE)
  if (length(meth_cols) == 0) stop("No *_meth columns found in region chunk.")
  samp_ids <- str_replace(meth_cols, "_meth$", "")

  # match *_cov
  cov_cols <- paste0(samp_ids, "_tot")
  have_cov <- cov_cols %in% names(dt_region)

  # coerce numeric
  dt_region[, (meth_cols) := lapply(.SD, as.numeric), .SDcols = meth_cols]
  if (any(have_cov)) {
    dt_region[, (cov_cols[have_cov]) := lapply(.SD, as.numeric), .SDcols = cov_cols[have_cov]]
  }

  # order CpGs by genomic start
  setorder(dt_region, start)

  # assemble matrices [rows = CpGs, cols = samples]
  # assemble matrices [rows = CpGs, cols = samples]
  Meth <- as.matrix(dt_region[, ..meth_cols])
  if (all(have_cov)) {
    Cov <- as.matrix(dt_region[, ..cov_cols])
  } else {
    warning("Coverage columns not found for all samples; using constant coverage = 20.")
    Cov <- matrix(20, nrow = nrow(Meth), ncol = ncol(Meth))
  }

  # ensure bounds & integer counts
  Meth[is.na(Meth)] <- 0
  Cov[is.na(Cov)]   <- 0
  Meth <- pmin(pmax(Meth, 0), 1)
  M_counts <- round(Meth * Cov)

  # *** IMPORTANT: set assay colnames to the bare sample IDs ***
  colnames(M_counts) <- samp_ids
  colnames(Cov)      <- samp_ids

  list(M = M_counts, Cov = Cov, sample_names = samp_ids)
}

## ---- Main loop over regions in this job --------------------------------
tstat_out_path <- file.path(PATH_output, sprintf("bsmooth_tstat_job_%d.rds", i_job_id))
dmr_out_path   <- file.path(PATH_output, sprintf("dmr_job_%d.tsv.gz", i_job_id))
log_summary    <- file.path(PATH_output, sprintf("summary_job_%d.tsv", i_job_id))

summary_conn <- file(log_summary, open = "wt")
writeLines(paste(
  "region_id", "data_chunk_id", "chr", "region_start", "region_end",
  "n_cpg", "n_samples", "n_case", "n_ctrl", "n_dmrs"
, sep = "\t"), con = summary_conn)

# container for all tstat objects if you want to keep them per-job
tstat_list <- list()

#!!!!!!! delte after test
i_region <- 57
i_region_info <- region_file[i_region, ]
i_chunk_id <- i_region_info$data_chunk_id


i_region <- i_region_id_vector[1]  # for testing
for (i_region in i_region_id_vector) {

  #!!!!! recover this after test run
  i_info     <- region_file[i_region, ]
  i_chunk_id <- i_info$data_chunk_id

  # load chunk!!!!
  chunk_fp <- file.path(PATH_wk, "dat/", sprintf("chunk_%04d.csv", i_chunk_id))
  #chunk_fp <- file.path(PATH_wk, "data/meth_split", sprintf("chunk_%04d.csv", i_chunk_id))

  if (!file.exists(chunk_fp)) {
    cat("Missing chunk:", chunk_fp, " — skipping region", i_region, "\n")
    next
  }

  i_chunk <- fread(chunk_fp)

  # subset to region
  i_region_dt <- i_chunk %>%
    filter(chr == i_info$chr,
           start >= i_info$region_start,
           start <= i_info$region_end) %>%
    as.data.table()

  if (nrow(i_region_dt) < min_cpgs) {
    cat("Region", i_region, "insufficient CpGs (", nrow(i_region_dt), ") — skip\n")
    next
  }

  # build bsseq matrices
  mats <- .build_bsseq_matrices(i_region_dt)
  M   <- mats$M
  Cov <- mats$Cov
  samp_ids <- mats$sample_names

  # intersect with phenotype and build groups
  ph_sub <- pheno_file %>% filter(ID %in% samp_ids)
  if (nrow(ph_sub) < 2) {
    cat("Region", i_region, "— too few samples with phenotype — skip\n")
    next
  }
  # match sample order
  keep <- samp_ids %in% ph_sub$ID
  if (!any(keep)) {
    cat("Region", i_region, "— no overlap between samples and phenotype — skip\n")
    next
  }

  # reorder to phenotype order for clarity
  idx <- match(ph_sub$ID, samp_ids)
  idx <- idx[!is.na(idx)]

  M2   <- M[, idx, drop = FALSE]
  Cov2 <- Cov[, idx, drop = FALSE]
  samp_ids2 <- samp_ids[idx]
  grp_vec <- ph_sub$AA_only[match(samp_ids2, ph_sub$ID)]
  grp_vec <- as.integer(grp_vec) # 1 = case, 0 = control

  # reorder phenotype rows to the exact column order of the matrices
  ph_sub <- ph_sub[match(samp_ids2, ph_sub$ID), , drop = FALSE]

  # sanity checks
  stopifnot(identical(colnames(M2), samp_ids2))
  stopifnot(identical(colnames(Cov2), samp_ids2))



  n_case <- sum(grp_vec == 1, na.rm = TRUE)
  n_ctrl <- sum(grp_vec == 0, na.rm = TRUE)
  if (n_case < 1 || n_ctrl < 1) {
    cat("Region", i_region, "— need at least 1 case and 1 control — skip\n")
    next
  }

  # build BSseq object for this region
  BS <- BSseq(M = M2,
              Cov = Cov2,
              chr = rep(i_info$chr, nrow(M2)),
              pos = i_region_dt$start[order(i_region_dt$start)],
              sampleNames = samp_ids2)


  # smooth
  mc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  if (is.na(mc) || mc < 1) mc <- 1
  BS_sm <- BSmooth(BS, ns = 70, h = 1000, maxGap = 1e08, verbose = FALSE)


  # group indices
  g1 <- which(grp_vec == 1)
  g2 <- which(grp_vec == 0)

  # t-statistics
  tstat <- BSmooth.tstat(BS_sm,
                         group1 = g1,
                         group2 = g2,
                         estimate.var = "group2",
                         local.correct = TRUE,
                         verbose = FALSE)

  # DMRs
  dmrs <- dmrFinder(tstat, cutoff = c(-4.6, 4.6))  # ~p<4e-6 two-sided if large df; tune as you wish
  n_dmrs <- if (is.null(dmrs) || nrow(dmrs) == 0) 0L else nrow(dmrs)

  # annotate and write out (append)
  if (n_dmrs > 0) {
    dmrs$region_id    <- i_region
    dmrs$data_chunk_id<- i_chunk_id
    dmrs$chr_region   <- i_info$chr
    dmrs$reg_start    <- i_info$region_start
    dmrs$reg_end      <- i_info$region_end
    dmrs$n_cpgs     <-  nrow(M2)
    dmrs$n_samples  <-  ncol(M2)
    
        # --- add ranking/plotting columns ---
    # Prefer areaStat if present; else fall back to maxStat; else NA
    if ("areaStat" %in% names(dmrs)) {
      dmrs$score <- abs(dmrs$areaStat)
    } else if ("maxStat" %in% names(dmrs)) {
      dmrs$score <- abs(dmrs$maxStat)
    } else {
      dmrs$score <- NA_real_
    }

    # Direction of effect: use sign of areaStat if available; else meanDiff; else maxStat
    dmrs$direction <- if ("areaStat" %in% names(dmrs)) {
      sign(dmrs$areaStat)
    } else if ("meanDiff" %in% names(dmrs)) {
      sign(dmrs$meanDiff)
    } else if ("maxStat" %in% names(dmrs)) {
      sign(dmrs$maxStat)
    } else {
      NA_integer_
    }

    # Midpoint for Manhattan-like plotting
    dmrs$midpoint <- floor((dmrs$start + dmrs$end) / 2)

    # Optional: region length (bp) and CpG count already present as 'width' / 'n' in bsseq output
    if (!("width" %in% names(dmrs)) && all(c("start","end") %in% names(dmrs))) {
      dmrs$width <- (dmrs$end - dmrs$start + 1L)
    }

    
    fwrite(dmrs,
           file = dmr_out_path,
           append = file.exists(dmr_out_path),
           sep = "\t")
  }

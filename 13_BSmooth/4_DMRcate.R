#!/usr/bin/env Rscript
# Detect DMRs related to allergic asthma (AA_only 0/1) using DMRcate (sequencing mode)
# Layout & IO mirror your BSmooth script; single-threaded & per-chunk/per-region friendly.

# ---- Memory tracker (optional) ----
mem_rss_mb <- function() {
  p <- readLines("/proc/self/status")
  as.numeric(sub(".*:\\s+([0-9]+).*","\\1", grep("^VmRSS:", p, value=TRUE)))/1024
}

suppressPackageStartupMessages({
  library(DMRcate)
  library(BiocParallel)
  library(bsseq)          # to build BSseq objects from counts
  library(edgeR)          # modelMatrixMeth
  library(data.table); data.table::setDTthreads(1)
  library(dplyr)
  library(stringr)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

# ---- Hard-stop parallelism ----
BiocParallel::register(BiocParallel::SerialParam())
options(mc.cores = 1)

# ---- Config ----
min_cpgs   <- 10
N_jobs     <- 300
PATH_wk    <- "~/scratch/UQAC/meth/"
#！！！
i_job = 1
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_scr11 <- file.path(PATH_wk, "scr/11_mgcv/")
PATH_out   <- file.path(PATH_wk, "results/13_dmrcate/test")
dir.create(PATH_out, recursive = TRUE, showWarnings = FALSE)

# ---- Args ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("No SLURM_ARRAY_TASK_ID provided.")
i_job <- as.integer(args[1])
cat("DMRcate job id:", i_job, "\n")

# ---- Inputs ----
region_file <- data.table::fread(file.path(PATH_scr11, "dat/region_file_1_chunk.csv"))
vec_list <- split(seq_len(nrow(region_file)),
                  cut(seq_len(nrow(region_file)), breaks = N_jobs, labels = FALSE))
i_region_ids <- vec_list[[i_job]]
if (length(i_region_ids) == 0 || all(is.na(i_region_ids))) {
  cat("No regions assigned to this job; exiting.\n"); q("no")
}

# Phenotype: expect 'pheno_file' with columns ID, AA_only (0/1)
load(file = file.path(PATH_scr11, "dat/bsmooth_only_pheno_bmi.RData"), verbose = TRUE)
pheno_file <- tibble::as_tibble(pheno_file) %>% dplyr::select(ID, AA_only)

# ---- Helpers ----
# Convert a region data.table to BSseq (counts-based)
build_bsseq_from_region <- function(dt_region) {
  meth_cols <- grep("_meth$", names(dt_region), value = TRUE)
  if (!length(meth_cols)) stop("No *_meth columns found.")
  samp_ids <- str_replace(meth_cols, "_meth$", "")
  cov_cols <- paste0(samp_ids, "_tot")

  # numeric coercion
  dt_region[, (meth_cols) := lapply(.SD, as.numeric), .SDcols = meth_cols]
  have_cov <- cov_cols %in% names(dt_region)
  if (any(have_cov)) {
    dt_region[, (cov_cols[have_cov]) := lapply(.SD, as.numeric), .SDcols = cov_cols[have_cov]]
  }

  setorder(dt_region, start)
  Meth <- as.matrix(dt_region[, ..meth_cols])
  Cov  <- if (all(have_cov)) as.matrix(dt_region[, ..cov_cols]) else {
    warning("Missing some coverage columns; using constant coverage = 20 for those.")
    C <- matrix(20, nrow = nrow(Meth), ncol = ncol(Meth))
    if (any(have_cov)) C[, have_cov] <- as.matrix(dt_region[, ..cov_cols[have_cov]])
    C
  }

  Meth[is.na(Meth)] <- 0
  Cov[ is.na(Cov) ] <- 0
  Meth <- pmin(pmax(Meth, 0), 1)              # if stored as proportions
  Mcnt <- round(Meth * Cov)                    # integer methylated counts

  colnames(Mcnt) <- samp_ids
  colnames(Cov)  <- samp_ids

  BSseq(
    M = Mcnt,
    Cov = Cov,
    chr = dt_region$chr,
    pos = dt_region$start,
    sampleNames = samp_ids
  )
}

# ---- Outputs ----
dmr_out   <- file.path(PATH_out, sprintf("dmrcate_job_%d.tsv", i_job))
sum_out   <- file.path(PATH_out, sprintf("summary_job_%d.tsv", i_job))
if (!file.exists(sum_out)) {
  cat("region_id","data_chunk_id","chr","region_start","region_end",
      "n_cpg","n_samples","n_case","n_ctrl","n_dmrs",
      sep="\t", file=sum_out); cat("\n", file=sum_out, append=TRUE)
}

# ---- Subset to regions owned by this job and split by chunk ----
regions_this_job <- region_file[i_region_ids]
regions_this_job$region_id <- i_region_ids
chunks <- split(regions_this_job, regions_this_job$data_chunk_id)

cat("** DMRcate sequencing mode **\n")


#!!!!
chunk_name <- names(chunks)[1]
for (chunk_name in names(chunks)) {

  chunk_regions <- chunks[[chunk_name]]
  chunk_id_num  <- as.integer(chunk_name)

  #!!!!
  chunk_fp <- file.path(PATH_wk, "dat/", sprintf("chunk_%04d.csv", chunk_id_num))

  #chunk_fp <- file.path(PATH_wk, "data/meth_split", sprintf("chunk_%04d.csv", chunk_id_num))
  if (!file.exists(chunk_fp)) { cat("Missing chunk:", chunk_fp, "\n"); next }

  # Only read needed columns for this chunk
  hdr <- data.table::fread(chunk_fp, nrows=0, showProgress=FALSE)
  all_cols <- names(hdr)
  file_sample_ids <- gsub("_meth$|_tot$", "", setdiff(all_cols, c("chr","start","end")))
  keep_ids <- intersect(file_sample_ids, pheno_file$ID)
  if (!length(keep_ids)) { cat("No overlapping samples in chunk", chunk_id_num, "\n"); next }
  sel_cols <- intersect(c("chr","start","end", paste0(keep_ids,"_meth"), paste0(keep_ids,"_tot")), all_cols)

  cat(sprintf("Reading chunk %s with %d selected columns (%d samples)\n",
              chunk_id_num, length(sel_cols), length(keep_ids)))
  i_chunk <- data.table::fread(chunk_fp, select=sel_cols, showProgress=FALSE)
  i_chunk[, chr := as.character(chr)]


  row_i = 1
  for (row_i in seq_len(nrow(chunk_regions))) {
    i_info <- chunk_regions[row_i, ]
    i_info$chr <- as.character(i_info$chr)

    # Slice region
    i_region <- i_chunk[chr == i_info$chr &
                        start >= i_info$region_start &
                        start <= i_info$region_end]

    if (nrow(i_region) < min_cpgs) {
      cat("Region", i_info$region_id, ": CpGs <", min_cpgs, "— skip\n"); next
    }

    # Build BSseq for this region
    bs <- build_bsseq_from_region(i_region)

    # Match phenotype and build design!!!!
    samp_ids <- colnames(SummarizedExperiment::assay(bs, "M"))


    ph_sub   <- pheno_file %>% filter(ID %in% samp_ids, !is.na(AA_only))
    idx      <- match(ph_sub$ID, samp_ids); idx <- idx[!is.na(idx)]
    if (length(idx) < 2) { cat("Region", i_info$region_id, ": <2 phenotyped samples — skip\n"); next }

    # Reorder BSseq to phenotype
    bs <- bs[, idx]
    grp <- as.integer(ph_sub$AA_only)  # 0/1

    # Ensure at least one case and one control
    n_case <- sum(grp == 1); n_ctrl <- sum(grp == 0)
    if (n_case < 1 || n_ctrl < 1) { cat("Region", i_info$region_id, ": need 1+ case & 1+ control — skip\n"); next }

    # ---- DMRcate sequencing path ----
    # Regular design (group indicator), then convert to methylation model
    design <- model.matrix(~ grp)
    colnames(design) <- c("Intercept","AA_only")
    rownames(design) <- colnames(bs)

    methdesign <- edgeR::modelMatrixMeth(design)  # per DMRcate sequencing vignette

    # Contrast: AA_only effect (the "AA_only" column from 'design')
    cont.mat <- limma::makeContrasts(AA_vs_ctrl = AA_only, levels = methdesign)

    ## ---- coverage prefilter at the data.table stage ----
    # choose a modest requirement: at least 'min_cov' in at least 'min_samples' samples
    min_cov <- 5
    min_samples <- max(2, floor(0.7 * length(keep_ids)))   # e.g., 70% of samples

    # Build a logical matrix for coverage columns
    cov_cols <- paste0(samp_ids, "_tot")    # or derive from your keep_ids
    cov_cols <- cov_cols[cov_cols %in% names(i_region)]
    if (length(cov_cols) > 0) {
    keep_row <- rowSums(as.matrix(i_region[, ..cov_cols]) >= min_cov, na.rm = TRUE) >= min_samples
    i_region <- i_region[keep_row]
    }

    # Need at least 2 CpGs for voom/DMRcate
    if (nrow(i_region) < 2) {
    cat("Region", i_info$region_id, ": < 2 CpGs after coverage filter — skip\n")
    next
    }

    # Allow partial coverage to avoid collapsing to 0–1 CpGs
    seq_annot <- sequencing.annotate(
    obj      = bs,
    methdesign      = methdesign,
    contrasts   = TRUE,
    cont.matrix = cont.mat,
    coef        = "AA_vs_ctrl",
    all.cov     = FALSE,      # <-- key change
    fdr         = 0.05        # screening FDR; loosen (e.g., 0.25) if still too sparse
    )


    # Region calling
    # 'lambda' is kernel bandwidth in bp (default 1000); 'C' scales bandwidth by effect;
    # 'min.cpgs' is minimum CpGs per region.
    dmrc <- dmrcate(seq_annot, lambda = 1000, C = 2, min.cpgs = 3)

    # Extract GRanges of DMRs (set genome if you want UCSC-style seqlevels)
    dmrs_gr <- extractRanges(dmrc, genome = "hg19")

    n_dmrs <- length(dmrs_gr)
    cat(sprintf("Region %d: DMRcate found %d DMRs\n", i_info$region_id, n_dmrs))

    # Write summary row
    cat(i_info$region_id, i_info$data_chunk_id, i_info$chr,
        i_info$region_start, i_info$region_end,
        nrow(i_region), ncol(bs), n_case, n_ctrl, n_dmrs,
        sep="\t", file=sum_out, append=TRUE); cat("\n", file=sum_out, append=TRUE)

    # Write DMRs (append, TSV)
    if (n_dmrs > 0) {
      # Add your region-ID context and simple score
      dmrt <- as.data.table(dmrs_gr)
      # DMRcate returns columns like: no.cpgs, minfdr, Stouffer, etc.
      if (!"no.cpgs" %in% names(dmrt)) dmrt[, no.cpgs := NA_integer_]
      if (!"minfdr"  %in% names(dmrt)) dmrt[, minfdr  := NA_real_]
      dmrt[, `:=`(
        region_id     = i_info$region_id,
        data_chunk_id = i_info$data_chunk_id,
        chr_region    = i_info$chr,
        reg_start     = i_info$region_start,
        reg_end       = i_info$region_end,
        n_cpgs_window = nrow(i_region),
        n_samples     = ncol(bs),
        score         = -log10(pmax(minfdr, .Machine$double.xmin)) # larger is stronger
      )]

      fwrite(dmrt,
             file = file.path(PATH_out, sprintf("dmrcate_job_%d.tsv", i_job)),
             sep = "\t", append = file.exists(file.path(PATH_out, sprintf("dmrcate_job_%d.tsv", i_job))),
             col.names = !file.exists(file.path(PATH_out, sprintf("dmrcate_job_%d.tsv", i_job))))
    }

    # cleanup
    rm(i_region, bs, methdesign, seq_annot, dmrc, dmrs_gr); gc(FALSE)
  } # for each region
  rm(i_chunk, hdr, all_cols, file_sample_ids, keep_ids, sel_cols, chunk_regions); gc(FALSE)
} # for each chunk

cat("Done. Outputs in: ", PATH_out, "\n")

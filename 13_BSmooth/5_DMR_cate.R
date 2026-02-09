#!/usr/bin/env Rscript
# Detect DMRs related to allergic asthma (AA_only 0/1) using DMRcate (sequencing mode)
# Layout & IO mirror your BSmooth script; single-threaded & per-chunk/per-region friendly.

safe_extract_gr <- function(dmrc_obj,
                            cutoff_fdr = 0.05,     # stricter than 0.05
                            min_cpgs   = 10,       # require ≥10 CpGs per DMR
                            min_abs_md = 0.05,     # min |mean methyl diff| (5 pp)
                            min_width  = 100) {    # min genomic width (bp)
  # Try official helper first
  gr_try <- try(DMRcate::extractRanges(dmrc_obj, genome = "hg19",
                                       cutoff = cutoff_fdr,
                                       min.cpgs = min_cpgs),
                silent = TRUE)
  if (!inherits(gr_try, "try-error")) {
    # Optionally filter further by effect size / width
    if (!is.null(S4Vectors::mcols(gr_try)$meandiff)) {
      gr_try <- gr_try[abs(S4Vectors::mcols(gr_try)$meandiff) >= min_abs_md]
    }
    gr_try <- gr_try[IRanges::width(gr_try) >= min_width]
    return(gr_try)
  }

  # ---- Fallback: DMResults -> GRanges ----
  if (methods::is(dmrc_obj, "DMResults")) {
    coords   <- methods::slot(dmrc_obj, "coord")
    no.cpgs  <- methods::slot(dmrc_obj, "no.cpgs")
    minfdr   <- methods::slot(dmrc_obj, "min_smoothed_fdr")
    stouffer <- if ("Stouffer" %in% slotNames(dmrc_obj)) methods::slot(dmrc_obj, "Stouffer") else NA_real_
    hmfdr    <- if ("HMFDR"    %in% slotNames(dmrc_obj)) methods::slot(dmrc_obj, "HMFDR")    else NA_real_
    fisher   <- if ("Fisher"   %in% slotNames(dmrc_obj)) methods::slot(dmrc_obj, "Fisher")   else NA_real_
    maxdiff  <- if ("maxdiff"  %in% slotNames(dmrc_obj)) methods::slot(dmrc_obj, "maxdiff")  else NA_real_
    meandiff <- if ("meandiff" %in% slotNames(dmrc_obj)) methods::slot(dmrc_obj, "meandiff") else NA_real_

    parse_coord <- function(z) {
      m <- regexec("^([^:]+):(\\d+)-(\\d+)$", z)
      r <- regmatches(z, m)[[1]]
      if (length(r) != 4L) return(c(NA, NA, NA))
      c(r[2], r[3], r[4])
    }
    parsed <- t(vapply(as.character(coords), parse_coord, FUN.VALUE = character(3)))
    chr   <- parsed[,1]
    start <- suppressWarnings(as.integer(parsed[,2]))
    end   <- suppressWarnings(as.integer(parsed[,3]))

    gr <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges   = IRanges::IRanges(start = start, end = end)
    )
    S4Vectors::mcols(gr)$no.cpgs  <- as.integer(no.cpgs)
    S4Vectors::mcols(gr)$minfdr   <- as.numeric(minfdr)
    S4Vectors::mcols(gr)$Stouffer <- as.numeric(stouffer)
    S4Vectors::mcols(gr)$HMFDR    <- as.numeric(hmfdr)
    S4Vectors::mcols(gr)$Fisher   <- as.numeric(fisher)
    S4Vectors::mcols(gr)$maxdiff  <- as.numeric(maxdiff)
    S4Vectors::mcols(gr)$meandiff <- as.numeric(meandiff)

    # ---- Hard filters here ----
    keep <- !is.na(S4Vectors::mcols(gr)$minfdr) & S4Vectors::mcols(gr)$minfdr <= cutoff_fdr
    keep <- keep & S4Vectors::mcols(gr)$no.cpgs >= min_cpgs
    if (!all(is.na(S4Vectors::mcols(gr)$meandiff))) {
      keep <- keep & abs(S4Vectors::mcols(gr)$meandiff) >= min_abs_md
    }
    keep <- keep & IRanges::width(gr) >= min_width
    return(gr[keep])
  }

  if (is.list(dmrc_obj) && !is.null(dmrc_obj$results)) {
    res <- dmrc_obj$results
    if (!NROW(res)) return(GenomicRanges::GRanges())
    cn <- tolower(colnames(res)); pick <- function(...) { keys <- c(...); i <- match(keys, cn, nomatch=0L); res[[ which(i != 0L)[1] ]] }
    seqv  <- pick("seqnames","chr","chrom","chromosome")
    start <- as.integer(pick("start","cstart","startpos","chromstart"))
    end   <- as.integer(pick("end","cend","endpos","chromend"))
    gr <- GenomicRanges::GRanges(seqnames = seqv, ranges = IRanges::IRanges(start, end))
    if ("no.cpgs" %in% colnames(res)) S4Vectors::mcols(gr)$no.cpgs <- res$no.cpgs
    if ("minfdr"  %in% colnames(res)) S4Vectors::mcols(gr)$minfdr  <- res$minfdr
    if ("meandiff" %in% colnames(res)) S4Vectors::mcols(gr)$meandiff <- res$meandiff

    keep <- rep(TRUE, length(gr))
    if ("minfdr" %in% colnames(res))  keep <- keep & res$minfdr <= cutoff_fdr
    if ("no.cpgs" %in% colnames(res)) keep <- keep & res$no.cpgs >= min_cpgs
    if ("meandiff" %in% colnames(res)) keep <- keep & abs(res$meandiff) >= min_abs_md
    keep <- keep & IRanges::width(gr) >= min_width
    return(gr[keep])
  }

  GenomicRanges::GRanges()
}




# ---- Memory tracker (optional) ----
mem_rss_mb <- function() {
  p <- readLines("/proc/self/status")
  as.numeric(sub(".*:\\s+([0-9]+).*","\\1", grep("^VmRSS:", p, value=TRUE)))/1024
}

suppressPackageStartupMessages({
  library(DMRcate)
  library(DMRcatedata)   # for example data; ensure installed
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
##!!!
#i_job = 1
#PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_scr11 <- file.path(PATH_wk, "scr/11_mgcv/")
#PATH_out   <- file.path(PATH_wk, "results/13_dmrcate/test")
PATH_out   <- file.path(PATH_wk, "results/13_dmrcate/stringent")

#dir.create(PATH_out, recursive = TRUE, showWarnings = FALSE)

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
#chunk_name <- names(chunks)[1]
for (chunk_name in names(chunks)) {

  chunk_regions <- chunks[[chunk_name]]
  chunk_id_num  <- as.integer(chunk_name)

  #!!!!
  #chunk_fp <- file.path(PATH_wk, "dat/", sprintf("chunk_%04d.csv", chunk_id_num))

  chunk_fp <- file.path(PATH_wk, "data/meth_split", sprintf("chunk_%04d.csv", chunk_id_num))
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
    #design <- model.matrix(~ grp)
    #colnames(design) <- c("Intercept","AA_only")
    #rownames(design) <- colnames(bs)

    #methdesign <- edgeR::modelMatrixMeth(design)  # per DMRcate sequencing vignette

    # Contrast: AA_only effect (the "AA_only" column from 'design')
    #cont.mat <- limma::makeContrasts(AA_vs_ctrl = AA_only, levels = methdesign)


    #!!!!
    design <- model.matrix(~ AA_only + Sex + AgeCalc + BMI + sv1 + sv2 + sv3 + sv4 + sv5, data = ph_sub)
    colnames(design)[2] <- "AA_only"   # ensure the second term is AA_only
    rownames(design) <- colnames(bs)

    methdesign <- edgeR::modelMatrixMeth(design)
    cont.mat   <- limma::makeContrasts(AA_vs_ctrl = AA_only, levels = methdesign)


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
    fdr         = 0.05       # screening FDR; loosen (e.g., 0.25) if still too sparse !!!
    )


    # --- EARLY EXIT if no individually significant CpGs (avoids dmrcate error) ---
    a <- try(seq_annot@ranges, silent = TRUE)   # or methods::slot(seq_annot, "annot")

    ## After: seq_annot <- sequencing.annotate(...)

    # --- Get per-CpG annotation from the CpGannotated object (S4) ---
    rng <- try(methods::slot(seq_annot, "ranges"), silent = TRUE)

    # If 'ranges' slot is missing or empty, skip
    if (inherits(rng, "try-error") || is.null(rng) || length(rng) == 0L) {
    cat("Region", i_info$region_id, ": no CpG annotation present — skip\n")
    cat(i_info$region_id, i_info$data_chunk_id, i_info$chr,
        i_info$region_start, i_info$region_end,
        nrow(i_region), ncol(bs), n_case, n_ctrl, 0L,
        sep = "\t", file = sum_out, append = TRUE); cat("\n", file = sum_out, append = TRUE)
    rm(rng, seq_annot, bs); gc(FALSE); next
    }

    # Pull metadata as data.frame
    md <- as.data.frame(S4Vectors::mcols(rng))
    cols <- names(md)

    print("cols")
    print(cols)
    
    print(md)

    print(md$is.sig)

    print(sum(md$is.sig, na.rm = TRUE) > 0L)

    # Prefer the package's own seed flag if present
    has_seed <- FALSE
    if ("is.sig" %in% cols) {
    has_seed <- sum(md$is.sig, na.rm = TRUE) > 0L
    } else {
    # Fallback: compute BH-FDR from available p-values
    p <- if ("ind.p" %in% cols) md$ind.p else if ("rawpval" %in% cols) md$rawpval else NA_real_
    if (is.numeric(p) && any(is.finite(p))) {
        q <- p.adjust(p, method = "BH")
        # Use the same screening FDR you passed to sequencing.annotate()
        has_seed <- sum(q <= 0.05, na.rm = TRUE) > 0L  # adjust 0.05 if you used another fdr
    }
    }

    if (!has_seed) {
    cat("Region", i_info$region_id, ": no individually significant CpGs at screening FDR — skip\n")
    cat(i_info$region_id, i_info$data_chunk_id, i_info$chr,
        i_info$region_start, i_info$region_end,
        nrow(i_region), ncol(bs), n_case, n_ctrl, 0L,
        sep = "\t", file = sum_out, append = TRUE); cat("\n", file = sum_out, append = TRUE)
    rm(rng, md, seq_annot, bs); gc(FALSE); next
    }

    # ... proceed to dmrcate(...)




    # --- Region calling (guarded) ---
    dmrc <- try(dmrcate(seq_annot, lambda = 1000, C = 2, min.cpgs = 10), silent = TRUE)
    if (inherits(dmrc, "try-error")) {
    cat("Region", i_info$region_id, "— dmrcate error:", as.character(dmrc), "\n")

    cat(i_info$region_id, i_info$data_chunk_id, i_info$chr,
        i_info$region_start, i_info$region_end,
        nrow(i_region), ncol(bs), n_case, n_ctrl, 0L,
        sep = "\t", file = sum_out, append = TRUE)
    cat("\n", file = sum_out, append = TRUE)

    rm(dmrc, seq_annot, bs); gc(FALSE)
    next
    }

    # Extract GRanges of DMRs (hg19 per your data build)
    #dmrs_gr <- extractRanges(dmrc, genome = "hg19")
    dmrs_gr <- safe_extract_gr(dmrc, cutoff_fdr = 0.05, min_cpgs = 10, min_abs_md = 0.05, min_width = 100)

    n_dmrs <- length(dmrs_gr)
    cat(sprintf("Region %d: DMRcate found %d DMRs\n", i_info$region_id, n_dmrs))

    # Write summary line
    cat(i_info$region_id, i_info$data_chunk_id, i_info$chr,
        i_info$region_start, i_info$region_end,
        nrow(i_region), ncol(bs), n_case, n_ctrl, n_dmrs,
        sep = "\t", file = sum_out, append = TRUE)
    cat("\n", file = sum_out, append = TRUE)

    # Write DMRs (append, TSV)
    if (n_dmrs > 0) {
    dmrt <- as.data.table(dmrs_gr)
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
        score         = -log10(pmax(minfdr, .Machine$double.xmin))
    )]

    fwrite(dmrt,
            file = file.path(PATH_out, sprintf("dmrcate_job_%d.tsv", i_job)),
            sep = "\t",
            append = file.exists(file.path(PATH_out, sprintf("dmrcate_job_%d.tsv", i_job))),
            col.names = !file.exists(file.path(PATH_out, sprintf("dmrcate_job_%d.tsv", i_job))))
    }

    # cleanup
    rm(i_region, bs, methdesign, seq_annot, dmrc, dmrs_gr); gc(FALSE)
  } # for each region
  rm(i_chunk, hdr, all_cols, file_sample_ids, keep_ids, sel_cols, chunk_regions); gc(FALSE)
} # for each chunk

cat("Done. Outputs in: ", PATH_out, "\n")

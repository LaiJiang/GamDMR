
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DMRcate)
  library(BiocParallel)
  library(bsseq)
  library(edgeR)
  library(limma)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(GenomicRanges)
  library(S4Vectors)
})

data.table::setDTthreads(1L)
BiocParallel::register(BiocParallel::SerialParam())
options(mc.cores = 1L)

min_cpgs <- 10L
min_cov <- 5
n_jobs <- 300L
PATH_wk <- path.expand(Sys.getenv("METH_BASE_DIR", "~/scratch/UQAC/meth"))
PATH_scr11 <- file.path(PATH_wk, "scr", "11_mgcv")
output_dir <- file.path(PATH_wk, "results", "13_dmrcate")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) stop("Usage: Rscript dmrcate_region_clean.R <job_id>")
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

env <- new.env(parent = emptyenv())
load(file.path(PATH_scr11, "dat", "bsmooth_pheno_only.RData"), envir = env)
if (!exists("pheno_file", envir = env, inherits = FALSE)) stop("pheno_file is missing.")
pheno_file <- as.data.table(get("pheno_file", envir = env))
rm(env)

design_vars <- c(
  "ID", "AA_only", "Sex", "AgeCalc", "Non.smoker",
  "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc", "BMI",
  "sv1", "sv2", "sv3", "sv4", "sv5"
)
missing_vars <- setdiff(design_vars, names(pheno_file))
if (length(missing_vars)) stop("Missing phenotype columns: ", paste(missing_vars, collapse = ", "))
pheno_file <- pheno_file[, ..design_vars]
pheno_file[, ID := as.character(ID)]

build_bsseq <- function(region_dt, sample_ids) {
  meth_cols <- paste0(sample_ids, "_meth")
  cov_cols <- paste0(sample_ids, "_tot")
  if (!all(meth_cols %in% names(region_dt)) || !all(cov_cols %in% names(region_dt))) return(NULL)

  x <- copy(region_dt)
  x[, (meth_cols) := lapply(.SD, as.numeric), .SDcols = meth_cols]
  x[, (cov_cols) := lapply(.SD, as.numeric), .SDcols = cov_cols]
  setorder(x, start)

  meth <- as.matrix(x[, ..meth_cols])
  cov <- as.matrix(x[, ..cov_cols])
  storage.mode(meth) <- "numeric"
  storage.mode(cov) <- "numeric"

  valid <- is.finite(cov) & cov >= min_cov & is.finite(meth)
  cov[!valid] <- 0
  meth[!valid] <- 0
  meth <- pmin(pmax(meth, 0), 1)
  mcnt <- round(meth * cov)

  keep <- rowSums(cov > 0) > 0
  x <- x[keep]
  cov <- cov[keep, , drop = FALSE]
  mcnt <- mcnt[keep, , drop = FALSE]
  if (nrow(x) < 2L) return(NULL)

  colnames(cov) <- sample_ids
  colnames(mcnt) <- sample_ids
  BSseq(
    M = mcnt,
    Cov = cov,
    chr = as.character(x$chr),
    pos = as.integer(x$start),
    sampleNames = sample_ids
  )
}

extract_dmrs <- function(obj) {
  gr <- tryCatch(
    DMRcate::extractRanges(obj, genome = "hg19", cutoff = 0.05, min.cpgs = 10L),
    error = function(e) GRanges()
  )
  if (!length(gr)) return(gr)
  if ("no.cpgs" %in% names(mcols(gr))) gr <- gr[mcols(gr)$no.cpgs >= 10L]
  if ("meandiff" %in% names(mcols(gr))) {
    gr <- gr[is.finite(mcols(gr)$meandiff) & abs(mcols(gr)$meandiff) >= 0.05]
  }
  gr
}

dmr_out <- file.path(output_dir, sprintf("dmrcate_job_%d.tsv", job_id))
summary_out <- file.path(output_dir, sprintf("summary_job_%d.tsv", job_id))
regions_this_job <- copy(region_file[region_ids])
chunks <- split(regions_this_job, regions_this_job$data_chunk_id)

summary_rows <- list()
dmr_rows <- list()

for (chunk_name in names(chunks)) {
  chunk_regions <- as.data.table(chunks[[chunk_name]])
  chunk_id <- as.integer(chunk_name)
  chunk_file <- file.path(PATH_wk, "data", "meth_split", sprintf("chunk_%04d.csv", chunk_id))
  if (!file.exists(chunk_file)) next

  header <- fread(chunk_file, nrows = 0L, showProgress = FALSE)
  sample_ids <- intersect(
    gsub("_meth$|_tot$", "", grep("_meth$|_tot$", names(header), value = TRUE)),
    pheno_file$ID
  )
  sample_ids <- sample_ids[
    paste0(sample_ids, "_meth") %in% names(header) &
      paste0(sample_ids, "_tot") %in% names(header)
  ]
  if (!length(sample_ids)) next

  selected <- c("chr", "start", "end", paste0(sample_ids, "_meth"), paste0(sample_ids, "_tot"))
  selected <- intersect(selected, names(header))
  chunk_dt <- fread(chunk_file, select = selected, showProgress = FALSE)
  chunk_dt[, chr := as.character(chr)]

  for (j in seq_len(nrow(chunk_regions))) {
    info <- chunk_regions[j]
    region_dt <- chunk_dt[
      chr == as.character(info$chr) &
        start >= info$region_start &
        start <= info$region_end
    ]
    if (uniqueN(region_dt$start) < min_cpgs) next

    ph <- pheno_file[ID %in% sample_ids]
    ph <- ph[complete.cases(ph[, ..design_vars])]
    if (nrow(ph) < 2L || uniqueN(ph$AA_only) < 2L) next

    bs <- tryCatch(build_bsseq(region_dt, ph$ID), error = function(e) NULL)
    if (is.null(bs) || nrow(bs) < 2L) next
    ph <- ph[match(colnames(bs), ID)]

    design <- model.matrix(
      ~ AA_only + Sex + AgeCalc + Non.smoker +
        EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + BMI +
        sv1 + sv2 + sv3 + sv4 + sv5,
      data = as.data.frame(ph)
    )
    colnames(design)[2L] <- "AA_only"
    rownames(design) <- colnames(bs)
    methdesign <- edgeR::modelMatrixMeth(design)
    colnames(methdesign) <- make.names(colnames(methdesign), unique = TRUE)
    if (!"AA_only" %in% colnames(methdesign)) next
    cont <- limma::makeContrasts(AA_vs_ctrl = AA_only, levels = methdesign)

    annot <- tryCatch(
      DMRcate::sequencing.annotate(
        obj = bs,
        methdesign = methdesign,
        contrasts = TRUE,
        cont.matrix = cont,
        coef = "AA_vs_ctrl",
        all.cov = FALSE,
        fdr = 0.05
      ),
      error = function(e) NULL
    )
    if (is.null(annot)) next

    rng <- tryCatch(methods::slot(annot, "ranges"), error = function(e) NULL)
    n_seed <- 0L
    if (!is.null(rng) && length(rng)) {
      md <- as.data.frame(mcols(rng))
      if ("is.sig" %in% names(md)) {
        n_seed <- sum(md$is.sig %in% TRUE, na.rm = TRUE)
      } else {
        pv <- if ("ind.p" %in% names(md)) md$ind.p else if ("rawpval" %in% names(md)) md$rawpval else numeric()
        if (length(pv)) n_seed <- sum(p.adjust(pv, "BH") <= 0.05, na.rm = TRUE)
      }
    }

    dmrs_gr <- GRanges()
    if (n_seed >= 1L) {
      dmrc <- tryCatch(
        DMRcate::dmrcate(annot, lambda = 1000, C = 2, min.cpgs = 10L),
        error = function(e) NULL
      )
      if (!is.null(dmrc)) dmrs_gr <- extract_dmrs(dmrc)
    }

    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      region_id = as.integer(info$region_id),
      data_chunk_id = as.integer(info$data_chunk_id),
      chr = as.character(info$chr),
      region_start = as.integer(info$region_start),
      region_end = as.integer(info$region_end),
      n_cpg = nrow(bs),
      n_samples = ncol(bs),
      n_case = sum(ph$AA_only == 1L),
      n_ctrl = sum(ph$AA_only == 0L),
      n_seed_cpgs = n_seed,
      n_dmrs = length(dmrs_gr)
    )

    if (length(dmrs_gr)) {
      d <- as.data.table(as.data.frame(dmrs_gr))
      d[, `:=`(
        chr = as.character(seqnames(dmrs_gr)),
        region_start = as.integer(start(dmrs_gr)),
        region_end = as.integer(end(dmrs_gr)),
        region_width = as.integer(width(dmrs_gr)),
        region_id = as.integer(info$region_id),
        data_chunk_id = as.integer(info$data_chunk_id),
        parent_chr = as.character(info$chr),
        parent_start = as.integer(info$region_start),
        parent_end = as.integer(info$region_end),
        n_samples = ncol(bs)
      )]
      dmr_rows[[length(dmr_rows) + 1L]] <- d
    }
  }
}

if (length(summary_rows)) {
  fwrite(rbindlist(summary_rows, fill = TRUE), summary_out, sep = "\t")
}
if (length(dmr_rows)) {
  fwrite(rbindlist(dmr_rows, fill = TRUE), dmr_out, sep = "\t")
}

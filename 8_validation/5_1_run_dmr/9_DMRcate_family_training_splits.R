
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DMRcate)
  library(BiocParallel)
  library(bsseq)
  library(edgeR)
  library(limma)
  library(data.table)
  library(GenomicRanges)
  library(S4Vectors)
})

data.table::setDTthreads(1L)
BiocParallel::register(BiocParallel::SerialParam())
options(mc.cores = 1L)

PATH_wk <- path.expand(Sys.getenv("METH_BASE_DIR", "~/scratch/UQAC/meth/"))
PATH_scr11 <- file.path(PATH_wk, "scr", "11_mgcv")
region_path <- file.path(PATH_wk, "results", "15_revision", "5_cv", "region_DMRcate_M12.csv")
split_path <- file.path(PATH_scr11, "dat", "100_family_training_splits.csv")
pheno_path <- file.path(PATH_scr11, "dat", "18_pheno_BMI.RData")
chunk_dir <- file.path(PATH_wk, "data", "meth_split")
PATH_out <- file.path(PATH_wk, "results", "15_revision", "5_cv", "3_dmrcate", "training_splits")
dir.create(PATH_out, recursive = TRUE, showWarnings = FALSE)

min_cpgs <- 10L
min_cov <- 5
N_jobs <- as.integer(Sys.getenv("N_JOBS", "900"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Pass SLURM_ARRAY_TASK_ID.")
job_id <- as.integer(args[1])
if (is.na(job_id)) stop("Invalid job ID.")

regions <- fread(region_path)
if (!"region_id" %in% names(regions)) regions[, region_id := .I]
if ("n_cpgs" %in% names(regions)) regions <- regions[n_cpgs >= min_cpgs]

groups <- split(
  seq_len(nrow(regions)),
  cut(seq_len(nrow(regions)), breaks = min(N_jobs, nrow(regions)), labels = FALSE)
)
if (job_id < 1L || job_id > length(groups)) quit(save = "no", status = 0L)
regions <- regions[groups[[job_id]]]

splits <- fread(split_path)
train_cols <- grep("^train_FID_[0-9]+$", names(splits), value = TRUE)
train_cols <- train_cols[order(as.integer(sub("^train_FID_", "", train_cols)))]

env <- new.env(parent = emptyenv())
load(pheno_path, envir = env)
if (!exists("pheno_file", envir = env, inherits = FALSE)) stop("pheno_file is missing.")
pheno <- as.data.table(get("pheno_file", envir = env))
rm(env)

design_vars <- c(
  "ID", "FID", "AA_only", "Sex", "AgeCalc", "Non.smoker",
  "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc", "BMI",
  "sv1", "sv2", "sv3", "sv4", "sv5"
)
missing <- setdiff(design_vars, names(pheno))
if (length(missing)) stop("Missing phenotype columns: ", paste(missing, collapse = ", "))
pheno <- pheno[, ..design_vars]
pheno[, `:=`(ID = as.character(ID), FID = as.character(FID))]

build_bs <- function(region_dt, ids) {
  mcols <- paste0(ids, "_meth")
  ccols <- paste0(ids, "_tot")
  if (!all(c(mcols, ccols) %in% names(region_dt))) return(NULL)

  x <- copy(region_dt)
  x[, (mcols) := lapply(.SD, as.numeric), .SDcols = mcols]
  x[, (ccols) := lapply(.SD, as.numeric), .SDcols = ccols]
  setorder(x, start)

  meth <- as.matrix(x[, ..mcols])
  cov <- as.matrix(x[, ..ccols])
  valid <- is.finite(cov) & cov >= min_cov & is.finite(meth)
  cov[!valid] <- 0
  meth[!valid] <- 0
  meth <- pmin(pmax(meth, 0), 1)
  M <- round(meth * cov)
  keep <- rowSums(cov > 0) > 0
  x <- x[keep]
  cov <- cov[keep, , drop = FALSE]
  M <- M[keep, , drop = FALSE]
  if (nrow(x) < 2L) return(NULL)

  colnames(M) <- ids
  colnames(cov) <- ids
  BSseq(M = M, Cov = cov, chr = as.character(x$chr), pos = x$start, sampleNames = ids)
}

summary_rows <- list()
dmr_rows <- list()

for (chunk_id in unique(regions$data_chunk_id)) {
  chunk_file <- file.path(chunk_dir, sprintf("chunk_%04d.csv", chunk_id))
  if (!file.exists(chunk_file)) next
  header <- fread(chunk_file, nrows = 0L)
  chunk_ids <- intersect(
    gsub("_meth$|_tot$", "", grep("_meth$|_tot$", names(header), value = TRUE)),
    pheno$ID
  )
  selected <- intersect(
    c("chr", "start", "end", paste0(chunk_ids, "_meth"), paste0(chunk_ids, "_tot")),
    names(header)
  )
  chunk_dt <- fread(chunk_file, select = selected)
  chunk_dt[, chr := as.character(chr)]

  for (j in which(regions$data_chunk_id == chunk_id)) {
    info <- regions[j]
    region_dt <- chunk_dt[
      chr == as.character(info$chr) &
        start >= info$region_start &
        start <= info$region_end
    ]
    if (uniqueN(region_dt$start) < min_cpgs) next

    for (sid in splits$splitID) {
      train_fids <- unique(as.character(unlist(splits[splitID == sid, ..train_cols], use.names = FALSE)))
      train_fids <- train_fids[!is.na(train_fids) & nzchar(train_fids)]
      ph <- pheno[FID %in% train_fids & ID %in% chunk_ids]
      ph <- ph[complete.cases(ph[, ..design_vars])]
      if (nrow(ph) < 2L || uniqueN(ph$AA_only) < 2L) next

      bs <- tryCatch(build_bs(region_dt, ph$ID), error = function(e) NULL)
      if (is.null(bs)) next
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
        sequencing.annotate(
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

      gr <- GRanges()
      if (n_seed >= 1L) {
        dmrc <- tryCatch(dmrcate(annot, lambda = 1000, C = 2, min.cpgs = 10L), error = function(e) NULL)
        if (!is.null(dmrc)) {
          gr <- tryCatch(
            extractRanges(dmrc, genome = "hg19", cutoff = 0.05, min.cpgs = 10L),
            error = function(e) GRanges()
          )
          if (length(gr) && "no.cpgs" %in% names(mcols(gr))) gr <- gr[mcols(gr)$no.cpgs >= 10L]
          if (length(gr) && "meandiff" %in% names(mcols(gr))) {
            gr <- gr[is.finite(mcols(gr)$meandiff) & abs(mcols(gr)$meandiff) >= 0.05]
          }
        }
      }

      summary_rows[[length(summary_rows) + 1L]] <- data.table(
        prep_job_id = job_id,
        region_id = as.integer(info$region_id),
        data_chunk_id = as.integer(info$data_chunk_id),
        chr = as.character(info$chr),
        region_start = as.integer(info$region_start),
        region_end = as.integer(info$region_end),
        splitID = as.integer(sid),
        status = "FIT_OK",
        n_train_FIDs_requested = length(train_fids),
        n_train_FIDs_observed = uniqueN(ph$FID),
        n_train_samples = nrow(ph),
        n_train_cases = sum(ph$AA_only == 1L),
        n_train_controls = sum(ph$AA_only == 0L),
        n_cpgs_after_coverage = nrow(bs),
        n_seed_cpgs = n_seed,
        n_dmrs = length(gr)
      )

      if (length(gr)) {
        d <- as.data.table(as.data.frame(gr))
        d[, `:=`(
          prep_job_id = job_id,
          region_id = as.integer(info$region_id),
          data_chunk_id = as.integer(info$data_chunk_id),
          parent_chr = as.character(info$chr),
          parent_start = as.integer(info$region_start),
          parent_end = as.integer(info$region_end),
          splitID = as.integer(sid),
          chr = as.character(seqnames(gr)),
          region_start = as.integer(start(gr)),
          region_end = as.integer(end(gr))
        )]
        dmr_rows[[length(dmr_rows) + 1L]] <- d
      }
    }
  }
}

summary_out <- file.path(PATH_out, sprintf("dmrcate_training_summary_job_%04d.tsv", job_id))
dmr_out <- file.path(PATH_out, sprintf("dmrcate_training_dmrs_job_%04d.tsv", job_id))
if (length(summary_rows)) fwrite(rbindlist(summary_rows, fill = TRUE), summary_out, sep = "\t")
if (length(dmr_rows)) fwrite(rbindlist(dmr_rows, fill = TRUE), dmr_out, sep = "\t")

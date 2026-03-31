#!/usr/bin/env Rscript

safe_extract_gr <- function(dmrc_obj,
                            cutoff_fdr = 0.01,
                            min_cpgs = 10,
                            min_abs_md = 0.05,
                            min_width = 250) {
  gr_try <- try(
    DMRcate::extractRanges(dmrc_obj, genome = "hg19", cutoff = cutoff_fdr, min.cpgs = min_cpgs),
    silent = TRUE
  )
  if (!inherits(gr_try, "try-error")) {
    if (!is.null(S4Vectors::mcols(gr_try)$meandiff)) {
      gr_try <- gr_try[abs(S4Vectors::mcols(gr_try)$meandiff) >= min_abs_md]
    }
    gr_try <- gr_try[IRanges::width(gr_try) >= min_width]
    return(gr_try)
  }

  if (methods::is(dmrc_obj, "DMResults")) {
    coords <- methods::slot(dmrc_obj, "coord")
    no.cpgs <- methods::slot(dmrc_obj, "no.cpgs")
    minfdr <- methods::slot(dmrc_obj, "min_smoothed_fdr")
    maxdiff <- if ("maxdiff" %in% slotNames(dmrc_obj)) methods::slot(dmrc_obj, "maxdiff") else NA_real_
    meandiff <- if ("meandiff" %in% slotNames(dmrc_obj)) methods::slot(dmrc_obj, "meandiff") else NA_real_

    parse_coord <- function(z) {
      m <- regexec("^([^:]+):(\\d+)-(\\d+)$", z)
      r <- regmatches(z, m)[[1]]
      if (length(r) != 4L) return(c(NA, NA, NA))
      c(r[2], r[3], r[4])
    }

    parsed <- t(vapply(as.character(coords), parse_coord, FUN.VALUE = character(3)))
    chr <- parsed[, 1]
    start <- suppressWarnings(as.integer(parsed[, 2]))
    end <- suppressWarnings(as.integer(parsed[, 3]))

    gr <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges = IRanges::IRanges(start = start, end = end)
    )
    S4Vectors::mcols(gr)$no.cpgs <- as.integer(no.cpgs)
    S4Vectors::mcols(gr)$minfdr <- as.numeric(minfdr)
    S4Vectors::mcols(gr)$maxdiff <- as.numeric(maxdiff)
    S4Vectors::mcols(gr)$meandiff <- as.numeric(meandiff)

    keep <- !is.na(S4Vectors::mcols(gr)$minfdr) & S4Vectors::mcols(gr)$minfdr <= cutoff_fdr
    keep <- keep & S4Vectors::mcols(gr)$no.cpgs >= min_cpgs
    if (!all(is.na(S4Vectors::mcols(gr)$meandiff))) {
      keep <- keep & abs(S4Vectors::mcols(gr)$meandiff) >= min_abs_md
    }
    keep <- keep & IRanges::width(gr) >= min_width
    return(gr[keep])
  }

  GenomicRanges::GRanges()
}

suppressPackageStartupMessages({
  library(DMRcate)
  library(BiocParallel)
  library(bsseq)
  library(edgeR)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

data.table::setDTthreads(1)
BiocParallel::register(BiocParallel::SerialParam())
options(mc.cores = 1)

min_cpgs <- 10
n_jobs <- 300

base_dir <- Sys.getenv("METH_BASE_DIR", unset = getwd())
script_dir <- file.path(base_dir, "11_mgcv")
output_dir <- file.path(base_dir, "results", "13_dmrcate")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript dmrcate_region_clean.R <job_id>")
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

load(file.path(script_dir, "dat", "bsmooth_pheno_only.RData"))
pheno_file <- tibble::as_tibble(pheno_file)

build_bsseq_from_region <- function(dt_region) {
  meth_cols <- grep("_meth$", names(dt_region), value = TRUE)
  if (length(meth_cols) == 0) stop("No *_meth columns found.")
  sample_ids <- str_replace(meth_cols, "_meth$", "")
  cov_cols <- paste0(sample_ids, "_tot")

  dt_region[, (meth_cols) := lapply(.SD, as.numeric), .SDcols = meth_cols]
  have_cov <- cov_cols %in% names(dt_region)
  if (any(have_cov)) {
    dt_region[, (cov_cols[have_cov]) := lapply(.SD, as.numeric), .SDcols = cov_cols[have_cov]]
  }

  setorder(dt_region, start)
  meth <- as.matrix(dt_region[, ..meth_cols])
  cov <- if (all(have_cov)) {
    as.matrix(dt_region[, ..cov_cols])
  } else {
    warning("Some coverage columns are missing; using constant coverage = 20 where needed.")
    cmat <- matrix(20, nrow = nrow(meth), ncol = ncol(meth))
    if (any(have_cov)) cmat[, have_cov] <- as.matrix(dt_region[, ..cov_cols[have_cov]])
    cmat
  }

  meth[is.na(meth)] <- 0
  cov[is.na(cov)] <- 0
  meth <- pmin(pmax(meth, 0), 1)
  mcnt <- round(meth * cov)

  colnames(mcnt) <- sample_ids
  colnames(cov) <- sample_ids

  BSseq(
    M = mcnt,
    Cov = cov,
    chr = dt_region$chr,
    pos = dt_region$start,
    sampleNames = sample_ids
  )
}

dmr_out <- file.path(output_dir, sprintf("dmrcate_job_%d.tsv", job_id))
summary_out <- file.path(output_dir, sprintf("summary_job_%d.tsv", job_id))
if (!file.exists(summary_out)) {
  cat(
    "region_id", "data_chunk_id", "chr", "region_start", "region_end",
    "n_cpg", "n_samples", "n_case", "n_ctrl", "n_dmrs",
    sep = "\t", file = summary_out, append = FALSE
  )
  cat("\n", file = summary_out, append = TRUE)
}

regions_this_job <- copy(region_file[region_ids])
regions_this_job$region_id <- region_ids
chunks <- split(regions_this_job, regions_this_job$data_chunk_id)

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

    i_region <- i_chunk[
      chr == i_info$chr & start >= i_info$region_start & start <= i_info$region_end
    ]

    if (nrow(i_region) < min_cpgs) next

    sample_ids <- gsub("_meth$", "", grep("_meth$", names(i_region), value = TRUE))
    ph_sub <- pheno_file %>% filter(ID %in% sample_ids, !is.na(AA_only))
    if (nrow(ph_sub) < 2) next

    bs <- tryCatch(build_bsseq_from_region(i_region), error = function(e) NULL)
    if (is.null(bs)) next

    bs_sample_ids <- colnames(SummarizedExperiment::assay(bs, "M"))
    ph_sub <- ph_sub %>% filter(ID %in% bs_sample_ids)
    idx <- match(ph_sub$ID, bs_sample_ids)
    idx <- idx[!is.na(idx)]
    if (length(idx) < 2) {
      rm(bs)
      gc(FALSE)
      next
    }

    bs <- bs[, idx]
    ph_sub <- ph_sub[match(colnames(bs), ph_sub$ID), , drop = FALSE]
    grp <- as.integer(ph_sub$AA_only)

    n_case <- sum(grp == 1, na.rm = TRUE)
    n_ctrl <- sum(grp == 0, na.rm = TRUE)
    if (n_case < 1 || n_ctrl < 1) {
      rm(bs)
      gc(FALSE)
      next
    }

    design <- model.matrix(
      ~ AA_only + Sex + AgeCalc + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + BMI + sv1 + sv2 + sv3 + sv4 + sv5,
      data = ph_sub
    )
    colnames(design)[2] <- "AA_only"
    rownames(design) <- colnames(bs)

    methdesign <- edgeR::modelMatrixMeth(design)
    colnames(methdesign) <- make.names(colnames(methdesign), unique = TRUE)
    if (!make.names("AA_only") %in% colnames(methdesign)) {
      rm(bs)
      gc(FALSE)
      next
    }

    cont_mat <- limma::makeContrasts(AA_vs_ctrl = AA_only, levels = methdesign)

    cov_cols <- paste0(colnames(bs), "_tot")
    cov_cols <- cov_cols[cov_cols %in% names(i_region)]
    if (length(cov_cols) > 0) {
      min_cov <- 5
      min_samples <- max(2, floor(0.7 * ncol(bs)))
      keep_row <- rowSums(as.matrix(i_region[, ..cov_cols]) >= min_cov, na.rm = TRUE) >= min_samples
      i_region <- i_region[keep_row]
      if (nrow(i_region) < 2) {
        rm(bs)
        gc(FALSE)
        next
      }
    }

    seq_annot <- tryCatch(
      sequencing.annotate(
        obj = bs,
        methdesign = methdesign,
        contrasts = TRUE,
        cont.matrix = cont_mat,
        coef = "AA_vs_ctrl",
        all.cov = FALSE,
        fdr = 0.05
      ),
      error = function(e) NULL
    )
    if (is.null(seq_annot)) {
      rm(bs)
      gc(FALSE)
      next
    }

    rng <- tryCatch(methods::slot(seq_annot, "ranges"), error = function(e) NULL)
    if (is.null(rng) || length(rng) == 0L) {
      cat(
        i_info$region_id, i_info$data_chunk_id, i_info$chr,
        i_info$region_start, i_info$region_end,
        nrow(i_region), ncol(bs), n_case, n_ctrl, 0L,
        sep = "\t", file = summary_out, append = TRUE
      )
      cat("\n", file = summary_out, append = TRUE)
      rm(bs, seq_annot)
      gc(FALSE)
      next
    }

    md <- as.data.frame(S4Vectors::mcols(rng))
    cols <- names(md)

    has_seed <- FALSE
    if ("is.sig" %in% cols) {
      has_seed <- sum(md$is.sig, na.rm = TRUE) >= 5
    } else {
      p <- if ("ind.p" %in% cols) md$ind.p else if ("rawpval" %in% cols) md$rawpval else NA_real_
      if (is.numeric(p) && any(is.finite(p))) {
        q <- p.adjust(p, method = "BH")
        has_seed <- sum(q <= 0.05, na.rm = TRUE) >= 5
      }
    }

    if (!has_seed) {
      cat(
        i_info$region_id, i_info$data_chunk_id, i_info$chr,
        i_info$region_start, i_info$region_end,
        nrow(i_region), ncol(bs), n_case, n_ctrl, 0L,
        sep = "\t", file = summary_out, append = TRUE
      )
      cat("\n", file = summary_out, append = TRUE)
      rm(bs, seq_annot, rng, md)
      gc(FALSE)
      next
    }

    dmrc <- tryCatch(dmrcate(seq_annot, lambda = 1000, C = 2, min.cpgs = 10), error = function(e) NULL)
    if (is.null(dmrc)) {
      cat(
        i_info$region_id, i_info$data_chunk_id, i_info$chr,
        i_info$region_start, i_info$region_end,
        nrow(i_region), ncol(bs), n_case, n_ctrl, 0L,
        sep = "\t", file = summary_out, append = TRUE
      )
      cat("\n", file = summary_out, append = TRUE)
      rm(bs, seq_annot, rng, md)
      gc(FALSE)
      next
    }

    dmrs_gr <- safe_extract_gr(dmrc, cutoff_fdr = 0.01, min_cpgs = 10, min_abs_md = 0.05, min_width = 250)
    n_dmrs <- length(dmrs_gr)

    cat(
      i_info$region_id, i_info$data_chunk_id, i_info$chr,
      i_info$region_start, i_info$region_end,
      nrow(i_region), ncol(bs), n_case, n_ctrl, n_dmrs,
      sep = "\t", file = summary_out, append = TRUE
    )
    cat("\n", file = summary_out, append = TRUE)

    if (n_dmrs > 0) {
      dmrt <- as.data.table(dmrs_gr)
      if (!"no.cpgs" %in% names(dmrt)) dmrt[, no.cpgs := NA_integer_]
      if (!"minfdr" %in% names(dmrt)) dmrt[, minfdr := NA_real_]
      dmrt[, `:=`(
        region_id = i_info$region_id,
        data_chunk_id = i_info$data_chunk_id,
        chr_region = i_info$chr,
        reg_start = i_info$region_start,
        reg_end = i_info$region_end,
        n_cpgs_window = nrow(i_region),
        n_samples = ncol(bs),
        score = -log10(pmax(minfdr, .Machine$double.xmin))
      )]

      fwrite(
        dmrt,
        file = dmr_out,
        sep = "\t",
        append = file.exists(dmr_out),
        col.names = !file.exists(dmr_out)
      )
    }

    rm(i_region, bs, methdesign, seq_annot, dmrc, dmrs_gr, rng, md)
    gc(FALSE)
  }

  rm(i_chunk, hdr, all_cols, file_sample_ids, keep_ids, sel_cols, chunk_regions)
  gc(FALSE)
}

cat("Done.\n")
cat("DMRs:", dmr_out, "\n")
cat("Summary:", summary_out, "\n")

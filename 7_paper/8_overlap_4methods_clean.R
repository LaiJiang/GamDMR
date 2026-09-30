suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(VennDiagram)
})

base_dir <- path.expand(Sys.getenv("METH_BASE_DIR", "~/scratch/UQAC/meth"))
paper_dir <- file.path(base_dir, "scr", "14_paper")
output_dir <- file.path(paper_dir, "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

somnibus_file <- file.path(base_dir, "results", "9_regional", "16_somnibus_regions.csv")
mgcv_rdata <- file.path(output_dir, "4_DMR_M1_M2_overlap.RData")
bsmooth_file <- file.path(output_dir, "6_BSmooth_DMRs.csv")
dmrcate_file <- file.path(output_dir, "7_DMRcate_DMRs.csv")

somnibus_dmrs <- fread(somnibus_file)
mgcv_env <- new.env(parent = emptyenv())
load(mgcv_rdata, envir = mgcv_env)
if (!exists("dmr_gene_file", envir = mgcv_env, inherits = FALSE)) stop("dmr_gene_file is missing.")
mgcv_dmrs <- as.data.table(get("dmr_gene_file", envir = mgcv_env))
rm(mgcv_env)
bsmooth_dmrs <- fread(bsmooth_file)
dmrcate_dmrs <- fread(dmrcate_file)

normalize_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)

to_gr <- function(df, chr_col, start_col, end_col) {
  x <- as.data.table(df)
  z <- data.table(
    chr = normalize_chr(x[[chr_col]]),
    start = suppressWarnings(as.integer(x[[start_col]])),
    end = suppressWarnings(as.integer(x[[end_col]]))
  )
  z <- unique(z[!is.na(chr) & nzchar(chr) & is.finite(start) & is.finite(end) & start <= end])
  GRanges(
    seqnames = z$chr,
    ranges = IRanges(start = z$start, end = z$end)
  )
}

gr_list <- list(
  SOMNiBUS = to_gr(somnibus_dmrs, "chr", "region_start", "region_end"),
  GAM_DMR = to_gr(mgcv_dmrs, "chr", "region_start", "region_end"),
  BSmooth = to_gr(bsmooth_dmrs, "chr", "region_start", "region_end"),
  DMRcate = to_gr(dmrcate_dmrs, "chr", "region_start", "region_end")
)

methods <- names(gr_list)
count_mat <- matrix(0L, nrow = length(methods), ncol = length(methods), dimnames = list(methods, methods))
percent_mat <- matrix(0, nrow = length(methods), ncol = length(methods), dimnames = list(methods, methods))
hit_tables <- list()

for (i in seq_along(methods)) {
  query <- gr_list[[i]]
  for (j in seq_along(methods)) {
    subject <- gr_list[[j]]
    if (i == j) {
      count_mat[i, j] <- length(query)
      percent_mat[i, j] <- if (length(query)) 100 else 0
      next
    }
    hits <- findOverlaps(query, subject, ignore.strand = TRUE)
    q <- unique(queryHits(hits))
    count_mat[i, j] <- length(q)
    percent_mat[i, j] <- if (length(query)) 100 * length(q) / length(query) else 0
    if (length(hits)) {
      qh <- queryHits(hits)
      sh <- subjectHits(hits)
      ov_start <- pmax(start(query)[qh], start(subject)[sh])
      ov_end <- pmin(end(query)[qh], end(subject)[sh])
      hit_tables[[paste(methods[i], methods[j], sep = "__")]] <- data.table(
        query_method = methods[i],
        subject_method = methods[j],
        query_index = qh,
        subject_index = sh,
        chr = as.character(seqnames(query)[qh]),
        overlap_start = ov_start,
        overlap_end = ov_end,
        overlap_width = ov_end - ov_start + 1L
      )
    }
  }
}

count_dt <- cbind(data.table(query_method = rownames(count_mat)), as.data.table(count_mat))
percent_dt <- cbind(data.table(query_method = rownames(percent_mat)), as.data.table(percent_mat))
fwrite(count_dt, file.path(output_dir, "8_overlap_matrix.csv"))
fwrite(percent_dt, file.path(output_dir, "8_overlap_percent.csv"))
if (length(hit_tables)) {
  fwrite(rbindlist(hit_tables, use.names = TRUE, fill = TRUE), file.path(output_dir, "8_direct_pairwise_hits.csv"))
}

gr_intersection <- function(x, y) {
  hits <- findOverlaps(x, y, ignore.strand = TRUE)
  if (!length(hits)) return(GRanges())
  qh <- queryHits(hits)
  sh <- subjectHits(hits)
  z <- GRanges(
    seqnames = as.character(seqnames(x)[qh]),
    ranges = IRanges(
      start = pmax(start(x)[qh], start(y)[sh]),
      end = pmin(end(x)[qh], end(y)[sh])
    )
  )
  reduce(z, ignore.strand = TRUE)
}

reduced <- lapply(gr_list, reduce, ignore.strand = TRUE)
common <- Reduce(gr_intersection, reduced)
common_df <- data.table(
  chr = as.character(seqnames(common)),
  start = start(common),
  end = end(common),
  width = width(common)
)
fwrite(common_df, file.path(output_dir, "8_common_regions_all_methods.csv"))

atoms <- disjoin(unlist(GRangesList(reduced), use.names = FALSE), ignore.strand = TRUE)
atom_ids <- paste0(seqnames(atoms), ":", start(atoms), "-", end(atoms))
set_list <- lapply(reduced, function(gr) atom_ids[countOverlaps(atoms, gr, ignore.strand = TRUE) > 0L])

venn.diagram(
  x = set_list,
  category.names = c("SOMNiBUS", "GAM-DMR", "BSmooth", "DMRcate"),
  filename = file.path(output_dir, "8_DMR_4methods_venn.png"),
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 600
)

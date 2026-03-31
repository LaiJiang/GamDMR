suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(VennDiagram)
})

base_dir <- Sys.getenv("METH_BASE_DIR", unset = getwd())
paper_dir <- file.path(base_dir, "14_paper")
output_dir <- file.path(paper_dir, "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

somnibus_file <- file.path(base_dir, "results", "9_regional", "16_somnibus_regions.csv")
mgcv_rdata <- file.path(paper_dir, "results", "4_DMR_M1_M2_overlap.RData")
bsmooth_file <- file.path(output_dir, "6_BSmooth_DMRs.csv")
dmrcate_file <- file.path(output_dir, "7_DMRcate_DMRs.csv")

somnibus_dmrs <- fread(somnibus_file)
load(mgcv_rdata)
mgcv_dmrs <- dmr_gene_file
bsmooth_dmrs <- fread(bsmooth_file)
dmrcate_dmrs <- fread(dmrcate_file)

to_gr <- function(df, chr_col, start_col, end_col) {
  GRanges(
    seqnames = as.character(df[[chr_col]]),
    ranges = IRanges(
      start = as.integer(df[[start_col]]),
      end = as.integer(df[[end_col]])
    )
  )
}

gr_somni <- to_gr(somnibus_dmrs, "chr", "region_start", "region_end")
gr_mgcv <- to_gr(mgcv_dmrs, "chr", "region_start", "region_end")
gr_bs <- to_gr(bsmooth_dmrs, "chr", "region_start", "region_end")
gr_dmrc <- to_gr(dmrcate_dmrs, "chr", "region_start", "region_end")

method_gr_list <- GRangesList(
  SOMNiBUS = gr_somni,
  GAM_DMR = gr_mgcv,
  BSmooth = gr_bs,
  DMRcate = gr_dmrc
)

all_gr <- unlist(method_gr_list, use.names = FALSE)
all_regions <- reduce(all_gr)

region_ids <- paste0(seqnames(all_regions), ":", start(all_regions), "-", end(all_regions))
names(all_regions) <- region_ids

set_list <- lapply(method_gr_list, function(gr) {
  hits <- findOverlaps(all_regions, gr)
  unique(queryHits(hits))
})
set_list_ids <- lapply(set_list, function(idx) region_ids[idx])

method_names <- names(set_list_ids)
overlap_mat <- matrix(
  0L,
  nrow = length(method_names),
  ncol = length(method_names),
  dimnames = list(method_names, method_names)
)

for (i in seq_along(method_names)) {
  for (j in seq_along(method_names)) {
    overlap_mat[i, j] <- length(intersect(set_list_ids[[i]], set_list_ids[[j]]))
  }
}

write.csv(as.data.frame(overlap_mat), file.path(output_dir, "8_overlap_matrix.csv"))

venn.diagram(
  x = set_list_ids,
  category.names = c("SOMNiBUS", "GAM-DMR", "BSmooth", "DMRcate"),
  filename = file.path(output_dir, "8_DMR_4methods_venn.png"),
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 600,
  compression = "lzw"
)

common_ids <- Reduce(intersect, set_list_ids)
common_region <- all_regions[names(all_regions) %in% common_ids]
common_region_df <- data.frame(
  chr = as.character(seqnames(common_region)),
  start = start(common_region),
  end = end(common_region)
)

write.csv(common_region_df, file.path(output_dir, "8_common_regions_all_methods.csv"), row.names = FALSE)

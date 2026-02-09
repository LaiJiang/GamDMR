PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_scr14 <- paste0(PATH_wk, "/scr/14_paper/")


library(data.table)
library(dplyr)



#first load somnibus DMR results 


somnibus_DMRs <- fread(paste0(PATH_wk,"/results/9_regional/16_somnibus_regions.csv"))


#then load mgcv regions
load( file = paste0(PATH_wk,"/scr/14_paper/results/4_DMR_M1_M2_overlap.RData"),verbose=TRUE)

MGCV_DMRs <- dmr_gene_file 


#then BSmooth DMRs
BSmooth_DMRs <- fread(paste0(PATH_scr14, "results/6_BSmooth_DMRs.csv"))


#then DMRcate DMRs
DMRcate_DMRs <- fread(paste0(PATH_scr14, "results/7_DMRcate_DMRs.csv"))

#############
library(GenomicRanges)
library(data.table)
library(VennDiagram)

## ------------------------------------------------------------
## Convert each DMR table to GRanges safely
## ------------------------------------------------------------

library(GenomicRanges)

method_gr_list <- GRangesList(
  SOMNiBUS = gr_somni,
  MGCV     = gr_mgcv,
  BSmooth  = gr_bs,
  DMRcate  = gr_dmrc
)

all_gr <- unlist(method_gr_list, use.names = FALSE)

all_regions <- reduce(all_gr)



# Give each disjoint segment a unique ID
region_ids <- paste0(seqnames(all_regions), ":", 
                     start(all_regions), "-", 
                     end(all_regions))
names(all_regions) <- region_ids

## ------------------------------------------------------------
## 3. For each method, get which disjoint segments it covers
## ------------------------------------------------------------
set_list <- lapply(method_gr_list, function(gr) {
  hits <- findOverlaps(all_regions, gr)
  unique(queryHits(hits))  # indices of all_regions that overlap this method
})

# Turn indices into actual region IDs for set operations
set_list_ids <- lapply(set_list, function(idx) region_ids[idx])
names(set_list_ids) <- names(method_gr_list)

## ------------------------------------------------------------
## 4. Pairwise overlap counts between methods
##    Here: count = number of disjoint loci covered by both methods
## ------------------------------------------------------------
method_names <- names(set_list_ids)
n_methods <- length(method_names)

overlap_mat <- matrix(0L, nrow = n_methods, ncol = n_methods,
                      dimnames = list(method_names, method_names))

for (i in seq_len(n_methods)) {
  for (j in seq_len(n_methods)) {
    overlap_mat[i, j] <- length(intersect(set_list_ids[[i]], set_list_ids[[j]]))
  }
}

cat("Pairwise overlap (number of disjoint loci covered by both methods):\n")
print(overlap_mat)

## ------------------------------------------------------------
## 5. Make a 4-way Venn diagram
##    Using the disjoint-loci-based sets above
## ------------------------------------------------------------
venn_file <- paste0(PATH_scr14, "results/DMR_4methods_venn.jpeg")

print(method_names)

method_names <- c("SOMNiBUS", "GAM-DMR", "BSmooth", "DMRcate")
venn.diagram(
  x = set_list_ids,
  category.names = method_names,
  filename = venn_file,
  output = TRUE,
  imagetype = "png",
  height = 3000, width = 3000, resolution = 600,
  compression = "lzw"
)

cat("Venn diagram saved to:", venn_file, "\n")



##########################################################
#show the DMR all methods agree:


## ------------------------------------------------------------
## 6. Extract the region shared by all 4 methods
## ------------------------------------------------------------

# region IDs (disjoint segments) shared by all 4 methods
common_ids <- Reduce(intersect, set_list_ids)

cat("Number of regions shared by all 4 methods:", length(common_ids), "\n")
print(common_ids)

# Subset the GRanges for that common region
common_region <- all_regions[names(all_regions) %in% common_ids]

# Convert to a small data.frame with chr / start / end
common_region_df <- data.frame(
  chr   = as.character(seqnames(common_region)),
  start = start(common_region),
  end   = end(common_region)
)

cat("Common DMR region shared by all 4 methods:\n")
print(common_region_df)

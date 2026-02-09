# 
#Concordane test II: Systematic quantification across all tested regions:



#this file collects the mode l result of individual CpG result
#and compare that with DMR regions.
library(dplyr)
library(data.table)
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
outdir <- paste0(PATH_wk, "results/11_mgcv/")
PATH_scr11 <- paste0(PATH_wk, "/scr/11_mgcv/")
#first load M1 model CpG postion results
load( paste0(PATH_wk,"/scr/12_cpg_DMR/results/5_eval.RData"), verbose=TRUE)
# M2_manhattan

#the updated M1 signals
M2_cpg_signals <- M2_manhattan[M2_manhattan$pval < 1e-5,]

print(dim(M2_cpg_signals))
#then load the M1 model genes


#M1: gene list 
cpgs_int <- M2_manhattan$CpG[M2_manhattan$pval < 1e-5]
cpgs_int <- cpgs_int[!duplicated(cpgs_int)]
#source(paste0(PATH_wk,"/scr/8_rerun/0_func_enrichment.R"))
#use the updated annotations
source(paste0(PATH_wk,"/scr/8_rerun/0_func_enrichment_updates.R"))
M2_genes<- gene_map


"VARS1" %in% M2_genes$SYMBOL



#M3: gene list : NA.
cpgs_int <- M3_manhattan$CpG[M3_manhattan$pval < 1e-5]
cpgs_int <- cpgs_int[!duplicated(cpgs_int)]
#source(paste0(PATH_wk,"/scr/8_rerun/0_func_enrichment.R"))
#use the updated annotations
source(paste0(PATH_wk,"/scr/8_rerun/0_func_enrichment_updates.R"))
M3_genes<- gene_map


save(M2_genes, file = paste0(PATH_wk,"/scr/12_cpg_DMR/results/6_M2_genes.RData"))

#load( file = paste0(PATH_wk,"/scr/8_rerun/results/6_GO_genes.RData"),verbose=TRUE)
# M2_genes


#now load DMR region information

DMRs <- fread(   file.path(outdir, "24_dmrs_STRICT.tsv"))

region_file <- fread(paste0(PATH_scr11, "dat/region_file_1_chunk.csv"))
region_file$region_id <- 1:nrow(region_file)
# Merge with results for only intersecting region_ids
DMRs <- merge(DMRs, region_file, by = "region_id", all.x = TRUE)
#added pvals, to chek these negative and 0 pvalues should be removed or not?
DMRs <- DMRs %>%
  dplyr::select(region_id, chr, region_start, region_end,pvals)



# 1) DMR gene list (from your pipeline; must contain 'genes_symbol' column)
dmr_gene_file  <- fread(paste0(PATH_wk, "results/11_mgcv/25_DMR_strict_genes.csv"))

##############################################################################################################
#first print the number of genes in single CpG, and DMR regions, seperately.
#and then the overlapping:

print("number of unique genes from CpG model M2:")
print(length(unique(M2_genes$SYMBOL)))

print("number of unique genes from DMR regions:")
print(length(unique(dmr_gene_file$genes_symbol)))

print("number of overlapping genes:")
print(length(intersect(unique(M2_genes$SYMBOL), unique(dmr_gene_file$genes_symbol))))

#####################################################################################################################################

#now for cpgs in M2_cpg_signals, how many are overlapping with at least one of the region in DMRs?
library(data.table)
library(stringr)

# Assume your objects already exist:
# DMRs            (columns: region_id, chr, region_start, region_end, pvals)
# M2_cpg_signals  (columns include: CpG, pval, coef, pval_old)

## ---- Prep DMR intervals ----
dmr_dt <- as.data.table(DMRs)[, .(chr = as.integer(chr),
                                  start = as.integer(region_start),
                                  end   = as.integer(region_end),
                                  region_id)]
setkey(dmr_dt, chr, start, end)

## ---- Parse CpG IDs into intervals ----
# CpG format examples: "1:11023676-NA", "1:25173094-25173095"
cpg_dt <- as.data.table(M2_cpg_signals)[, .(CpG, pval, coef, pval_old)]

# Split "chr:pos1-pos2"
chr_part <- tstrsplit(cpg_dt$CpG, ":", fixed = TRUE)
pos_part <- tstrsplit(chr_part[[2]], "-", fixed = TRUE)

# Clean to integers; treat "NA" as missing
chr_num <- suppressWarnings(as.integer(chr_part[[1]]))
pos1    <- suppressWarnings(as.integer(pos_part[[1]]))
pos2    <- suppressWarnings(as.integer(ifelse(pos_part[[2]] %in% c("NA",""), NA, pos_part[[2]])))

# Build CpG intervals: if pos2 is NA -> single-point interval
cpg_iv <- data.table(
  idx   = seq_len(nrow(cpg_dt)),
  chr   = chr_num,
  start = ifelse(is.na(pos2), pos1, pmin(pos1, pos2)),
  end   = ifelse(is.na(pos2), pos1, pmin(pos1, pos2))
)

# Drop rows with malformed coordinates
cpg_iv <- cpg_iv[!is.na(chr) & !is.na(start) & !is.na(end)]

setkey(cpg_iv, chr, start, end)

## ---- Overlap via non-equi interval join ----
# foverlaps requires both tables keyed on (chr, start, end)
hits <- foverlaps(cpg_iv, dmr_dt, nomatch = 0L)

# Mark overlaps back on the original order
overlap_idx <- unique(hits$idx)
overlap_flag <- integer(nrow(cpg_dt))
overlap_flag[overlap_idx] <- 1L

# Add the new column to M2_cpg_signals (aligned by original row order)
M2_cpg_signals$overlap_DMR <- overlap_flag

## Optional: peek
# table(M2_cpg_signals$overlap_DMR)
# head(M2_cpg_signals[, c("CpG","overlap_DMR")])
##################################################################
library(data.table)

# --- assuming previous parsing already done ---
# dmr_dt : chr, start, end, region_id
# cpg_iv : chr, start, end, idx   (CpG intervals)

# Make sure keys are set
setkey(dmr_dt, chr, start, end)
setkey(cpg_iv, chr, start, end)

# Find overlaps: DMRs vs CpGs
hits2 <- foverlaps(dmr_dt, cpg_iv, nomatch = 0L)

# Which DMR rows have at least one CpG hit
overlap_dmr_ids <- unique(hits2$region_id)

# Add column to original DMRs data.table
DMRs$overlap_CpG <- as.integer(DMRs$region_id %in% overlap_dmr_ids)

## Check
table(DMRs$overlap_CpG)
head(DMRs[, c("region_id","chr","region_start","region_end","overlap_CpG")])

#now attach the  columns of chr, region_start, region_end, pvals, overlap_CpG to dmr_gene_file by matching region_id
dmr_gene_file <- merge(dmr_gene_file, DMRs[, .(region_id, chr, region_start, region_end, pvals, overlap_CpG)], by = "region_id", all.x = TRUE)

save(M2_cpg_signals, DMRs, cpg_iv, dmr_gene_file, dmr_dt, file = paste0(PATH_wk,"/scr/14_paper/results/4_overlap_M2.RData"))

#try, m1, m2, m3
#ty bfore and after introducing BMI adjustment 

dmr_gene_file[dmr_gene_file$overlap_CpG==1,]
M2_cpg_signals[M2_cpg_signals$overlap_DMR==1,]


show(cpg_iv[cpg_iv$chr==1,])

show(dmr_dt[dmr_dt$chr==1,])

#the final ovelapping results saved in two dataframes: 
#dmr_gene_file
# M2_cpg_signals

#now compare with M1 results: 
DMR_M2 <- DMRs


load( file = paste0(PATH_wk,"/scr/12_cpg_DMR/results/6_overlap_results.RData"),verbose=TRUE)

DMR_M1 <- DMRs


DMR_M1[DMR_M1$overlap_CpG==1,]
DMR_M2[DMR_M2$overlap_CpG==1,]

save(DMR_M1, DMR_M2, dmr_gene_file,  file = paste0(PATH_wk,"/scr/14_paper/results/4_DMR_M1_M2_overlap.RData"))

DMR_M2_hits <- DMR_M2[DMR_M2$overlap_CpG==1,]


dmr_gene_file[dmr_gene_file$region_id %in% DMR_M2_hits$region_id,]


sum(dmr_gene_file$genes_symbol %in% c("IRAK2", "SPRY2"))
sum(dmr_gene_file$genes_symbol == "SPRY2")

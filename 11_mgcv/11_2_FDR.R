#this script we caluate FDR adjusted pvalues, and select the DMRs. check their biological meanings

# Load libraries
library(data.table)
library(dplyr)
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_scr11 <- paste0(PATH_wk, "scr/11_mgcv/")

# Read the results table
results_all_batch <- fread(paste0(PATH_wk, "results/11_mgcv/results_all_batch.txt"))  # Replace with actual path

#remove the negative p-values
results_all_batch <- results_all_batch[results_all_batch$`s(start):AA_only` >= 0, ]



# Extract raw p-values
raw_pvals <- results_all_batch$`s(start):AA_only`



# Compute FDR-adjusted p-values using Benjamini-Hochberg
fdr_pvals <- p.adjust(raw_pvals, method = "BH")

# Add to original data
results_all_batch$pval_FDR   <- fdr_pvals

# Save updated table
#fwrite(results_all_batch, "path/to/results_with_FDR.csv")  # Replace with desired path

# Load region info
region_file <- fread(paste0(PATH_scr11, "dat/region_file_1_chunk.csv"))
region_file$region_id <- 1:nrow(region_file)
# Merge with results for only intersecting region_ids
results_merged <- merge(results_all_batch,region_file, by = "region_id", all.x = TRUE)

#first select the DMRs based on FDR p-value
dmrs_selected <- results_merged[results_merged$pval_FDR < 0.05, ]


#now we first extract the gene region of these DMRs, check their genes, if they make sense. 
#for those biological meaningful DMRs, whats their R2, pvalues, ...etc.

# Save the selected DMRs
fwrite(dmrs_selected, paste0(PATH_wk, "results/11_mgcv/11_2_dmrs_selected.csv"))
#########################################################################################################

#now extract the genes

DMRs <- dmrs_selected %>%
  dplyr::select(region_id, chr, region_start, region_end)


  source("C:/Per/LaiJiang/Project/UQAC/meth/scr/9_regional/0_func_genes_16.R")


write.csv(df_annot, file = paste0(PATH_wk, "/results/11_mgcv/11_2_DMR_genes.csv"), row.names = FALSE)
#write.csv(df_annot, file = paste0(PATH_wk, "/scr/9_regional/results/DMR_genes.csv"), row.names = FALSE)

df_annot <- read.csv(paste0(PATH_wk, "/results/11_mgcv/11_2_DMR_genes.csv"), stringsAsFactors = FALSE)
genes_int <- c("IL33","PTPRE", "ORMDL3;GSDMB", "IL1RL1", "TSLP", "IL13", "IL4","IL4R","FCER1A", "ATXN7L1;CDHR3","GATA3", "LRP1;STAT6",
"IRAK3","RORA","IKZF3","ADAM33","RUNX1","DPP10")

df_annot[which(df_annot$genes_symbol %in% genes_int),]

AA_region_id <- df_annot$region_id[which(df_annot$genes_symbol %in% genes_int)]



dmrs_AA <-  dmrs_selected[dmrs_selected$region_id %in% AA_region_id, ]
    
#now manually compare the pvalues and R2 OF dmrs_AA among dmrs_selected

#step 0: remove these regsion pvalue of s(start) is strong. it means the methylation is flat.
dmrs_AA$"s(start)"
min(dmrs_selected$"s(start)")

dmrs_AA$"s(start):AA_only"
plot(dmrs_selected$"s(start):AA_only")

#step 1: remove edf_s(start) < 1, Regions where the smooth was effectively removed (no variance)
# edf_s(start) < 1


plot(-log10((dmrs_selected$pval_FDR)),pch=19)
abline(h=-log10(dmrs_AA$pval_FDR),col="red")
#

#pvalue are not useful to enrich AA genes, they could also be 0 Pvalue.
dmrs_selected$pval_FDR[dmrs_selected$"s(start):AA_only"==0]
dmrs_AA$region_id[dmrs_AA$"s(start):AA_only"==0]

results_all_batch <- read.table(paste0(PATH_wk, "results/11_mgcv/results_all_batch.txt"), 
                                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

results_all_batch[results_all_batch$region_id %in% dmrs_AA$region_id[dmrs_AA$"s(start):AA_only"==0],]


#step 2: remove the pvalues < 0. which means the pvalue is not reliable
pvals <- results_all_batch$`s.start..AA_only`
sum(pvals <0)

dmrs_AA$"edf_s(start):AA_only"

plot(((dmrs_selected$"edf_s(start):AA_only")),pch=19)
abline(h=dmrs_AA$"edf_s(start):AA_only",col="red")

#step 3: remove the edf_s(start):AA_only < 1, which means the smooth term is not significant


dmrs_AA$"n_cpgs"
plot(((dmrs_selected$"R2")),pch=19)
abline(h=dmrs_AA$"R2",col="red")


#step 4: R2 filter by 0.2, which is the default  threshold in bumphunter, DMRcate, or metilene pipleines using GAM,
#or 0.6, to extract only strong and confient associations

plot(((dmrs_selected$"n_cpgs")),pch=19)
abline(h=dmrs_AA$"n_cpgs",col="red")

#step 5: N_samples > 300. to make sure the sample size is large enough to detect the effect

#step 6: optional, filter by N_cpgs >40
mean(dmrs_AA$"n_cpgs">40)
mean(dmrs_selected$"n_cpgs">40)

#step 1: remove edf_s(start) < 1, Regions where the smooth was effectively removed (no variance)
#step 2: remove the pvalues < 0. which means the pvalue is not reliable
#step 3: remove the edf_s(start):AA_only < 1, which means the smooth term is not significant
#step 4: R2 filter by 0.2, which is the default  threshold in bumphunter, DMRcate, or metilene pipleines using GAM,
#or 0.6, to extract only strong and confient associations   
#step 5: N_samples > 300. to make sure the sample size is large enough to detect the effect
#step 6: optional, filter by N_cpgs >40
dmrs_selected_strict <- dmrs_selected %>%
  dplyr::filter(
    `edf_s(start)` >= 1,
    `edf_s(start):AA_only` >= 1,
    pval_FDR < 0.05,
    `R2` > 0.6,
    N_samples > 300
  )

plot(((dmrs_selected_strict$"edf_s(start):AA_only")),pch=19)
abline(h=dmrs_AA$"edf_s(start):AA_only",col="red")

mean(dmrs_AA$"pval_FDR"<0.01)
mean(dmrs_selected_strict$"pval_FDR"<0.01)


mean(dmrs_AA$"pval_FDR"<5e-8)
mean(dmrs_selected_strict$"pval_FDR"<5e-8)


dmrs_selected_strict$hits <- dmrs_selected_strict$region_id %in% dmrs_AA$region_id

#rank the DMRs based on FDR p-value
dmrs_selected_strict$rank <- rank(dmrs_selected_strict$pval_FDR, ties.method = "first")

dmrs_selected_strict$rank[dmrs_selected_strict$hits==1]

sum(dmrs_selected_strict$pval_FDR<1e-5)
sum(dmrs_AA$pval_FDR<1e-5)


plot(((dmrs_selected_strict$"edf_s(start):AA_only")),dmrs_selected_strict$"N_cpgs",pch=19)


colnames(dmrs_selected_strict)

dmrs_selected_strict_merge <- dmrs_selected_strict %>% select(region_id,chr,region_start, region_end)

dmrs_selected_strict_merge$region_id[1:10]
colnames(dmrs_selected_strict_merge)


dmrs_selected_strict$region_size <- dmrs_selected_strict$region_end - dmrs_selected_strict$region_start + 1
dmrs_AA$region_size <- dmrs_AA$region_end - dmrs_AA$region_start + 1
###########
# Define function to plot p-values
#now plot manhattan plot for M1, M2 and M3
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(gtools)  # for mixedsort
plot_manhattan <- function(data, pval_col, title) {
  ggplot(data, aes(x = pos_cum, y = -log10(get(pval_col)), color = chr)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_color_manual(values = rep(c("grey30", "steelblue"), length.out = length(unique(data$chr)))) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      title = title,
      x = "Genomic Position",
      y = expression(-log[10](p))
    )
}

# Plot and save
p1 <- plot_manhattan(dmrs_selected_strict, "pval_FDR", "Manhattan Plot - GAM model")

print(p1)



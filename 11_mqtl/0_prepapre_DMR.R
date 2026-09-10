#prepare DMR regions for mqtl results. 
#note:
#DMR file should contain DMR_region_id DMR_chr DMR_region_start DMR_region_end columns



PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_scr14 <- paste0(PATH_wk, "/scr/14_paper/")


PATH_mqtl_results <- paste0(PATH_wk, "results/15_revision/8_mqtl/DMR_input/")
library(data.table)
library(dplyr)
#collect the M1, M2, M3 results after rerun




dmr_gene_file <- fread(file = paste0(PATH_scr14, "results/9_table_2.csv"))

#rename column chr to DMR_chr
gam_mqtl <- dmr_gene_file %>%
  rename(DMR_chr = chr, DMR_region_start = region_start, DMR_region_end = region_end) %>%#add a name DMR_region_id column as 1:nrow(dmr_gene_file)
  mutate(DMR_region_id = 1:nrow(dmr_gene_file))

#save gam_mqtl to a csv file
fwrite(gam_mqtl, file = paste0(PATH_mqtl_results, "gam_dmr_input.csv"))

########################################################################
#supplementary table 3: somnibus results of the selected cpg sets and genes

somnibus_DMRs <- fread(paste0(PATH_wk,"/results/9_regional/16_somnibus_regions.csv"))

#remove column pval_adjusted and rename pval to pvalue_somnibus
somnibus_DMRs <- somnibus_DMRs %>% 
  dplyr::select(-pval_adjusted) %>%
  rename(pvalue_somnibus = pval)

somnibus_mqtl <- somnibus_DMRs %>%
  rename(DMR_chr = chr, DMR_region_start = region_start, DMR_region_end = region_end) %>%
  mutate(DMR_region_id = 1:nrow(somnibus_DMRs))

fwrite(somnibus_mqtl, file = paste0(PATH_mqtl_results, "somnibus_dmr_input.csv"))

#########
#supplementary table 4: BSmooth results of the selected cpg sets and genes
BSmooth_DMRs <- fread( file = paste0(PATH_scr14, "results/6_BSmooth_DMRs.csv"))

#remove columns: region_id  data_chunk_id n_samples n_case n_ctrl n_dmrs
BSmooth_DMRs <- BSmooth_DMRs %>%
  dplyr::select(-region_id, -data_chunk_id, -n_samples, -n_case, -n_ctrl, -n_dmrs)


#rename columns: chr to DMR_chr, region_start to DMR_region_start, region_end to DMR_region_end
BSmooth_mqtl <- BSmooth_DMRs %>%
  rename(DMR_chr = chr, DMR_region_start = region_start, DMR_region_end = region_end) %>%
  mutate(DMR_region_id = 1:nrow(BSmooth_DMRs))

fwrite(BSmooth_mqtl, file = paste0(PATH_mqtl_results, "BSmooth_dmr_input.csv"))
###############
#supplementary table 5: DMRcate results of the selected cpg sets and genes
DMRcate_DMRs <- fread( file = paste0(PATH_scr14, "results/7_DMRcate_DMRs.csv"))
#remove column: data_chunk_id n_samples n_case n_ctrl n_dmrs region_id   
DMRcate_DMRs <- DMRcate_DMRs %>%
  dplyr::select(-data_chunk_id, -n_samples, -n_case, -n_ctrl, -n_dmrs, -region_id)

#rename columns: chr to DMR_chr, region_start to DMR_region_start, region_end to DMR_region_end
DMRcate_mqtl <- DMRcate_DMRs %>%
  rename(DMR_chr = chr, DMR_region_start = region_start, DMR_region_end = region_end) %>%
  mutate(DMR_region_id = 1:nrow(DMRcate_DMRs))
fwrite(DMRcate_mqtl, file = paste0(PATH_mqtl_results, "DMRcate_dmr_input.csv"))
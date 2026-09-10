########################################################################
# Merge DMRs from all four regional methods
########################################################################
#prepare DMR regions for mqtl results. 
#note:
#DMR file should contain DMR_region_id DMR_chr DMR_region_start DMR_region_end columns



PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_scr14 <- paste0(PATH_wk, "/scr/14_paper/")


PATH_mqtl_results <- paste0(PATH_wk, "results/15_revision/8_mqtl/DMR_input/")
library(data.table)
library(dplyr)
#collect the M1, M2, M3 results after rerun

PATH_DMR_report_results <- paste0(PATH_wk, "results/15_revision/9_dmr_report/DMR_input/")

dir.create(PATH_DMR_report_results, recursive = TRUE, showWarnings = FALSE)

dmr_gene_file <- fread(file = paste0(PATH_scr14, "results/9_table_2.csv"))

#rename column chr to DMR_chr
gam_mqtl <- dmr_gene_file %>%
  rename(DMR_chr = chr, DMR_region_start = region_start, DMR_region_end = region_end) %>%#add a name DMR_region_id column as 1:nrow(dmr_gene_file)
  mutate(DMR_region_id = 1:nrow(dmr_gene_file))

#save gam_mqtl to a csv file

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



library(data.table)
library(dplyr)

# ----------------------------------------------------------------------
# 1. Add method name to each DMR table
# ----------------------------------------------------------------------

gam_all <- gam_mqtl %>%
  mutate(method = "GAM-DMR")

somnibus_all <- somnibus_mqtl %>%
  mutate(method = "SOMNiBUS")

BSmooth_all <- BSmooth_mqtl %>%
  mutate(method = "BSmooth")

DMRcate_all <- DMRcate_mqtl %>%
  mutate(method = "DMRcate")


# ----------------------------------------------------------------------
# 2. Stack all DMRs into one long-format table
#
# fill = TRUE is important because the four methods have some
# method-specific columns that are not identical.
# ----------------------------------------------------------------------

all_DMRs_long <- rbindlist(
  list(
    gam_all,
    somnibus_all,
    BSmooth_all,
    DMRcate_all
  ),
  use.names = TRUE,
  fill = TRUE
)

# Put the most useful columns first
setcolorder(
  all_DMRs_long,
  c(
    "method",
    "DMR_region_id",
    "DMR_chr",
    "DMR_region_start",
    "DMR_region_end",
    setdiff(
      names(all_DMRs_long),
      c(
        "method",
        "DMR_region_id",
        "DMR_chr",
        "DMR_region_start",
        "DMR_region_end"
      )
    )
  )
)

# Add genomic width
all_DMRs_long[
  ,
  DMR_width_bp := DMR_region_end - DMR_region_start + 1
]
#remove DMR_region_id column 
all_DMRs_long <- all_DMRs_long %>%
    dplyr::select(-DMR_region_id)

#only keep method, DMR_chr, DMR_region_start, DMR_region_end, DMR_width_bp columns
all_DMRs_long <- all_DMRs_long %>%
    dplyr::select(method, DMR_chr, DMR_region_start, DMR_region_end, DMR_width_bp)

# Save
fwrite(
  all_DMRs_long,
  file = paste0(PATH_DMR_report_results, "all_DMRs_long.csv")
)



########################################################################
# 4. Quick summaries
########################################################################

cat("\nNumber of DMRs by method:\n")
print(table(all_DMRs_long$method))

cat("\nTotal method-specific DMR rows:",
    nrow(all_DMRs_long), "\n")



##############
#
scp /mnt/c/Per/LaiJiang/Project/UQAC/meth/results/15_revision/9_dmr_report/DMR_input/all_DMRs_long.csv laj773@rorqual.calculquebec.ca:/home/laj773/scratch/UQAC/meth/results/15_revision/9_dmr_report/DMR_input/all_DMRs_long.csv

scp -r /mnt/c/Per/LaiJiang/Project/UQAC/meth/scr/15_revision/9_dmr_report/* \
laj773@rorqual.calculquebec.ca:/home/laj773/scratch/UQAC/meth/scr/15_revision/9_dmr_report/
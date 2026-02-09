#evaluate the cpg spacing startegy on a data chunk .

i_chunk_id <- 227 #testing data chunk id
###########################################################################
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_results <- paste0(PATH_wk,"results/10_sanity/")
PATH_scr6 <- paste0(PATH_wk,"scr/6_beluga/")
PATH_scr8 <- paste0(PATH_wk,"scr/8_rerun/")
PATH_scr9 <- paste0(PATH_wk, "scr/9_regional/")
PATH_scr10 <- paste0(PATH_wk, "scr/10_sanity/")
PATH_scr11 <- paste0(PATH_wk, "scr/11_mgcv/")

#PATH_scr10 <- "/mnt/c/Per/LaiJiang/Project/UQAC/meth/scr/10_sanity/"


library(dplyr)
Sys.sleep(0.1)

library(data.table)
Sys.sleep(0.1)

library(stringr)
Sys.sleep(0.1)

library(ggplot2)
Sys.sleep(0.1)

library(tidyr)
Sys.sleep(0.1)

library(mgcv)
Sys.sleep(0.1)



#now load that data chunk
#i_chunk <- fread(paste0(PATH_wk, "results/10_sanity/chunk_", sprintf("%04d", i_chunk_id), ".csv"))
i_chunk <- fread(paste0(PATH_wk, "dat/chunk_", sprintf("%04d", i_chunk_id), ".csv"))

i_cpg_info <- i_chunk %>%
  select(chr, start) %>%
  distinct() %>%
  rename(position = start)
#######################
#######################
#load the function of split data by cpg spacing
source(paste0(PATH_scr11, "0_split_spacing.R"))

cpg_regions <- splitCpGsBySpacing(i_cpg_info, gap = 250, min_cpgs = 10, max_cpgs = 2000) 

print(cpg_regions[1,])

fwrite(cpg_regions, paste0(PATH_wk, "dat/region_spacing/spacing_", sprintf("%04d", i_chunk_id), ".csv"), col.names = FALSE)

if(FALSE){
#visulaize the region spit

library(ggplot2)

# Plot each region as a horizontal line segment
ggplot1 <- ggplot(cpg_regions, aes(x = region_start, xend = region_end, y = chr, yend = chr)) +
  geom_segment(linewidth = 2, color = "steelblue") +
  theme_minimal() +
  labs(
    title = "CpG Regions by Genomic Coordinates",
    x = "Genomic Position (bp)",
    y = "Chromosome"
  ) +
  theme(panel.grid.major.y = element_blank())

#save ggplot
ggsave(filename = paste0(PATH_wk, "results/11_mgcv/2_cpg_regions_by_spacing.jpg"), plot = ggplot1, width = 10, height = 5)

}
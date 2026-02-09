#now this start over, with a given gene region, gene start and end, and methylation data on that region.

#we split the gene region into sub-regions baed on cpg spacing, and run the GAM model on each sub-region to test if the interaction between AA and methylation is significant.


#start with the example in 2_col.R:

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

#look at the first gene: PTPRE
i_gene <- 1

cpg_spacing <- 250 #in bp

i_results_dat <- NULL
###########################################################################
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_results <- paste0(PATH_wk,"results/10_sanity/")
PATH_scr6 <- paste0(PATH_wk,"scr/6_beluga/")
PATH_scr8 <- paste0(PATH_wk,"scr/8_rerun/")
PATH_scr9 <- paste0(PATH_wk, "scr/9_regional/")
PATH_scr10 <- paste0(PATH_wk, "scr/10_sanity/")


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

#load top 10 genes of well-established regions
gene_table <- read.csv(file=paste0(PATH_scr10,"dat/asthma_gene_ranking_with_gene_ranges.csv"), header=TRUE)

i_region <- gene_table[i_gene,]
i_position <- i_region$SNP.Position..hg19.
i_gene_range <- i_region$Full_Gene_Range_hg19

#extract the chromosome and position
chr_pos <- strsplit(i_position, ":")[[1]]; chr_pos <- gsub(",", "", chr_pos); chr_pos <- gsub("chr", "", chr_pos) #remove the comma


clean_range <- str_replace_all(i_gene_range, ",", "")
start_pos <- as.numeric(str_extract(clean_range, "(?<=:)[0-9]+"))
end_pos <- as.numeric(str_extract(clean_range, "(?<=–)[0-9]+"))


#the summar table of the region of interest
i_position <- data.frame(
    Chromosome = chr_pos[1],
    Position = as.numeric(chr_pos[2]),
    Gene = i_region$Gene.s.,
    GWAS_pval = i_region$GWAS.p.value,
    gene_start = start_pos,
    gene_end = end_pos
)

#now look for the data chunk that contains this region ?
#1. load the dat chunk first row info
##########################################################################################
#1737 data chunks
chunk_info <- read.table(file=paste0(PATH_wk,"results/10_sanity/first_row_chunk.txt"))

colnames(chunk_info) <- c("Chromosome", "Position")
chunk_info$chunk_id <- 1:nrow(chunk_info)
#2. find the chunk that contains this region,i.e. the top row that the chromosome and position match

func_chunk_id <- function(chr_input, pos_input){
i_chunk_id <- chunk_info %>%
  filter(Chromosome == chr_input & Position < pos_input) %>%
  pull(chunk_id) %>%
  max()
i_chunk_id
}


################################
# Add phenotype info (AA status)
#load the pheno_file from the created data
load(file=paste0(PATH_scr6,"4_1_data.RData"), verbose=TRUE)


#head(pheno_file)

i_pheno_file = pheno_file %>% rename(sample=ID) 


##########################################################################################
#which data chunk to load?
print("the required data chunk for GWAS SNP is:")
print(paste0("chunk_",func_chunk_id(i_position$Chromosome , i_position$Position) ))


print("the required data chunk for gene start is:")
print(paste0("chunk_",func_chunk_id(i_position$Chromosome , i_position$gene_start) ))

print("the required data chunk for gene end is:")
print(paste0("chunk_",func_chunk_id(i_position$Chromosome , i_position$gene_end) ))
##########################################################################################

#!!!!!!!
#which region to loook at based on the gene range?
i_chunk_id <- func_chunk_id(i_position$Chromosome , i_position$Position) 
#it seems the GWA SNP is 100 kb after the gene start and 100 kb before the gene end
#so we will show the whole span of the gene 

#now load that data chunk
i_chunk <- fread(paste0(PATH_wk, "results/10_sanity/chunk_", sprintf("%04d", i_chunk_id), ".csv"))

#######################
#first extract the region of interest in the data chunk
i_region_data <- i_chunk %>%
  filter(chr == i_position$Chromosome & 
         start >= i_position$gene_start & 
         start <= i_position$gene_end)
##################################################################
#split the region into sub-regions based on cpg spacing
#assume the cpg spacing is 1000 bp, and we want to split the region into sub-regions of 1000 bp each


library(data.table)

split_region_by_distance <- function(i_region_data, gap_threshold = 250) {
  # Ensure data.table format
  dt <- as.data.table(i_region_data)
  
  # Order by genomic position
  dt <- dt[order(chr, start)]

  # Compute gap between consecutive CpGs
  dt[, gap := c(0, diff(start))]

  # Mark new subregion when gap > threshold
  dt[, subregion := cumsum(gap > gap_threshold)]

  # Collapse to one row per subregion: get start and end
  subregion_ranges <- dt[, .(
    subregion_start = min(start),
    subregion_end = max(end)
  ), by = subregion]

  return(subregion_ranges[, .(subregion_start, subregion_end)])
}

subregions <- split_region_by_distance(i_region_data, gap_threshold = cpg_spacing)
##################################################################
#now look at the 1st sub-region 
j_region_index <- 1

print("total subregions:")
print(nrow(subregions))

for(j_region_index in 1:nrow(subregions)){

print(paste0("processing subregion ", j_region_index, " of ", nrow(subregions)))

j_subregion <- subregions[j_region_index,]

# Reshape to long format if needed
j_subregion_data <- i_region_data %>%
  filter(chr == i_position$Chromosome & 
         start >= j_subregion$subregion_start & 
         start <= j_subregion$subregion_end)

if( length(unique(j_subregion_data$start))>20){
meth_long <- j_subregion_data %>%
  select(start, matches("_meth$")) %>%
  pivot_longer(-start, names_to = "sample", values_to = "meth") %>% #remove the NA rows
  na.omit() %>% #remove the "_meth" from the sample 
  mutate(sample = str_replace(sample, "_meth", "")) #remove the "_meth" from the sample names


#now attach the i_pheno_file to the i_meth_long
#library(dplyr)

i_meth_long_cov <- meth_long %>%
  left_join(i_pheno_file, by = "sample") %>%
  na.omit() # remove rows with NA values


  


# Optional: transform methylation !!! try different transofrmations later !!!!
i_meth_long_cov$meth_arcsin <- asin(sqrt(i_meth_long_cov$meth))


#fix the column type bug
for (v in colnames(i_meth_long_cov)) {
  # If it's matrix/array in training data, convert to numeric
  if (is.matrix(i_meth_long_cov[[v]]) || is.array(i_meth_long_cov[[v]])) {
    i_meth_long_cov[[v]] <- as.numeric(i_meth_long_cov[[v]])
  }
}



gam_model <- mgcv::gam(meth_arcsin ~ 
                   s(start, bs = "cs") +                   # smooth baseline
                   s(start, by = AA_only, bs = "cs") +     # group-specific smooth
                   AA_only +                               # fixed group effect
                   AgeCalc + Sex + Non.smoker +            # fixed covariates
                   EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5 +
                   s(FID, bs = "re"),                      # random intercept for FID
                 data = i_meth_long_cov,
                 method = "REML")

#interpretation in the documents
# Extract GAM summary
gam_summary <- summary(gam_model)
# Get smooth term table
smooth_table <- gam_summary$s.table
# Extract p-value for s(start):AA_only
pval_smooth_AA <- smooth_table["s(start):AA_only", "p-value"]
# View

j_results_dat <- data.frame(
    gene_name = i_position$Gene,
    N_subregions = nrow(subregions),
    subregion_index = j_region_index,
    subregion_chr= j_subregion_data$chr[1],
    subregion_start = j_subregion$subregion_start,
    subregion_end = j_subregion$subregion_end,
    pval_smooth_AA = pval_smooth_AA,
    stringsAsFactors = FALSE
)

i_results_dat <- rbind(i_results_dat, j_results_dat)

}

}

#save the results
write.csv(i_results_dat, file=paste0(PATH_results,"3_all_gam_subregion_results.csv"), row.names = FALSE)

#i_results_dat <- read.csv(file=paste0(PATH_results,"3_all_gam_subregion_results.csv"), header=TRUE)
#test run somnibus on a single region
#myabe also bump hunter? https://bioconductor.org/packages/release/bioc/html/bumphunter.html

#following such results: https://www.bioconductor.org/packages/release/bioc/vignettes/SOMNiBUS/inst/doc/SOMNiBUS.html


library(dplyr)
library(SOMNiBUS)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No command-line arguments supplied. Please pass the SLURM_ARRAY_TASK_ID.")
}


if(FALSE){
# Convert the first argument to an integer
task_id <- as.integer(args[1])
cat("SLURM_ARRAY_TASK_ID is:", task_id, "\n")
}

#first data chunk
task_id <- 690
###########################################################################
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_scr6 <- paste0(PATH_wk,"scr/6_beluga/")
PATH_scr8 <- paste0(PATH_wk,"scr/8_rerun/")
PATH_scr9 <- paste0(PATH_wk, "scr/9_regional/")

#first to load one splited data C:\Per\LaiJiang\Project\UQAC\meth\dat\chunk_017.csv
meth_file_loc <- paste0(PATH_wk,"scr/6_beluga/chunk_",sprintf("%04d",task_id),".csv")

#location to save reults
results_loc <-   paste0(PATH_scr9,"results/AA_",sprintf("%04d",task_id))

###########################################################################
meth_file <- read.csv(file=meth_file_loc, header=TRUE, sep="\t")

#methfile contains the JobID!!!
#the saving lcoation also contains the JobID!!!

#load the pheno_file from the created data: load pheno_file
load(file=paste0(PATH_scr6,"4_1_data.RData"), verbose=TRUE)


#step1: convert the data to a long format
source(paste0(PATH_scr9,"0_func_somnibus_input.R"))
#
somnibus_input_meth <- convert_to_somnibus_format(meth_file)
rm(meth_file)

#step 2: attach the covaraites onto the somnibus input
# Filter datasets to only keep common IDs
common_ids <- intersect(somnibus_input_meth$ID, pheno_file$ID)

somnibus_input_filtered <- somnibus_input_meth %>%
  filter(ID %in% common_ids) %>%
  left_join(pheno_file, by = "ID") %>%
  select(-FID) %>%
  mutate(
    Meth_Counts = round(Meth_Counts * Total_Counts),
    Meth_Counts = as.numeric(Meth_Counts),
    Total_Counts = as.numeric(Total_Counts)
  )

# Check the result
#head(somnibus_input_filtered)

#filtering missing cpgs
somnibus_input <- na.omit(somnibus_input_filtered[somnibus_input_filtered$Total_Counts != 0, ])

rm(somnibus_input_filtered)

#step: run Somnibus analysis
#we use the default region split method:  which splits methylation data into regions based on the spacing of CpGs.
outs <- runSOMNiBUS(dat = somnibus_input,
                    split = list(approach = "region", gap = 500),
                    n.k = rep(10,14), p0 = 0.003, p1 = 0.9, 
                    min.cpgs = 10, max.cpgs = 2000, verbose = TRUE)

###############################################################################
#################################!!!!!!

###########################################

###########################################

#attach this line to a file
save( outs, file=paste0(results_loc,"_somnbibus.RData"))

    

######visulaization gaps


# Get unique CpG positions
cpg_positions <- somnibus_input %>%
  distinct(Position) %>%
  arrange(Position)

# Calculate spacing between consecutive CpGs
cpg_positions <- cpg_positions %>%
  mutate(gap = c(NA, diff(Position)))

library(ggplot2)


#save the plot as png
# Visualize CpG spacing


# Visualize CpG spacing: the choice of gap == 500 is from here:
p<-ggplot(cpg_positions, aes(x = gap)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  labs(title = "Histogram of CpG Spacing", x = "Gap between CpGs (bp, log scale)", y = "Frequency")

# Save the plot
ggsave(p, filename = paste0(PATH_scr9,"results/cpg_spacing.png"), width = 8, height = 6)

  
# Example for gap = 1e3 (1,000 bp)
regions_1k <- splitDataByRegion(dat = somnibus_input, gap = 100, min.cpgs = 10, max.cpgs = 2000, verbose = FALSE)
length(regions_1k)  # Number of regions

# Example for gap = 1e4 (10,000 bp)
regions_10k <- splitDataByRegion(dat = somnibus_input, gap = 500, min.cpgs = 10, max.cpgs = 2000, verbose = FALSE)
length(regions_10k)

#Check number of CpGs per region
summary(sapply(regions_1k, nrow))    # For gap = 1k
summary(sapply(regions_10k, nrow))   # For gap = 10k

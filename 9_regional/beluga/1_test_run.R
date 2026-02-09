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


if(TRUE){
# Convert the first argument to an integer
task_id <- as.integer(args[1])
cat("SLURM_ARRAY_TASK_ID is:", task_id, "\n")
}

#first data chunk!!!
#task_id <- 690
###########################################################################

PATH_wk <- "~/scratch/UQAC/meth/"

PATH_scr6 <- paste0(PATH_wk,"scr/6_whole/")

PATH_scr8 <- paste0(PATH_wk,"scr/8_rerun/")
PATH_scr9 <- paste0(PATH_wk, "scr/9_regional/")

#first to load one splited data C:\Per\LaiJiang\Project\UQAC\meth\dat\chunk_017.csv
meth_file_loc <- paste0(PATH_wk,"data/meth_split/chunk_",sprintf("%04d",task_id),".csv")

#location to save reults
results_loc <-   paste0(PATH_wk,"results/somnibus/AA_",sprintf("%04d",task_id))


###########################################################################
meth_file <- read.csv(file=meth_file_loc, header=TRUE, sep="\t")

#methfile contains the JobID!!!
#the saving lcoation also contains the JobID!!!

#load the pheno_file from the created data
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
                    split = list(approach = "region", gap = 5e4),
                    n.k =  rep(10,14), p0 = 0.003, p1 = 0.9, 
                    min.cpgs = 50, max.cpgs = 2000, verbose = TRUE)

###############################################################################
#################################!!!!!!

###########################################

###########################################

#attach this line to a file
save( outs, file=paste0(results_loc,"_somnbibus.RData"))

    
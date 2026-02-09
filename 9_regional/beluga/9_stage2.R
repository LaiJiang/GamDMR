#this file select the regions of interst 
library(data.table)
library(dplyr)
library(stringr)
library(SOMNiBUS)

PATH_wk <- "~/scratch/UQAC/meth/"


PATH_scr6 <- paste0(PATH_wk,"scr/6_whole/")

PATH_scr8 <- paste0(PATH_wk,"scr/8_rerun/")
PATH_scr9 <- paste0(PATH_wk, "scr/9_regional/")

#load the selected cpgs from the previous step
sel_cpgs <- read.table(
  file = paste0(PATH_scr9, "/data/4_sel_region.txt")
)


# extract only chr and start from the selected cpgs
sel_cpgs_chr_start <- gsub("-.*$", "", sel_cpgs[[1]])

#fix the selected covaraites from the 5_sel_cov.R

#conlcudion: first stage we only need to include age, response, and Non-smoker, and maybe EOSINOpc, LYMPHOpc.

#next: for each chunk, filter the selected cpgs (4_1.R) in the data chunk. THen use 25kb/50kb window to select cpg regions. Then rerun the SOMNIBUS only on the selected region of cpgs.

# List of features to check
ft_names <- c("Intercept","meth_AA", "AgeCalc", "Sex", "Non.smoker", 
              "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc", 
              "sv1", "sv2", "sv3", "sv4", "sv5")

#sel_covariate <- c("AgeCalc", "meth_AA","Non.smoker", "EOSINOpc", "LYMPHOpc")
#sel_covariate <- c("AgeCalc", "AA_only", "Non.smoker", "EOSINOpc", "LYMPHOpc")
sel_covariate <- c( "AA_only")

######################################################################################
######################################################################################
######################################################################################
#now load a data chunk


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

#first data chunk
#task_id <- 690
###########################################################################


#first to load one splited data C:\Per\LaiJiang\Project\UQAC\meth\dat\chunk_017.csv
meth_file_loc <- paste0(PATH_wk,"data/meth_split/chunk_",sprintf("%04d",task_id),".csv")

#location to save reults
results_loc <-   paste0(PATH_wk,"results/somnibus_stage2/AA_",sprintf("%04d",task_id))

###########################################################################
meth_file <- read.csv(file=meth_file_loc, header=TRUE, sep="\t")

#methfile contains the JobID!!!
#the saving lcoation also contains the JobID!!!

#load the pheno_file from the created data: load pheno_file
load(file=paste0(PATH_scr6,"4_1_data.RData"), verbose=TRUE)

chunk_cpg_ids <- paste0(meth_file$chr, ":", meth_file$start)

overlap_cpgs <- intersect(sel_cpgs_chr_start, chunk_cpg_ids)




######################################################################
#if no cpgs worth looking at, then skip the chunk
if(length(overlap_cpgs) ==0 ){

    print("No CpGs pass the univariate filtering step in this chunk.")
}



######################################################################
#if there are cpgs worth looking at, then proceed to the next step
if(length(overlap_cpgs) >0 ){


#first find regions of  interest
source(paste0(PATH_scr9,"0_func_region_6.R"))

#then create long format of the somnibus data input 


#step1: convert the data to a long format
source(paste0(PATH_scr9,"0_func_somnibus_input.R"))
#
somnibus_input_meth <- convert_to_somnibus_format(meth_file)
rm(meth_file)

#step 2: attach the covaraites onto the somnibus input
# Filter datasets to only keep common IDs

#note we only want the selected covariates at stage 1
pheno_file <- pheno_file %>%
  select(ID, FID, all_of(sel_covariate))



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

rm(somnibus_input_meth)
# Check the result
#head(somnibus_input_filtered)

#filtering missing cpgs
somnibus_input <- na.omit(somnibus_input_filtered[somnibus_input_filtered$Total_Counts != 0, ])

rm(somnibus_input_filtered)

#step: run Somnibus analysis
#we use the default region split method:  which splits methylation data into regions based on the spacing of CpGs.
#estimate n.k by the default method
n_k_dim <- max(3, as.integer(length(unique(somnibus_input$Position)) / 20))


outs <- runSOMNiBUS(dat = somnibus_input,
                    split = list(approach = "region", gap = 250),
                    n.k = rep(n_k_dim,length(sel_covariate)+1),
                     p0 = 0.003, p1 = 0.9, 
                    min.cpgs = 51, max.cpgs = 2000, verbose = TRUE)

###############################################################################
#################################!!!!!!

###########################################

###########################################

#attach this line to a file
save( outs, file=paste0(results_loc,"_somnibus.RData"))

    



}



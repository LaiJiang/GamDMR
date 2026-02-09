
#look at the first gene: PTPRE
i_gene <- 11

cpg_spacing <- 250 #in bp

i_results_dat <- NULL
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

#load top 10 genes of well-established regions
#gene_table <- read.csv(file=paste0(PATH_scr10,"dat/5_col_asthma_gene_table.csv"), header=TRUE) #old one

#updated gene table from script below
gene_table <- read.csv(file=paste0(PATH_scr11,"dat/asthma_gene_ranking_with_gene_ranges.csv"), header=TRUE)


###################################################################################################
#lets update the gene table with gene coordinates from Ensembl GRCh37
#only need to run once.
if(FALSE){

library(biomaRt)

# Connect to Ensembl GRCh37
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
                      host = "grch37.ensembl.org")

# Query for genes
genes <- gene_table$Gene.s.
attrs <- c("external_gene_name","chromosome_name","start_position","end_position")
gene_coords <- getBM(attributes = attrs, filters = "external_gene_name",
                     values = genes, mart = ensembl)
#attach the chromosome name to the gene table
gene_table <- gene_table %>%
  left_join(gene_coords, by = c("Gene.s." = "external_gene_name"))
#assign new range:
# Full_Gene_Range_hg19 = paste0("chr",chromosome_name,":",
#                               start_position,"–",end_position)
gene_table$Full_Gene_Range_hg19 <- paste0("chr", gene_table$chromosome_name, ":",
                                           gene_table$start_position, "–", gene_table$end_position)


#save gene_table with updated gene coordinates !!!
write.csv(gene_table, file=paste0(PATH_scr11,"dat/asthma_gene_ranking_with_gene_ranges.csv"), row.names = FALSE)

}
###################################################################################################

i_region <- gene_table[i_gene,]
i_position <- i_region$SNP.Position..hg19.
i_gene_range <- i_region$Full_Gene_Range_hg19

#extract the chromosome and position of GWAS SNP, but this one do not have one
#chr_pos <- strsplit(i_position, ":")[[1]]; chr_pos <- gsub(",", "", chr_pos); chr_pos <- gsub("chr", "", chr_pos) #remove the comma


clean_range <- str_replace_all(i_gene_range, ",", "")
start_pos <- as.numeric(str_extract(clean_range, "(?<=:)[0-9]+"))
end_pos <- as.numeric(str_extract(clean_range, "(?<=–)[0-9]+"))

#chromosome of the gene region
chr_pos <- as.numeric(gsub("chr","",strsplit(clean_range, ":")[[1]][1]))

#the summar table of the region of interest
i_position <- data.frame(
    Chromosome = chr_pos,
    #Position = as.numeric(chr_pos[2]),
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
#print("the required data chunk for GWAS SNP is:")
#print(paste0("chunk_",func_chunk_id(i_position$Chromosome , i_position$Position) ))


print("the required data chunk for gene start is:")
print(paste0("chunk_",func_chunk_id(i_position$Chromosome , i_position$gene_start) ))

print("the required data chunk for gene end is:")
print(paste0("chunk_",func_chunk_id(i_position$Chromosome , i_position$gene_end) ))
##########################################################################################


#!!!!!!!download data in local terminal

# scp laj773@beluga.computecanada.ca:/home/laj773/scratch/UQAC/meth/data/meth_split/chunk_0227.csv /mnt/c/Per/LaiJiang/Project/UQAC/meth/dat/


#which region to loook at based on the gene range?
i_chunk_id <- func_chunk_id(i_position$Chromosome , i_position$gene_start) 
#it seems the GWA SNP is 100 kb after the gene start and 100 kb before the gene end
#so we will show the whole span of the gene 

#now load that data chunk
#i_chunk <- fread(paste0(PATH_wk, "results/10_sanity/chunk_", sprintf("%04d", i_chunk_id), ".csv"))
i_chunk <- fread(paste0(PATH_wk, "dat/chunk_", sprintf("%04d", i_chunk_id), ".csv"))

#######################
#######################


#first extract the region of interest in the data chunk
i_region_data <- i_chunk %>%
  filter(chr == i_position$Chromosome & 
         start >= i_position$gene_start & 
         start <= i_position$gene_end)
###########################################################################
###########################################################################

#now look at the region of interest
# Reshape to long format if needed
library(dplyr)
library(tidyr)
library(stringr)

library(data.table)
library(stringr)

# Make sure your data is a data.table
setDT(i_region_data)

# Melt the data
meth_long <- melt(i_region_data,
                  id.vars = "start",
                  measure.vars = patterns("_meth$"),
                  variable.name = "sample",
                  value.name = "meth",
                  na.rm = TRUE)

# Clean the sample names
meth_long[, sample := str_replace(as.character(sample), "_meth", "")]

# Add phenotype info (AA status)
#load the pheno_file from the created data
load(file=paste0(PATH_scr6,"4_1_data.RData"), verbose=TRUE)
pheno_file <- as_tibble(pheno_file)

# Assuming you have pheno_file with AA column
meth_long <- meth_long %>%
  left_join(pheno_file %>% dplyr::select(ID, AA_only), by = c("sample" = "ID")) %>%
  na.omit()


# Quick plot
ggplot1 <- ggplot(meth_long, aes(x = start, y = meth, color = factor(AA_only))) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
  theme_minimal() + #add vertical line for i_position$Position
  geom_vline(xintercept = i_position$Position, linetype = "dashed", color = "red") +
  labs(title = paste("Methylation in region:", i_position$Gene),
       x = "Genomic Position",
       y = "Methylation Level",
       color = "AA Status") 


# Save the plot
ggsave(filename = paste0(PATH_wk, "results/11_mgcv/1_cor_ggplot1.png"), plot = ggplot1, width = 10, height = 6)


#save meth_long
save(i_position, meth_long, file = paste0(PATH_wk, "results/10_sanity/6_fp_1", i_gene, ".RData"))


#####################################2_col.R replication on FP gene###################

#after 1_test.R get the raw methylation data.
#we can see that there are some delicate signals in the localized region

#why some univaraite models are not significant
#why somnibus is not signficiant 

#step1: first verify if any of these cpgs in the region are significant with raw group comparision



#save meth_long

# Load required package
library(dplyr)

# Function to run Wilcoxon test for each CpG site
cpg_test_results <- meth_long %>%
  group_by(start) %>%
  summarise(
    pvalue = tryCatch(
      wilcox.test(meth ~ AA_only)$p.value,
      error = function(e) NA_real_  # handle any errors due to low sample size
    ),
    nsample = n()
  ) %>%
  ungroup()

# View top results
head(cpg_test_results)



library(ggplot2)

# Add -log10(pvalue) column
cpg_test_results$logp <- -log10(cpg_test_results$pvalue)

#sopped here.
# Plot
ggplot2<- ggplot(cpg_test_results, aes(x = start, y = logp)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  theme_minimal() +
  labs(
    title = paste0("CpG Association P-values in Region: ", i_position$Gene),
    x = "Genomic Position (start)",
    y = expression(-log[10](p))
  ) + #add vertical line for i_position$Position
  geom_vline(xintercept = i_position$Position, linetype = "dashed", color = "red") 

# Save the plot
ggsave(filename = paste0(PATH_wk, "results/11_mgcv/1_corr_ggplot2.jpg"), plot = ggplot2, width = 10, height = 6)


####################################################
# step2: how the somnibus results, do any of these regions even have results?

# Load somnibus results
#evalaute the region-specific results from 13_col.R


all_df <- read.csv( file = paste0(PATH_wk,"/results/9_regional/16_eval_dat.csv"))

#
i_df <- all_df %>% filter(chr==i_position$Chromosome) %>% #region_start < i_position$Position & region_end > i_position$Position)
  filter((region_start < i_position$gene_start & region_end > i_position$gene_start) | 
           (region_start < i_position$gene_end & region_end > i_position$gene_end) |
           (region_start > i_position$gene_start & region_end < i_position$gene_end))

print(i_df)

#somnibus ignored this region because there are very few cpgs in the region, and the cpgs are too sparse to be considered as a region.




#####################################################################

#for univarite analysis, do we find any signals in this region?

# Load the univariate results
#this file select the regions of interst 
library(data.table)
library(dplyr)
library(stringr)

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

#before rerun, the first stage results
AA_simple_results <- fread(paste0(PATH_wk,"/scr/8_rerun/dat/rerun_simple.txt"), header=FALSE)
colnames(AA_simple_results) <- c("Chunk_index","cpg_index","CpG",
                                  "M1_coef","M1_pval","M1_F_sdv",
                                  "M2_coef","M2_pval","M2_F_sdv",
                                  "M3_coef","M3_pval","M3_F_sdv" ,
                                  "M1_lambda","M3_lambda")                                
                                  
#remove dupliated rows
AA_simple_results <- AA_simple_results[!duplicated(AA_simple_results), ]

#replace the NA in M2_coef as 0 
AA_simple_results$M2_coef[is.na(AA_simple_results$M2_coef)] <- 0
AA_simple_results$M2_pval[is.na(AA_simple_results$M2_pval)] <- 1

min(AA_simple_results$M1_pval)
min(AA_simple_results$M3_pval)


AA_simple_results_chr <- AA_simple_results %>%
#cpg contain i_position$Chromosome and start between i_position$gene_start and i_position$gene_end
  filter(str_detect(CpG, paste0(i_position$Chromosome, ":")) )


rm(AA_simple_results)
#extract the cpgs are in the region of interest
AA_simple_results_chr <- AA_simple_results_chr %>%
  filter(str_extract(CpG, "(?<=:)[0-9]+") %>% as.numeric() >= i_position$gene_start & 
         str_extract(CpG, "(?<=:)[0-9]+") %>% as.numeric() <= i_position$gene_end)


#now plot manhattan plot for M1, M2 and M3
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(gtools)  # for mixedsort

# Convert to data.table
DT <- as.data.table(AA_simple_results_chr)

# Correctly extract chr, start, and end from CpG
DT[, c("chr", "start", "end") := tstrsplit(CpG, "[:-]", fixed = FALSE)]
DT[, start := as.numeric(start)]
DT[, chr := gsub("^chr", "", chr)]  # optional, in case "chr" prefix exists
DT[, chr := factor(chr, levels = mixedsort(unique(chr)))]

# Order by chromosome and position
setorder(DT, chr, start)

# Compute cumulative position for Manhattan plot
chr_lengths <- DT[, .(chr_len = max(start, na.rm = TRUE)), by = chr]
chr_lengths[, chr_start := cumsum(shift(chr_len, fill = 0))]

DT <- merge(DT, chr_lengths[, .(chr, chr_start)], by = "chr", all.x = TRUE)
DT[, pos_cum := start + chr_start]

# Define function to plot p-values
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
p1 <- plot_manhattan(DT, "M1_pval", "Manhattan Plot - Model 1")
p2 <- plot_manhattan(DT, "M2_pval", "Manhattan Plot - Model 2")
p3 <- plot_manhattan(DT, "M3_pval", "Manhattan Plot - Model 3")

# Display plots
print(p1)
print(p2)
print(p3)

# Optionally save
ggsave(paste0(PATH_wk, "results/11_mgcv/manhattan_M1.jpg"), p1, width = 10, height = 5)
ggsave(paste0(PATH_wk, "results/11_mgcv/manhattan_M2.jpg"), p2, width = 10, height = 5)
ggsave(paste0(PATH_wk, "results/11_mgcv/manhattan_M3.png"), p3, width = 10, height = 5)
#M3 no signals at all all pvalues = 1
#M1 and M2 have some signals,but not significant min-pvalue = 0.001
print(AA_simple_results_chr[which.min(AA_simple_results_chr$M1_pval),])



print(AA_simple_results_chr[which.min(AA_simple_results_chr$M2_pval),])


#would that possible to poll that weak signals into a strong regional signal, with mgcv method ?
#########################################################################################################
#########################################################################################################
#########################################################################################################

#now attaches covariates information to the methylation data
PATH_scr6 <- paste0(PATH_wk,"scr/6_beluga/")

load(file=paste0(PATH_scr6,"4_1_data.RData"), verbose=TRUE)

head(pheno_file)

i_pheno_file = pheno_file %>% rename(sample=ID) %>% 
  select(-AA_only) 

#now attach the i_pheno_file to the i_meth_long
library(dplyr)

#note: we first run i_meth_long = meth_long for all region. did not find signals
#we then zoom into the last part of region.



i_meth_long <- meth_long[which(meth_long$start > 119290000), ] #select the first 8 cpgs for testing


# 1. Sort by start position
sorted_starts <- i_meth_long_cov %>%
  arrange(start) %>%
  pull(start)
cpg_spacing <- diff(sorted_starts)

#the cpg signal is at 119303306, so it is between the range (119299905,XXX)

i_meth_long <- meth_long[which(meth_long$start > 119299905), ] #select the first 8 cpgs for testing

i_meth_long_cov <- i_meth_long %>%
  left_join(i_pheno_file, by = "sample") %>%
  na.omit() # remove rows with NA values


  

library(mgcv)

# Optional: transform methylation !!! try different transofrmations later
i_meth_long_cov$meth_arcsin <- asin(sqrt(i_meth_long_cov$meth))


#fix the column type bug
for (v in colnames(i_meth_long_cov)) {
  # If it's matrix/array in training data, convert to numeric
  if (is.matrix(i_meth_long_cov[[v]]) || is.array(i_meth_long_cov[[v]])) {
    i_meth_long_cov[[v]] <- as.numeric(i_meth_long_cov[[v]])
  }
}


k_cpg <- max(5, floor(length(unique(i_meth_long_cov$start)) /20))
k_cpg <- min(k_cpg, 10) # cap at 10 to avoid overfitting
gam_model <- mgcv::gam(meth_arcsin ~ 
                   s(start, bs = "cs",k=k_cpg) +                   # smooth baseline
                   s(start, by = AA_only, bs = "cs", k = k_cpg) +     # group-specific smooth
                   AA_only +                               # fixed group effect
                   AgeCalc + Sex + Non.smoker +            # fixed covariates
                   EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5 +
                   s(FID, bs = "re"),                      # random intercept for FID
                 data = i_meth_long_cov,
                 method = "REML")

#interpretation in the documents

#now plot the ftitted curves by AA status: 
source(paste0(PATH_scr10,"/0_func_plot_gam.R"))
ggsave(filename = paste0(PATH_wk, "results/11_mgcv/1_col_gam_fitted_curves.jpg"), plot=ggplot_gam,width = 10, height = 5)  


#now plot the difference between AA=1 and AA=0 curves
source(paste0(PATH_scr10,"/0_gam_region.R"))
ggsave(filename = paste0(PATH_wk, "results/11_mgcv/1_col_gam_region_diff.jpg"), plot=ggplot_gam_region, width = 10, height = 5)

####################################################################################

#vconcluion: the regional results are not significant, which shows the mgcv model thinks that the single cpg signal (from model 1 and model 2)
#below is not strong enough to support a regional signal. !!!! mgcv pval = 0.03.



#         Chunk_index cpg_index                    CpG  M1_coef      M1_pval
#1: AA_0227_simple.txt      1728 10:119303306-119303307 5.533025 1.447412e-06
#    M1_F_sdv    M2_coef      M2_pval   M2_F_sdv M3_coef M3_pval  M3_F_sdv
#       <num>      <num>        <num>      <num>   <num>   <num>     <num>
#1: 0.2128317 0.06373638 6.425879e-07 0.02345242       0       1 0.1365992
#   M1_lambda M3_lambda
#       <int>     <int>
#1:        10        10



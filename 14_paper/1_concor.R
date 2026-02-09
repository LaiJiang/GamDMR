
#prepare materials for the Nov_17_2025_regional_1.doc file on google drive

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
load( file = paste0(PATH_wk,"/scr/12_cpg_DMR/results/5_eval.RData"),verbose=TRUE)

sum(M2_manhattan$pval == min(M2_manhattan$pval, na.rm=TRUE), na.rm=TRUE)

#the cpg of interest that has the minimum pvalue in M2 model is:
cpg_int <- M2_manhattan[which(M2_manhattan$pval == min(M2_manhattan$pval, na.rm=TRUE)),]



PATH_scr6 <- paste0(PATH_wk,"scr/6_beluga/")
PATH_scr8 <- paste0(PATH_wk,"scr/8_rerun/")
PATH_scr9 <- paste0(PATH_wk, "scr/9_regional/")
PATH_scr10 <- paste0(PATH_wk, "scr/10_sanity/")


library(dplyr)
Sys.sleep(0.1)

library(data.table)
Sys.sleep(0.1)

library(stringr)

chr_pos <- strsplit(cpg_int$CpG, ":")[[1]][1]
start_pos <- as.numeric(str_extract(cpg_int$CpG, "(?<=:)[0-9]+"))

#note end_pos can be NA sometines
end_pos <- as.numeric(str_extract(cpg_int$CpG, "(?<=–)[0-9]+"))

#the summar table of the region of interest
i_position <- data.frame(
    Chromosome = as.numeric(chr_pos),
    Position = as.numeric(start_pos),
    GWAS_pval = cpg_int$pval
)

#now look for the data chunk that contains this region ?
#1. load the dat chunk first row info

#the gene CNOT3  contain the region:  chr19:54,137,749-54,155,681

i_gene <- "CNOT3"
i_position$gene_start <- 54137749
i_position$gene_end <- 54155681
i_position$Gene <- i_gene
#we manually set the gene end to curate a proper region rather than two separate regions.
i_position$gene_end <- 54151900
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
##########################################################################################
#which data chunk to load?
print("the required data chunk for GWAS SNP is:")
print(paste0("chunk_",func_chunk_id(i_position$Chromosome , i_position$Position) ))

##########################################################################################

#!!!!!!!
#which region to loook at based on the gene range?
i_chunk_id <- func_chunk_id(i_position$Chromosome , i_position$Position) 
#it seems the GWA SNP is 100 kb after the gene start and 100 kb before the gene end
#so we will show the whole span of the gene 

#now load that data chunk !!!!!!
i_chunk <- fread(paste0(PATH_wk, "dat/chunk_", sprintf("%04d", i_chunk_id), ".csv"))
#######################
#first extract the region of interest in the data chunk thats overlapping with this gene

i_region_data <- i_chunk %>%
  filter(chr == i_position$Chromosome & 
         start >= i_position$gene_start & 
         start <= i_position$gene_end)

#now look at the region of interest
# Reshape to long format if needed
library(ggplot2)
library(tidyr)

meth_long <- i_region_data %>%
  select(start, matches("_meth$")) %>%
  pivot_longer(-start, names_to = "sample", values_to = "meth") %>% #remove the NA rows
  na.omit() %>% #remove the "_meth" from the sample 
  mutate(sample = str_replace(sample, "_meth", "")) #remove the "_meth" from the sample names



# Add phenotype info (AA status)
#load the pheno_file from the created data
load(file=paste0(PATH_scr6,"4_1_data.RData"), verbose=TRUE)


# Assuming you have pheno_file with AA column
meth_long <- meth_long %>%
  left_join(pheno_file %>% select(ID, AA_only), by = c("sample" = "ID")) %>% na.omit()



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

# Quick plot
ggplot1 <- ggplot(meth_long, aes(x = start, y = meth, color = factor(AA_only))) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
  theme_minimal() +
  # Arrow instead of vertical line
  geom_segment(aes(x = i_position$Position, 
                   xend = i_position$Position,
                   y = 0.4, 
                   yend = 0.6),
               arrow = arrow(length = unit(0.15, "cm")),
               color = "red",
               linewidth = 0.8) +
  labs(title = paste("Methylation in region:", i_position$Gene),
       x = "Genomic Position",
       y = "Methylation Level",
       color = "AA Status")


# Save the plot
ggsave(filename = paste0(PATH_wk, "results/14_paper/1_meth_region_", i_gene, ".jpeg"), plot = ggplot1, width = 10, height = 6)


#save meth_long
save(i_position, meth_long, file = paste0(PATH_wk, "results/14_paper/meth_dat_", i_gene, ".RData")) 


#check the data on the i_Position$Position
meth_at_snp <- meth_long %>%
  filter(start == i_position$Position)  

  #now plot boxplot of meth by AA_only at this position
ggplot2 <- ggplot(meth_at_snp, aes(x = factor(AA_only), y = meth, fill = factor(AA_only))) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = paste("Methylation at SNP position:", i_position$Position),
       x = "AA Status",
       y = "Methylation Level",
       fill = "AA Status")  

  ggsave(filename = paste0(PATH_wk, "results/14_paper/1_boxplot_", i_gene, ".jpeg"), plot = ggplot2, width = 10, height = 6)
     


#for meth_at_snp, how many AA_only ==1 and AA_only ==0
table(meth_at_snp$AA_only,meth_at_snp$meth == 1)

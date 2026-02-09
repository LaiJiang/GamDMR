#test run analyis on 1 region 

#first look up the data chunk of interest

#look at the first gene: PTPRE
i_gene <- 1
###########################################################################
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_scr6 <- paste0(PATH_wk,"scr/6_beluga/")
PATH_scr8 <- paste0(PATH_wk,"scr/8_rerun/")
PATH_scr9 <- paste0(PATH_wk, "scr/9_regional/")
PATH_scr10 <- paste0(PATH_wk, "scr/10_sanity/")


library(dplyr)
Sys.sleep(0.1)

library(data.table)
Sys.sleep(0.1)

library(stringr)

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


# Save the plot
ggsave(filename = paste0(PATH_wk, "results/10_sanity/meth_region_", i_gene, ".png"), plot = ggplot1, width = 10, height = 6)


#save meth_long
save(i_position, meth_long, file = paste0(PATH_wk, "results/10_sanity/meth_dat_", i_gene, ".RData"))

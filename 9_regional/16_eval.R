#collect the chr information and evaluate somnibus regions
library(dplyr)
library(stringr)

#evalaute the region-specific results from 13_col.R
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

all_df <- read.csv( file = paste0(PATH_wk,"/results/9_regional/14_eval_pvals.csv"))

#

#load the chrosomome information 
chr_info <- read.csv(file = paste0(PATH_wk,"/results/9_regional/15_extract_chr.csv"))

chr_info <- chr_info %>%
  rowwise() %>%
  mutate(
    chr_freq = names(which.max(table(strsplit(chr, ":")[[1]]))),
    chr_no = length(unique(strsplit(chr, ":")[[1]]))
  ) %>%
  ungroup()

#check if every region has a unique chromosome
unique(chr_info$chr_no)

#Yes! now remove the column: chr start chr_no
chr_info <- chr_info %>%
  select(region_id, chr_freq) %>% unique() %>% rename(chr = chr_freq)


#now attach chr_info to all_df
all_df <- all_df %>%
  left_join(chr_info, by = "region_id")

  write.csv(all_df, file = paste0(PATH_wk,"/results/9_regional/16_eval_dat.csv"), row.names = FALSE)

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
#table 1: X-axis: number of significant cpgs, y-axis: number of overlapping cpgs with somnibus rsults

all_df <- read.csv( file = paste0(PATH_wk,"/results/9_regional/16_eval_dat.csv"))


library(data.table)

load( file = paste0(PATH_wk,"/scr/8_rerun/results/4_eval.RData"),verbose=TRUE)



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



# Extract significant CpG indices for each model:
sig_model1 <- M1_manhattan$CpG[which(M1_manhattan$pval<1e-5)]
sig_model2 <- AA_simple_results$CpG[which(AA_simple_results$M2_pval<1e-5)]
sig_model3 <- M3_manhattan$CpG[which(M3_manhattan$pval<1e-5)]

####################################################################################
####################################################################################
####################################################################################
#now generate the list of significant cpg selections from each model:

Model1_res <- M1_manhattan %>% select(CpG,pval) %>%
  filter(pval < 1e-5) %>%
  arrange(pval) %>% mutate(model = "Model1") 


Model2_res <- AA_simple_results %>% filter(M2_pval < 1e-5) %>%
  select(CpG, M2_pval) %>%
  arrange(M2_pval) %>% rename(pval=M2_pval)%>% mutate(model = "Model2") 

Model3_res <- M3_manhattan %>% select(CpG,pval) %>%
  filter(pval < 1e-5) %>%
  arrange(pval) %>% mutate(model = "Model3") 


somnibus_res <- all_df %>%
  filter(pval_adjusted < 1e-5) %>%
  select( chr, region_start, region_end, pval_adjusted) %>%
  mutate(model = "Somnibus") %>% mutate(CpG = paste0(chr, ":", region_start, "-", region_end)) %>%
  mutate(pval = pval_adjusted) %>%
  select(CpG, pval, model) %>%
  arrange(pval) 

#now combine the results from all models
all_models_res <- bind_rows(Model1_res, Model2_res, Model3_res, somnibus_res) 

#write to  a csv file
write.csv(all_models_res, file = paste0(PATH_wk,"/results/9_regional/16_eval_sig_cpgs_all_Models.csv"), row.names = FALSE)
####################################################################################
#now a supplementary for all cpgs in the models 


Model1_res <- M1_manhattan %>% select(CpG,pval) %>%
  arrange(pval) %>% mutate(model = "Model1") 


Model2_res <- AA_simple_results %>% 
  select(CpG, M2_pval) %>%
  arrange(M2_pval) %>% rename(pval=M2_pval)%>% mutate(model = "Model2") 

Model3_res <- M3_manhattan %>% select(CpG,pval) %>%
  arrange(pval) %>% mutate(model = "Model3") 


somnibus_res <- all_df %>%
  select( chr, region_start, region_end, pval_adjusted) %>%
  mutate(model = "Somnibus") %>% mutate(CpG = paste0(chr, ":", region_start, "-", region_end)) %>%
  mutate(pval = pval_adjusted) %>%
  select(CpG, pval, model) %>%
  arrange(pval) 

#now combine the results from all models
all_models_res <- bind_rows(Model1_res, Model2_res, Model3_res, somnibus_res) 

#write to  a csv file
write.csv(all_models_res, file = paste0(PATH_wk,"/results/9_regional/16_eval_all_cpgs_all_Models.csv"), row.names = FALSE)
####################################################################################
####################################################################################
######
#for each sig_model1, how many cpgs are overlapping with the all_df?

func_somnibus_cpg_overlap<- function(cpg_list){

# Step 1: Parse sig_model1 into chr and start
sig_model1_df <- data.frame(
  raw = cpg_list,
  stringsAsFactors = FALSE
) %>%
  mutate(
    chr = as.integer(str_extract(raw, "^[^:]+")),
    start = as.integer(str_extract(raw, "(?<=:)[0-9]+"))
  )

# Step 2: For each CpG, check if it overlaps any region in all_df
overlaps <- sig_model1_df %>%
  rowwise() %>%
  mutate(
    overlap = any(
      all_df$chr == chr &
      start >= all_df$region_start &
      start <= all_df$region_end
    )
  ) %>%
  ungroup()  %>%
  filter(overlap) %>%
  select(raw, chr, start, overlap)

# Step 3: Count number of overlaps
n_overlaps <- sum(overlaps$overlap)

return(overlaps)
}



model1_overlap <- func_somnibus_cpg_overlap(sig_model1)
model2_overlap <- func_somnibus_cpg_overlap(sig_model2)
model3_overlap <- func_somnibus_cpg_overlap(sig_model3)


print(intersect(model1_overlap$raw, model2_overlap$raw))


 all_df <-  read.csv( file = paste0(PATH_wk,"/results/9_regional/16_eval_dat.csv"))


print(sum(all_df$pval_adjusted<1e-5))


DMRs <- all_df %>%
  filter(pval_adjusted < 1e-5) %>%
  dplyr::select(region_id, chr, region_start, region_end)


  source("C:/Per/LaiJiang/Project/UQAC/meth/scr/9_regional/0_func_genes_16.R")


#select the DMRs with smallest 100 pval_adjusted
DMRs <- all_df %>%
  arrange(pval_adjusted) %>%
  mutate(rank = row_number()) 

DMRs <- DMRs %>%
  filter(rank <= 100) %>% 
    dplyr::select(region_id, chr, region_start, region_end)



  source("C:/Per/LaiJiang/Project/UQAC/meth/scr/9_regional/0_func_genes_16.R")


# How many unique genes?
length(unique(df_annot$genes_symbol))
####################################################################################


####################################################################################

all_cpgs<- AA_simple_results$CpG
rm(AA_simple_results)


library(GenomicRanges)
library(data.table)



# Prepare CpGs
# Split the all_cpgs string into chr, start, end
cpg_split <- tstrsplit(all_cpgs, "[:-]", type.convert = TRUE)
cpg_dt <- data.table(chr = as.integer(cpg_split[[1]]),
                     start = as.integer(cpg_split[[2]]),
                     end = as.integer(cpg_split[[3]]))

#restrict cpg_dt to the chr in somnibus_res_region$chr
cpg_dt <- cpg_dt[chr %in% unique(somnibus_res_region$chr)]    


# Prepare CpGs as GRanges
gr_cpgs <- GRanges(seqnames = paste0("chr", cpg_dt$chr),
                   ranges = IRanges(start = cpg_dt$start, end = cpg_dt$start+1))



# Prepare regions

#convert somnibus_res into chr start end format


somnibus_res_region <- all_df %>%
  filter(pval_adjusted < 0.01) %>%
  select( chr, region_start, region_end, pval_adjusted) %>%
  mutate(pval = pval_adjusted) %>%
  arrange(pval)


write.csv(somnibus_res_region, file = paste0(PATH_wk,"/results/9_regional/16_somnibus_regions.csv"), row.names = FALSE)


somnibus_res_region$chr <- as.integer(somnibus_res_region$chr)  # make sure it's integer
gr_regions <- GRanges(seqnames = paste0("chr", somnibus_res_region$chr),
                      ranges = IRanges(start = somnibus_res_region$region_start,
                                       end = somnibus_res_region$region_end))

hits <- findOverlaps(gr_cpgs, gr_regions)
overlapping_cpgs <- all_cpgs[queryHits(hits)]

# Output
length(overlapping_cpgs)  # how many overlaps found
head(overlapping_cpgs)    # example CpGs


somnibus_cpgs<- overlapping_cpgs

#remove all characters after the first "-"
somnibus_cpgs <- sub("-.*", "", somnibus_cpgs)

sig_model1 <- sub("-.*", "", sig_model1)
sig_model2 <- sub("-.*", "", sig_model2)
sig_model3 <- sub("-.*", "", sig_model3)


intersect(somnibus_cpgs, sig_model1)
intersect(somnibus_cpgs, sig_model2)
intersect(somnibus_cpgs, sig_model3)

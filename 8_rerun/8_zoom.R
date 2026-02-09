#this file zoom in to the HLA -MHC region, and check the pvalues for the cpgs in this region.



# 1. Load libraries
#library(clusterProfiler)
#library(org.Hs.eg.db)
#library(DOSE)  # sometimes needed for GSEA


PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(qqman)

library(VennDiagram)
####################################################################################
mhc_start <- 28477797
mhc_end   <- 33448354
####################################################################################
#Model 2 first 

#only run once: extract M2 cpgs 


##########################################
#temporary unmask for qqqplot !!!!!!!!!!!

print("read file")
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



#generate CpG_pval data frame and run gene annotations 
M2_manhattan <- AA_simple_results[,c("Chunk_index","CpG","M2_pval", "M2_coef" , "M2_F_sdv")]
CpG_pval <- M2_manhattan %>%  dplyr::rename(pval = `M2_pval`)

rm(M2_manhattan)


#!!!!!!!!!!!
if(TRUE){


#####################################################
#extract MHC region 

library(dplyr)
library(tidyr)

# Define the MHC interval (hg19 coordinates; adjust if you’re on hg38)

print("extraction start")

mhc_cpgs <- CpG_pval %>%
  # 0) Keep only chr‐6 entries via fast string‐filter
  filter(startsWith(CpG, "6:")) %>%
  
  # 1) Now split into chr, start, end
  separate(CpG, into = c("chr", "pos"), sep = ":",    remove = FALSE) %>%
  separate(pos, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  
  # 2) Filter on your MHC coordinates
  filter(start >= mhc_start, start <= mhc_end) 

#write the mhc_cpgs into a csv file
write.csv(mhc_cpgs, paste0(PATH_wk,"/scr/8_rerun/dat/8_zoom.csv"), row.names = FALSE)


############################

mhc_cpgs <- fread(paste0(PATH_wk,"/scr/8_rerun/dat/8_zoom.csv"), header=TRUE)
plot_name <- "Model2"
source(paste0(PATH_wk,"/scr/8_rerun/0_func_mhc_plot.R"))
}

#draw qqplot 
#draw qqplot
manhattan_dat <- CpG_pval
plot_name <- "Model2"
source(paste0(PATH_wk,"/scr/8_rerun/0_func_QQ_plot.R"))


# Convert p-values to chi-square statistics (df = 1)
chisq <- qchisq(1 - na.omit(manhattan_dat$pval), df = 1)

# Genomic inflation λ = median(observed chisq) / median(expected chisq)
lambda_gc <- median(chisq) / qchisq(0.5, df = 1)
lambda_gc



####################################################################################
####################################################################################
#Model 1 

#noew results of model 1 and model 3 (M1_manhattan,M3_manhattan)
load( paste0(PATH_wk,"/scr/8_rerun/results/4_eval.RData"), verbose=TRUE)
mhc_cpgs <- fread(paste0(PATH_wk,"/scr/8_rerun/dat/8_zoom.csv"), header=TRUE)

#input cogs_int and model_name
CpG_pval <- M1_manhattan #data set containg CpG and pval


M1_mhc_cpgs <- CpG_pval %>%
  # 0) Keep only chr‐6 entries via fast string‐filter
  filter(startsWith(CpG, "6:")) %>%
  
  # 1) Now split into chr, start, end
  separate(CpG, into = c("chr", "pos"), sep = ":",    remove = FALSE) %>%
  separate(pos, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  
  # 2) Filter on your MHC coordinates
  filter(start >= mhc_start, start <= mhc_end) 

#load M2 results: the mhc_cpgs 
missing_cpgs <- mhc_cpgs$CpG[!(mhc_cpgs$CpG %in% M1_mhc_cpgs$CpG)]

new_rows <- mhc_cpgs %>%
  filter(CpG %in% missing_cpgs) %>%
  transmute(
    CpG,
    chr      = as.character(chr),   # ← convert here
    start,
    end,
    pval  = 1,
    coef  = 0,
    pval_old = 1    # or drop this column if you don’t need it
  )
  # 3. Bind them onto your original M1_mhc_cpgs
M1_mhc_cpgs_filled <- bind_rows(M1_mhc_cpgs, new_rows)

#now plot the M1_mhc_cpgs_filled
mhc_cpgs <- M1_mhc_cpgs_filled
plot_name <- "Model1"
source(paste0(PATH_wk,"/scr/8_rerun/0_func_mhc_plot.R"))

M1_mhc_cpgs_filled[M1_mhc_cpgs_filled$CpG == "6:31745741-NA",]

#draw qqplot, note we need the complete manhattan plot pvalues 

M1_manhattan_old <- AA_simple_results[,c("Chunk_index","CpG","M1_pval", "M1_coef" , "M1_F_sdv")]
M1_manhattan_old$M1_pval[match(M1_manhattan$CpG,M1_manhattan_old$CpG)] <- M1_manhattan$pval
M1_manhattan_old <- M1_manhattan_old %>%  dplyr::rename(pval = `M1_pval`)
manhattan_dat <- M1_manhattan_old
plot_name <- "Model1"
source(paste0(PATH_wk,"/scr/8_rerun/0_func_QQ_plot.R"))




####################################################################################
####################################################################################
#Model 3

#noew results of model 1 and model 3 (M1_manhattan,M3_manhattan)
load( paste0(PATH_wk,"/scr/8_rerun/results/4_eval.RData"), verbose=TRUE)
mhc_cpgs <- fread(paste0(PATH_wk,"/scr/8_rerun/dat/8_zoom.csv"), header=TRUE)

#input cogs_int and model_name
CpG_pval <- M3_manhattan #data set containg CpG and pval


M3_mhc_cpgs <- CpG_pval %>%
  # 0) Keep only chr‐6 entries via fast string‐filter
  filter(startsWith(CpG, "6:")) %>%
  
  # 1) Now split into chr, start, end
  separate(CpG, into = c("chr", "pos"), sep = ":",    remove = FALSE) %>%
  separate(pos, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  
  # 2) Filter on your MHC coordinates
  filter(start >= mhc_start, start <= mhc_end) 

#load M2 results: the mhc_cpgs 
missing_cpgs <- mhc_cpgs$CpG[!(mhc_cpgs$CpG %in% M3_mhc_cpgs$CpG)]

new_rows <- mhc_cpgs %>%
  filter(CpG %in% missing_cpgs) %>%
  transmute(
    CpG,
    chr      = as.character(chr),   # ← convert here
    start,
    end,
    pval  = 1,
    coef  = 0,
    pval_old = 1    # or drop this column if you don’t need it
  )
  # 3. Bind them onto your original M1_mhc_cpgs
M3_mhc_cpgs_filled <- bind_rows(M3_mhc_cpgs, new_rows)

#now plot the M1_mhc_cpgs_filled
mhc_cpgs <- M3_mhc_cpgs_filled
plot_name <- "Model3"
source(paste0(PATH_wk,"/scr/8_rerun/0_func_mhc_plot.R"))

M3_mhc_cpgs_filled[M3_mhc_cpgs_filled$CpG == "6:31745741-NA",]

M3_mhc_cpgs_filled[M3_mhc_cpgs_filled$pval<1e-5,]
#no need for qqplot, since there are only 5535 cpgs are selected.
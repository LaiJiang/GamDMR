#after 1_test.R get the raw methylation data.
#we can see that there are some delicate signals in the localized region

#why some univaraite models are not significant
#why somnibus is not signficiant 

#step1: first verify if any of these cpgs in the region are significant with raw group comparision


i_gene <- 1

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

#save meth_long
load( file = paste0(PATH_wk, "results/10_sanity/meth_dat_", i_gene, ".RData"),verbose=TRUE)

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
ggsave(filename = paste0(PATH_wk, "results/10_sanity/2_colregion_", i_gene, ".jpg"), plot = ggplot2, width = 10, height = 6)


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
#this is actually the most left region. 
#now we zoom into the region and plot the metylation data

i_meth_long <- meth_long %>%
  filter(start >= i_df$region_start & start <= i_df$region_end)




# Quick plot
ggplot3 <- ggplot(i_meth_long, aes(x = start, y = meth, color = factor(AA_only))) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
  theme_minimal() + #add vertical line for i_position$Position
  labs(title = "",
       x = "Genomic Position",
       y = "Methylation Level",
       color = "AA Status") 


# Save the plot
ggsave(filename = paste0(PATH_wk, "results/10_sanity/2_col_meth_dat_", i_gene, ".png"), plot = ggplot3, width = 10, height = 6)


library(ggplot2)

ggplot3 <- ggplot(i_meth_long, aes(x = start, y = meth, color = factor(AA_only))) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
  theme_minimal() +
  labs(
    title = "",
    x = "Genomic Position",
    y = "Methylation Level",
    color = "AA Status"
  ) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15)
  )

# Save the plot
ggsave(filename = paste0(PATH_wk, "results/10_sanity/2_col_meth_dat_", i_gene, ".jpg"),
       plot = ggplot3, width = 10, height = 6)


############
#grant figure updates:



#check the report_materials for the conclusion and next steps for somnibus.

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
ggsave(paste0(PATH_wk, "results/10_sanity/manhattan_M1.jpg"), p1, width = 10, height = 5)
ggsave(paste0(PATH_wk, "results/10_sanity/manhattan_M2.jpg"), p2, width = 10, height = 5)
ggsave(paste0(PATH_wk, "results/10_sanity/manhattan_M3.png"), p3, width = 10, height = 5)
#M3 no signals at all all pvalues = 1
#M1 and M2 have some signals,but not significant min-pvalue = 0.001
print(min(AA_simple_results_chr$M1_pval))

print(min(AA_simple_results_chr$M2_pval))


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

#now plot the ftitted curves by AA status: 
source(paste0(PATH_scr10,"/0_func_plot_gam.R"))


#now plot the difference between AA=1 and AA=0 curves
source(paste0(PATH_scr10,"/0_gam_region.R"))

####################################################################################

#now test on the sub-region of interest 

# Get sub-region bounds
peak_region <- range(candidate_regions$start, na.rm = TRUE)
peak_start <- peak_region[1]
peak_end   <- peak_region[2]

cat("Peak region: ", peak_start, "-", peak_end, "\n")



# Subset to only CpGs within peak region
# Expand by ±500 bp or ±N positions as needed
window_pad <- 500
peak_start_window <- peak_start - window_pad
peak_end_window   <- peak_end + window_pad

meth_peak <- i_meth_long_cov %>%
  filter(start >= peak_start_window & start <= peak_end_window)

length(unique(meth_peak$start))  # Should now be >1

# Refit GAM only on peak region
gam_peak <- mgcv::gam(
  meth_arcsin ~ 
    s(start, bs = "cs") + 
    s(start, by = AA_only, bs = "cs") +
    AA_only +
    AgeCalc + Sex + Non.smoker +
    EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc +
    sv1 + sv2 + sv3 + sv4 + sv5 +
    s(FID, bs = "re"),
  data = meth_peak,
  method = "REML"
)


summary(gam_peak)

#In this ~1 kb region, methylation patterns vary significantly between individuals with vs. without allergic asthma (AA) — not due to global methylation shifts, 
#but due to differences in local CpG-specific patterns across the region.

#This result supports regional association, and justifies: Reporting this region as a differentially methylated region (DMR).

#Optionally validating with a region-based permutation test.






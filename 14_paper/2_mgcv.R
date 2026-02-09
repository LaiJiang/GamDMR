#this file collects the region  from 1_concor.R and evaluate the mgcv  results on this region.


i_gene <- "CNOT3"


PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

#save meth_long
load( file = paste0(PATH_wk, "results/14_paper/meth_dat_", i_gene, ".RData"),verbose=TRUE)

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
ggsave(filename = paste0(PATH_wk, "results/14_paper/2_mgcv_", i_gene, ".jpg"), plot = ggplot2, width = 10, height = 6)


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
#there are actually no somnibus results for this region.

if(FALSE){
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
ggsave(filename = paste0(PATH_wk, "results/14_paper/2_col_meth_dat_", i_gene, ".png"), plot = ggplot3, width = 10, height = 6)


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

}
############






#would that possible to poll that weak signals into a strong regional signal, with mgcv method ?

#!!!!
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
i_meth_long <- meth_long 

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

#!!!!
#now plot the ftitted curves by AA status: 
source(paste0(PATH_scr10,"/0_func_plot_gam.R"))

ggsave(filename = paste0(PATH_wk, "results/14_paper/2_col_gam_fitted_curves.jpg"), plot=ggplot_gam,width = 10, height = 5)  

#now plot the difference between AA=1 and AA=0 curves
source(paste0(PATH_scr10,"/0_gam_region.R"))


ggsave(filename = paste0(PATH_wk, "results/14_paper/2_col_gam_region_diff.jpg"), plot=ggplot_gam_region, width = 10, height = 5)

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






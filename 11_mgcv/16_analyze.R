#this file collect the second run of effect size calculation, and compare to the first run.
#and then filter by effect size and edf, and then detect DMR regions.

#this file collects the result of both batches. and record them.

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_results <- paste0(PATH_wk, "results/11_mgcv/")

#note this is different now!!!
PATH_results_B1 <- paste0(PATH_results, "/spacing/B1")
PATH_results_B2 <- paste0(PATH_results, "/spacing/B2/")
#PATH_results_B3 <- paste0(PATH_results, "/spacing/B3/")


library(data.table)
library(stringr)
################################################################################
#load function to collect results!!!
source(paste0(PATH_wk, "scr/11_mgcv/0_load_results_updated.R"))
results_B1 <- func_col_batch(PATH_results_B1)
results_B2 <- func_col_batch(PATH_results_B2)
#results_B3 <- func_col_batch(PATH_results_B3)


#merge the three results by job_id, but keep only intersecting job_ids
results_all_batch <- rbind(results_B1, results_B2)

print(dim(results_all_batch))

#remove the job_id column
results_all_batch <- results_all_batch[, -which(colnames(results_all_batch) == "job_id")]

#find  duplicated rows based on region_id
results_all_batch <- results_all_batch[!duplicated(results_all_batch), ]

print(dim(results_all_batch))

#save results_all_batch !!!
write.table(results_all_batch, 
            file = paste0(PATH_wk, "results/11_mgcv/16_results_all_batch.txt"), 
            row.names = FALSE, col.names = TRUE, sep = "\t")

##################################################
results_all_batch <- fread(paste0(PATH_wk, "results/11_mgcv/16_results_all_batch.txt"), 
                                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#now evalute the qqplot of the p-values in column s(start):AA_only
# Load libraries
library(ggplot2)

# Extract p-values
pvals <- results_all_batch$`s(start):AA_only`

edf <- results_all_batch$"edf_s(start):AA_only"
#assign the 0 value to the minimum nonzero value
min_nonzero <- min(pvals[pvals > 0], na.rm = TRUE)

pvals[pvals == 0] <- min_nonzero/2   # Assign a small value to zero p-values

#caveats: some pvalue ==0, some pvalue < 0. 
#we need to handle these cases, and decide if we want to remove them based on biological interpretation.!!!!!!

#remove the negative p-values
results_all_batch <- results_all_batch[results_all_batch$`s(start):AA_only` >= 0, ]
pvals <- pvals[pvals >= 0]

# Remove NAs or invalid values (e.g., negative, >1)
print(sum(is.na(pvals)))
print(sum((pvals<=0)))
print(sum((pvals>1)))


# Compute expected -log10(p-values) under uniform distribution
expected <- -log10(ppoints(length(pvals)))
observed <- -log10(sort(pvals))

# Make data frame for plotting
qq_data <- data.frame(expected = expected, observed = observed)

# Plot QQ plot
ggplot1 <- ggplot(qq_data, aes(x = expected, y = observed)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(
    title = "QQ Plot of smoothed AA_only effect p-values",
    x = "Expected -log10(p)",
    y = "Observed -log10(p)"
  ) +
  theme_minimal(base_size = 14)


  #save ggplot1
ggsave(filename = paste0(PATH_wk, "results/11_mgcv/16_qqplot_pvals_s_start_AA_only.jpg"), 
       plot = ggplot1, width = 8, height = 6, dpi = 300)

####
#check the edf of the s(start):AA_only, if they are strong enough
edf <- results_all_batch$"edf_s(start):AA_only"
plot(-log10(pvals),edf)

#check the R2
R2 <- results_all_batch$"R2"
plot(-log10(pvals),R2)


#check the maximum difference between the two groups
max_diff <- results_all_batch$"max_diff"
plot(-log10(pvals), max_diff)


#check the median difference between the two groups
mean_diff <- results_all_batch$"mean_diff"
plot(-log10(pvals), mean_diff)

#generate a data table containing the region_id, ,mean_diff, R2, edf, pvals
results_all_batch_filter <- results_all_batch %>%
  select(region_id, 
         edf= "edf_s(start):AA_only",
         pvals = `s(start):AA_only`, R2, N_cpgs, 
         mean_diff)

#add FRD adjusted p-values

# Compute FDR-adjusted p-values using Benjamini-Hochberg
fdr_pvals <- p.adjust(results_all_batch_filter$pvals, method = "BH")

results_all_batch_filter <- results_all_batch_filter %>%
    mutate(pval_FDR = fdr_pvals) %>% filter(pvals < 1e-5)

#write to a csv file
write.table(results_all_batch_filter, 
            file = paste0(PATH_wk, "results/11_mgcv/16_results_all_batch_filter.txt"), 
            row.names = FALSE, col.names = TRUE, sep = "\t")


####This script is a whatI run on a computer cluster to obtain region wise pvalues for testing the regional effect of AA on methylation.
#I have collected the results from the cluster and put the result in 16_results_all_batch_filter.txt.
#now write R script to filter DMRs based on the p-values, edf, R2 and mean_diff. 
#also give your justficiation for the filtering criteria.
################################################################################################
################################################################################################
################################################################################################
################################################################################################
#now filter by 
#FDR < 0.01 (for s(start):AA_only): ensures the regional shape difference by AA is robust after multiple testing across many regions.

#edf ≥ 1: guards against degenerate/over-penalized smooths; if edf is ~0, the AA-specific smooth contributed nothing. (You can raise to 1.5–2 if you want to emphasize clearly non-linear departures.)

#R² ≥ 0.30 (strict set): asks for a region whose methylation pattern is reasonably well captured by the model. This reduces “statistically significant but weakly explained” hits.

#mean_diff threshold:

#Strict (absolute): ≥ 0.01 on the model response scale you used (arcsin–sqrt of β). This flags regions with non-trivial average separation between AA vs control smooths across CpGs.

#Empirical (percentile): Because the arcsin scale isn’t always intuitive, the top 25% approach keeps the strongest quarter of effects among significant regions, adapting to your data’s natural spread.

dt <- results_all_batch_filter
# Expecting columns like:
# region_id, edf, pvals (for s(start):AA_only), R2, N_cpgs, mean_diff, FDR_pvals OR pval_FDR
# Standardize the FDR column name if both exist
if ("pval_FDR" %in% names(dt)) dt$FDR <- dt$pval_FDR
if (!"FDR" %in% names(dt) && "FDR_pvals" %in% names(dt)) dt$FDR <- dt$FDR_pvals

# ---- Basic sanity filters you likely want everywhere ----
# 1) Significant smooth interaction after multiple testing
# 2) Some minimum basis complexity (edf) so the smooth isn't effectively linear/no effect
# 3) Enough CpGs for a reliable regional fit
library(data.table)
library(dplyr)

base <- dt %>%
  filter(!is.na(FDR), FDR < 0.01,
         !is.na(edf), edf >= 1,     # edf ~ 1 = ~linear; keep >=1 to avoid degenerate fits
         !is.na(N_cpgs), N_cpgs >= 10)

# ---- Strict, fixed-threshold filter ----
# Rationale: demand decent fit + non-trivial effect size
strict <- base %>%
  filter(!is.na(R2), R2 >= 0.30,     # good regional fit
         !is.na(mean_diff), mean_diff >= 0.01)  # on the model’s response scale (arcsin-sqrt beta)

# ---- Empirical, data-driven filter (useful if mean_diff scaling is hard to interpret) ----
# Keep top quartile of mean_diff among significant hits that passed base criteria
q_md <- quantile(base$mean_diff, probs = 0.75, na.rm = TRUE)
empirical <- base %>%
  mutate(mean_diff_rank = percent_rank(mean_diff)) %>%
  filter(mean_diff >= q_md | mean_diff_rank >= 0.75)

# ---- Save outputs ----
outdir <-  paste0(PATH_wk, "results/11_mgcv/")

fwrite(base,     file.path(outdir, "16_dmrs_base.tsv"), sep = "\t")
fwrite(strict,   file.path(outdir, "16_dmrs_STRICT.tsv"), sep = "\t")
fwrite(empirical,file.path(outdir, "16_dmrs_EMPIRICAL.tsv"), sep = "\t")

# ---- Quick summary table to compare filters ----
summary_tbl <- tibble::tibble(
  total_input      = nrow(dt),
  base_kept        = nrow(base),
  strict_kept      = nrow(strict),
  empirical_kept   = nrow(empirical),
  mean_diff_q75    = as.numeric(q_md)
)
print(summary_tbl)
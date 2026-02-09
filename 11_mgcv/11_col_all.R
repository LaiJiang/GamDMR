#this file collects the result of all three batches. and record them.

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_results <- paste0(PATH_wk, "results/11_mgcv/")

PATH_results_B1 <- paste0(PATH_results, "/spacing/")
PATH_results_B2 <- paste0(PATH_results, "/spacing/B2/")
PATH_results_B3 <- paste0(PATH_results, "/spacing/B3/")


library(data.table)
library(stringr)
################################################################################
#load function to collect results
source(paste0(PATH_wk, "scr/11_mgcv/0_load_results.R"))
results_B1 <- func_col_batch(PATH_results_B1)
results_B2 <- func_col_batch(PATH_results_B2)
results_B3 <- func_col_batch(PATH_results_B3)


#merge the three results by job_id, but keep only intersecting job_ids
results_all_batch <- rbind(results_B1, results_B2, results_B3)


#remove the job_id column
results_all_batch <- results_all_batch[, -which(colnames(results_all_batch) == "job_id")]

#find  duplicated rows based on region_id
results_all_batch <- results_all_batch[!duplicated(results_all_batch), ]


#save results_all_batch
write.table(results_all_batch, 
            file = paste0(PATH_wk, "results/11_mgcv/results_all_batch.txt"), 
            row.names = FALSE, col.names = TRUE, sep = "\t")

##################################################
results_all_batch <- fread(paste0(PATH_wk, "results/11_mgcv/results_all_batch.txt"), 
                                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#now evalute the qqplot of the p-values in column s(start):AA_only
# Load libraries
library(ggplot2)

# Extract p-values
pvals <- results_all_batch$`s.start..AA_only`

edf <- results_all_batch$"edf_s.start..AA_only"
#assign the 0 value to the minimum nonzero value
min_nonzero <- min(pvals[pvals > 0], na.rm = TRUE)

pvals[pvals == 0] <- min_nonzero/2   # Assign a small value to zero p-values

#caveats: some pvalue ==0, some pvalue < 0. 
#we need to handle these cases, and decide if we want to remove them based on biological interpretation.!!!!!!

#remove the negative p-values
results_all_batch <- results_all_batch[results_all_batch$`s(start):AA_only` >= 0, ]

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
ggsave(filename = paste0(PATH_wk, "results/11_mgcv/11_qqplot_pvals_s_start_AA_only.jpg"), 
       plot = ggplot1, width = 8, height = 6, dpi = 300)

####
#check the edf of the s(start):AA_only, if they are strong enough
edf <- results_all_batch$"edf_s.start..AA_only"
plot(-log10(pvals),edf)

#check the R2
R2 <- results_all_batch$"R2"
plot(-log10(pvals),R2)
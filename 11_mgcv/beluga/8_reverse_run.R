
# This script runs an MGCV model on a single region of interest seperatelly

#i_region <- 57

#30% cpgs < 15
#10% cpgs <= 10
#1.4% cpgs <= 9

min_cpgs <- 10

N_jobs <- 300 # Number of jobs to split the regions into, this needs to be consistent with the SLURM_ARRAY_TASK_ID in the SLURM script

# Define paths
#PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
PATH_wk <- "~/scratch/UQAC/meth/"

PATH_results <- paste0(PATH_wk, "results/10_sanity/")
PATH_scr6 <- paste0(PATH_wk, "scr/6_whole/")
PATH_scr11 <- paste0(PATH_wk, "scr/11_mgcv/")
PATH_output <- paste0(PATH_wk, "results/11_mgcv/B2/")


# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No command-line arguments supplied. Please pass the SLURM_ARRAY_TASK_ID.")
}

# Convert the first argument to an integer
i_job_id <- as.integer(args[1])
cat("Job ID is:", i_job_id, "\n")



# Ensure output directory exists
if (!dir.exists(PATH_output)) dir.create(PATH_output, recursive = TRUE)

# Load libraries
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(tidyr)
library(mgcv)

# Load region info
region_file <- fread(paste0(PATH_scr11, "dat/region_file_1_chunk.csv"))


#split array 1:1737 into 999 arrays
vec_list <- split(1:nrow(region_file), cut(seq_along(1:nrow(region_file)), breaks = N_jobs, labels = FALSE))


#find the missing jobs ids
job_id_missing <- as.numeric(read.table(paste0(PATH_scr11, "dat/job_ids_missing_finished.txt"), header = FALSE, stringsAsFactors = FALSE)$V1)

original_job_id <- job_id_missing[i_job_id]

i_region_id_vector <- rev(vec_list[[original_job_id]])


for(i_region in i_region_id_vector) {

i_region_info <- region_file[i_region, ]
i_chunk_id <- i_region_info$data_chunk_id

# Load the corresponding data chunk
i_chunk <- fread(paste0(PATH_wk, "data/meth_split/chunk_", sprintf("%04d", i_chunk_id), ".csv"))

# Extract region of interest
i_region_data <- i_chunk %>%
  filter(chr == i_region_info$chr,
         start >= i_region_info$region_start,
         start <= i_region_info$region_end)

if (nrow(i_region_data) >= min_cpgs) {

  setDT(i_region_data)

  # Convert _meth columns to numeric
  meth_cols <- grep("_meth$", names(i_region_data), value = TRUE)
  i_region_data[, (meth_cols) := lapply(.SD, as.numeric), .SDcols = meth_cols]

  # Melt to long format
  meth_long <- melt(i_region_data,
                    id.vars = "start",
                    measure.vars = patterns("_meth$"),
                    variable.name = "sample",
                    value.name = "meth",
                    na.rm = TRUE)
  meth_long[, sample := str_replace(as.character(sample), "_meth", "")]

  # Load phenotype data
  load(file = paste0(PATH_scr6, "4_1_data.RData"), verbose = TRUE)
  pheno_file <- as_tibble(pheno_file)

  # Join phenotype data
  meth_long <- meth_long %>%
    left_join(pheno_file %>% select(ID, AA_only), by = c("sample" = "ID")) %>%
    na.omit()

  i_pheno_file <- pheno_file %>% rename(sample = ID) %>% select(-AA_only)

  i_meth_long_cov <- meth_long %>%
    left_join(i_pheno_file, by = "sample") %>%
    na.omit()

  if (nrow(i_meth_long_cov) < min_cpgs) {
  cat("Skipping region", i_region, "after join due to insufficient CpGs.\n")
  next
  }

  # Arcsin transform
  i_meth_long_cov$meth_arcsin <- asin(sqrt(i_meth_long_cov$meth))

  # Fix type bugs
  for (v in colnames(i_meth_long_cov)) {
    if (is.matrix(i_meth_long_cov[[v]]) || is.array(i_meth_long_cov[[v]])) {
      i_meth_long_cov[[v]] <- as.numeric(i_meth_long_cov[[v]])
    }
  }

  # Define basis dimension
  k_cpg <- max(5, min(10, floor(length(unique(i_meth_long_cov$start)) / 20)))

  # Fit GAM model
  gam_model <- mgcv::gam(meth_arcsin ~ 
    s(start, bs = "cs", k = k_cpg) +
    s(start, by = AA_only, bs = "cs", k = k_cpg) +
    AA_only + AgeCalc + Sex + Non.smoker +
    EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc +
    sv1 + sv2 + sv3 + sv4 + sv5 +
    s(FID, bs = "re"),
    data = i_meth_long_cov,
    method = "REML")

  # Extract results
  feature_list <- c("(Intercept)", "AA_only", "AgeCalc", "Sex", "Non.smoker",
                    "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc", 
                    "sv1", "sv2", "sv3", "sv4", "sv5")
  feature_smooth_terms <- c("s(start)", "s(start):AA_only", "s(FID)")

  param_coef <- summary(gam_model)$p.table
  smooth_terms <- summary(gam_model)$s.table

  present_features <- intersect(feature_list, rownames(param_coef))
  result_pvalue_fixed_effect <- rep(1, length(feature_list))
  names(result_pvalue_fixed_effect) <- feature_list
  result_pvalue_fixed_effect[present_features] <- param_coef[present_features, "Pr(>|t|)"]

  present_smooth <- intersect(feature_smooth_terms, rownames(smooth_terms))
  result_smooth_terms <- rep(1, length(feature_smooth_terms))
  names(result_smooth_terms) <- feature_smooth_terms
  result_smooth_terms[present_smooth] <- smooth_terms[present_smooth, "p-value"]

  edf_smooth_terms <- rep(0, length(feature_smooth_terms))
  names(edf_smooth_terms) <- feature_smooth_terms
  edf_smooth_terms[present_smooth] <- smooth_terms[present_smooth, "edf"]

  # Model fit statistics
  result_stats <- c(result_pvalue_fixed_effect,
                    result_smooth_terms,
                    edf_smooth_terms,
                    "R2" = summary(gam_model)$r.sq,
                    "AIC" = AIC(gam_model),
                    "Deviance_explained" = summary(gam_model)$dev.expl,
                    "REML" = gam_model$gcv.ubre,
                    "N_cpgs" = length(unique(i_meth_long_cov$start)),
                    "N_samples" = length(unique(i_meth_long_cov$sample)),
                    "data_chunk_id" = i_chunk_id,
                    "region_id" = i_region)


   #important !!!!
   #Look for edf > 0.5 or 1 before trusting the s(start):AA_only p-value.

  # Convert to data.frame (single-row, unnamed)
  result_stats_df <- as.data.frame(t(as.numeric(result_stats)))

  # Save output
  # Append to job-specific .txt file
  fwrite(
  result_stats_df,
  file = paste0(PATH_output, "8_results_job_", original_job_id, ".txt"),
  append = TRUE,
  sep = "\t",
  col.names = FALSE
  )

}

else {
  cat("Skipping region", i_region, "due to insufficient CpGs.\n")
}
}
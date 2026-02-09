#collect the existing results, and see if anything is missing

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
PATH_dat <- file.path(PATH_wk, "results/11_mgcv/spacing/")
PATH_scr11 <- paste0(PATH_wk, "scr/11_mgcv/")
library(data.table)
library(stringr)

all_files <- list.files(PATH_dat, pattern = "\\.txt$", full.names = TRUE)

results_col <- NULL 
for(loop_index in 1:length(all_files)){

i_file <- all_files[loop_index]

i_file_dat <- read.table(i_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#this job ID in the file name
i_job_id <- as.numeric(str_extract(i_file, "(?<=7_results_job_)\\d+(?=\\.txt)"))


feature_fixed <- c("(Intercept)", "AA_only", "AgeCalc", "Sex", "Non.smoker",
                    "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc", 
                    "sv1", "sv2", "sv3", "sv4", "sv5")
feature_smooth_terms <- c("s(start)", "s(start):AA_only", "s(FID)")
feature_edf <- paste0("edf_",c("s(start)", "s(start):AA_only", "s(FID)"))
# Model fit statistics
features_column_names <- c(feature_fixed,
                    feature_smooth_terms,
                    feature_edf,
                    "R2" ,
                    "AIC"  ,
                    "Deviance_explained" ,
                    "REML" ,
                    "N_cpgs" ,
                    "N_samples",
                    "data_chunk_id",
                    "region_id")

colnames(i_file_dat) <- features_column_names


#check the last line
last_line <- tail(i_file_dat, 1)


#this shows the edf  versus the p-value, which means the p-value is reliable 
plot(i_file_dat$"edf_s(start):AA_only", i_file_dat$"s(start):AA_only")


#plot qqplot of the p-values
qqnorm(i_file_dat$"s(start):AA_only")
qqline(i_file_dat$"s(start):AA_only")

#the last line of the region_id procssed by this job
i_last_region_id <- last_line$region_id


results_col <- rbind(results_col, 
                      c(i_job_id, i_last_region_id))


}


###now we compare with the true last lines 

# Load region info
region_file <- fread(paste0(PATH_scr11, "dat/region_file_1_chunk.csv"))

N_jobs = 300
#split array 1:1737 into 300 arrays
vec_list <- split(1:nrow(region_file), cut(seq_along(1:nrow(region_file)), breaks = N_jobs, labels = FALSE))

result_true <- NULL 
for(i_job_id in 1:N_jobs){

i_region_id_vector <- vec_list[[i_job_id]]

#the last vaue of i_region_id_vector
i_last_region_id <- tail(i_region_id_vector, 1)

result_true <- rbind(result_true, 
                      c(i_job_id, i_last_region_id))
}


colnames(result_true) <- c("job_id", "last_region_id")
colnames(results_col) <- c("job_id", "col_last_region_id")


#merge the two results
results_merged <- merge(result_true, results_col, by = "job_id")

#extract the rows where the last_region_id is not equal to col_last_region_id
results_mismatch <- results_merged[results_merged$last_region_id != results_merged$col_last_region_id, ]        


#calculate the differences between the last_region_id and col_last_region_id
results_mismatch$diff <- as.numeric(results_mismatch$last_region_id) - as.numeric(results_mismatch$col_last_region_id)  

results_mismatch[which.max(results_mismatch$diff), ]


#####################################
#alternatively, collect the log files to see if the *.out file contains the word Finished
#also, check if the *.err file contains the word error

log_dir <- "C:/Per/LaiJiang/Project/UQAC/meth/results/11_mgcv/logs"


# List files
out_files <- list.files(log_dir, pattern = "\\.out$", full.names = TRUE)
err_files <- list.files(log_dir, pattern = "\\.err$", full.names = TRUE)

# Initialize result holders
out_missing_finished <- c()
err_with_errors <- c()
err_time_limit <- c()
# Check .out files for "Finished"
for (f in out_files) {
  lines <- readLines(f, warn = FALSE)
  if (!any(grepl("Finished", lines))) {
    out_missing_finished <- c(out_missing_finished, basename(f))
  }
}

# Check .err files for "error" / "Error" / "ERROR"
for (f in err_files) {
  lines <- readLines(f, warn = FALSE)
  if (any(grepl("(?i)\\berror\\b", lines))) {
    err_with_errors <- c(err_with_errors, basename(f))
  }
}

# Output summary
cat("*.out files missing 'Finished':\n")
print(out_missing_finished)

cat("\n*.err files containing 'error' variants:\n")
print(err_with_errors)




# Check .err files for "DUE TO TIME LIMIT" 
for (f in err_files) {
  lines <- readLines(f, warn = FALSE)
    if (any(grepl("(?i)DUE TO TIME LIMIT", lines))) {
    err_time_limit <- c(err_time_limit, basename(f))
    }
}

length(intersect(err_time_limit, err_with_errors))
length(err_time_limit)
length(err_with_errors)

#extract the job IDs from the out_missing_finished names
job_ids_missing_finished <- as.numeric(sub(".*_(\\d+)\\.out", "\\1", out_missing_finished))

job_ids_error <- as.numeric(sub(".*_(\\d+)\\.err", "\\1", err_with_errors))


#the job IDs that are in err_with_errors but not out_missing_finished 
job_ids_error_finished <- setdiff(job_ids_error,job_ids_missing_finished)

results_merged[results_merged$job_id%in%job_ids_error_finished,]

#save job_ids_missing_finished
write.table(job_ids_missing_finished, 
            file = paste0(PATH_scr11, "dat/job_ids_missing_finished.txt"), 
            row.names = FALSE, col.names = FALSE)


#save results_merged
write.table(results_merged, 
            file = paste0(PATH_scr11, "dat/results_merged.txt"), 
            row.names = FALSE, col.names = TRUE, sep = "\t")


#
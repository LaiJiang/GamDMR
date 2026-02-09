#next: collect batch 2 results (from BMI runs), and compare to the results_merged to see if still any region is missing.


#collect the existing results, and see if anything is missing

#  scp laj773@rorqual.calculquebec.ca:/home/laj773/scratch/UQAC/meth/results/11_mgcv/BMI_B2/*.txt /mnt/c/Per/LaiJiang/Project/UQAC/meth/results/11_mgcv/BMI_B2/

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
PATH_dat <- file.path(PATH_wk, "results/11_mgcv/BMI_B2/")
PATH_scr11 <- paste0(PATH_wk, "scr/11_mgcv/")
library(data.table)
library(stringr)

all_files <- list.files(PATH_dat, pattern = "\\.txt$", full.names = TRUE)

results_col <- NULL 
for(loop_index in 1:length(all_files)){

i_file <- all_files[loop_index]

i_file_dat <- read.table(i_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#this job ID in the file name
i_job_id <- as.numeric(str_extract(i_file, "(?<=8_results_job_)\\d+(?=\\.txt)"))


feature_fixed <- c("(Intercept)", "AA_only", "AgeCalc", "Sex", "Non.smoker",
                    "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc", 
                    "sv1", "sv2", "sv3", "sv4", "sv5", "BMI")
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
first_line <- head(i_file_dat, 1)

#this shows the edf  versus the p-value, which means the p-value is reliable 
plot(i_file_dat$"edf_s(start):AA_only", i_file_dat$"s(start):AA_only")


#plot qqplot of the p-values
qqnorm(i_file_dat$"s(start):AA_only")
qqline(i_file_dat$"s(start):AA_only")

#the last line of the region_id procssed by this job
i_last_region_id <- last_line$region_id

i_first_region_id <- first_line$region_id

results_col <- rbind(results_col, 
                      c(i_job_id, i_first_region_id, i_last_region_id))


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
i_first_region_id <- i_region_id_vector[1]

result_true <- rbind(result_true, 
                      c(i_job_id, i_first_region_id,i_last_region_id))
}


colnames(result_true) <- c("job_id", "first_region_id","last_region_id")
colnames(results_col) <- c("job_id","col_first_region_id", "col_last_region_id")


#merge the two results
results_merged <- merge(result_true, results_col, by = "job_id")


#save results_merged
write.table(results_merged, 
            file = paste0(PATH_scr11, "dat/23_results_merged_B2.txt"), 
            row.names = FALSE, col.names = TRUE, sep = "\t")

###################################################################
#collect results of 8_col.R and 9_col.R, see if any region is missing
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
PATH_dat <- file.path(PATH_wk, "results/11_mgcv/BMI_B2/")
PATH_scr11 <- paste0(PATH_wk, "scr/11_mgcv/")
library(data.table)
library(stringr)
library(dplyr)

results_8_col <- fread(paste0(PATH_scr11, "dat/21_results_merged.txt"))
results_9_col <- fread(paste0(PATH_scr11, "dat/23_results_merged_B2.txt"))


#merge the two results by job_id, but keep only intersecting job_ids
results_merged <- merge(results_8_col, results_9_col, by = "job_id", suffixes = c("_B1", "_B2"))

#the Batch 3 regions need to be run again
B3_regions <- results_merged[which(results_merged$col_last_region_id_B1 < results_merged$col_last_region_id_B2),] 



B3_regions <- B3_regions %>% #keep only job_id, col_last_region_id_B1,col_last_region_id_B2
    select(job_id, col_last_region_id_B1, col_last_region_id_B2) %>% 
    rename(job_id = job_id,
           first_region = col_last_region_id_B1,
           end_region = col_last_region_id_B2)

print(B3_regions)

#no need to run batch 3, as the regions are already processed in batch 2!!!!

if(FALSE){

    
#save job_ids_missing_finished
write.table(B3_regions, 
            file = paste0(PATH_scr11, "dat/B3_job_region_ids_missing_finished.txt"), 
            row.names = FALSE, col.names = FALSE)

}
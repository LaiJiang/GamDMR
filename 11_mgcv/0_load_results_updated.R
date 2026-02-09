#this code chunk loads the results from the previous runs


func_col_batch<- function(PATH_dat){

all_files <- list.files(PATH_dat, pattern = "\\.txt$", full.names = TRUE)

results_col <- NULL 
for(loop_index in 1:length(all_files)){

i_file <- all_files[loop_index]

i_file_dat <- read.table(i_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#this job ID in the file name
#i_job_id <- as.numeric(str_extract(i_file, "(?<=7_results_job_)\\d+(?=\\.txt)"))
i_job_ids <- as.numeric(str_extract(i_file, "(?<=_results_job_)\\d+(?=\\.txt)"))


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
                    "region_id",
                    "max_diff" ,
                    "mean_diff" )

colnames(i_file_dat) <- features_column_names


i_file_dat$job_id <- i_job_ids

results_col <- rbind(results_col, 
                     i_file_dat)


}

results_col

}

#collect the somnibus results from Beluga and evaluate them.
#this is the end-result of running all scripts on Beluga


PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

#repo to load results
PATH_somnibus_stage1 <- paste0(PATH_wk,"/results/somnibus_all/somnibus/")
PATH_somnibus_stage2 <- paste0(PATH_wk,"/results/somnibus_all/somnibus_stage2/")
PATH_somnibus_stage3 <- paste0(PATH_wk,"/results/somnibus_all/somnibus_stage3/results/")


#repo to outut results
file_output_stage1 <- paste0(PATH_wk,"/results/9_regional/13_col_stage1.csv")
file_output_stage2 <- paste0(PATH_wk,"/results/9_regional/13_col_stage2.csv")
file_output_stage3 <- paste0(PATH_wk,"/results/9_regional/13_col_stage3.csv")


##################################################################
####define a function to collect the results from the somnibus runns
#input: PATH_somnibus_results
#input: file_output: where to output the results
func_collect_somnibus <- function(PATH_somnibus_beluga, file_output){

for(i_file in list.files(PATH_somnibus_beluga)){

load(paste0(PATH_somnibus_beluga, i_file),verbose = TRUE)

n_regions <- length(outs)

j_region = 1

if(n_regions > 0){

for(j_region in 1:n_regions){

    region_name <- names(outs)[j_region]
    #print(paste0("Processing region: ", region_name))
    
    # Extract the result for the current region
    res <- outs[[region_name]]
    
    # Check if the result is valid
    if (is.null(res)) {
        next
    }
    
    # Print summary of the result
    # print(summary(res))
    
    # Extract positions and covariate effects
    pos <- res$uni.pos

    if(FALSE){
    pos <- res$uni.pos
    beta_AA <- res$Beta.out[,2]
    se_AA <- res$SE.out[,2]
    
    # Plot the smoothed effect of AA_only
    plot(pos, beta_AA, type = "l", lwd = 2, col = "blue",
         ylab = "Smoothed effect of AA_only", xlab = "Genomic position",
         main = paste("Smoothed effect of AA_only on methylation in", region_name))
    
    # Add confidence interval
    lines(pos, beta_AA + 1.96 * se_AA, col = "blue", lty = 2)
    lines(pos, beta_AA - 1.96 * se_AA, col = "blue", lty = 2)
    abline(h = 0, col = "gray50", lty = 3)
    }


    pvals <- res$reg.out
    
    pval_AA_only <- pvals[2,3]

    #output: pval_AA_only, chunk_id, region_id
    chunk_id <- as.numeric(sub(".*_(\\d+)_.*", "\\1", i_file))

    region_id <- j_region

    region_start <- min(pos)
    region_end <- max(pos)
   
    ij_output <- c(
        chunk_id ,
        region_id ,
        region_name ,
        region_start ,
        region_end ,
        pval_AA_only   )
    #append output to the file_output txt file
    write.table(t(ij_output), file = file_output, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
}
}

}

}

##################################################################

func_collect_somnibus(PATH_somnibus_stage1, file_output_stage1)
func_collect_somnibus(PATH_somnibus_stage2, file_output_stage2)


#stage 3 is different, it has no region names, so we need to handle it differently

PATH_somnibus_beluga = PATH_somnibus_stage3
file_output = file_output_stage3

i_file = list.files(PATH_somnibus_beluga)[1]

for(i_file in list.files(PATH_somnibus_beluga)){

load(paste0(PATH_somnibus_beluga, i_file),verbose = TRUE)

n_regions <- length(outs)

j_region = 1

if(n_regions > 0){

    #no region names in stage 3, so we will use a placeholder
    region_name <- NA
    #print(paste0("Processing region: ", region_name))
    
    # Extract the result for the current region
    res <- outs
    
    # Check if the result is valid
    if (is.null(res)) {
        next
    }
    
    # Print summary of the result
    # print(summary(res))
    
    # Extract positions and covariate effects
    pos <- res$uni.pos

    if(FALSE){
    pos <- res$uni.pos
    beta_AA <- res$Beta.out[,2]
    se_AA <- res$SE.out[,2]
    
    # Plot the smoothed effect of AA_only
    plot(pos, beta_AA, type = "l", lwd = 2, col = "blue",
         ylab = "Smoothed effect of AA_only", xlab = "Genomic position",
         main = paste("Smoothed effect of AA_only on methylation in", region_name))
    
    # Add confidence interval
    lines(pos, beta_AA + 1.96 * se_AA, col = "blue", lty = 2)
    lines(pos, beta_AA - 1.96 * se_AA, col = "blue", lty = 2)
    abline(h = 0, col = "gray50", lty = 3)
    }


    pvals <- res$reg.out
    
    pval_AA_only <- pvals[2,3]

    #output: pval_AA_only, chunk_id, region_id
    chunk_id <- data_chunk_id

    region_id <- as.numeric(sub(".*_(\\d+).*", "\\1", i_file))

    region_start <- min(pos)
    region_end <- max(pos)
   
    ij_output <- c(
        chunk_id ,
        region_id ,
        region_name ,
        region_start ,
        region_end ,
        pval_AA_only   )
    #append output to the file_output txt file
    write.table(t(ij_output), file = file_output, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)

}

}

##################################################################
##################################################################
##################################################################

#now collect all csv files
stage1_df <- read.csv(file_output_stage1, header = FALSE)
stage2_df <- read.csv(file_output_stage2, header = FALSE)
stage3_df <- read.csv(file_output_stage3, header = FALSE)

#combine all dataframes
all_df <- rbind(stage1_df, stage2_df, stage3_df)
#rename columns
colnames(all_df) <- c("chunk_id", "region_id", "region_name", "region_start", "region_end", "pval_AA_only")

#remove duplicated rows
all_df <- all_df[!duplicated(all_df), ]

#assign a region ID as chunk_id + region_id
all_df$region_id <- paste0(all_df$chunk_id, "_", all_df$region_id) 

#remove the column chunk_id region_name
all_df$region_name <- NULL
all_df$chunk_id <- NULL

#save the file to csv 
write.csv(all_df, file = paste0(PATH_wk,"/results/9_regional/13_col_combined.csv"), row.names = FALSE)
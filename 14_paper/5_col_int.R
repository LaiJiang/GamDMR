#collect samples of inter-cpg correlation distributions

if(FALSE){

    library(data.table)
    library(dplyr)
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

PATH_chunk <- paste0(PATH_wk, "dat/")
PATH_scr11 <- paste0(PATH_wk, "/scr/11_mgcv/")

chunk1  <-  fread(paste0(PATH_chunk,"chunk_0001.csv"))

chunk1_min <- min(chunk1$start)
chunk1_max <- max(chunk1$start)


load(file = paste0(PATH_wk,"/scr/14_paper/results/4_DMR_M1_M2_overlap.RData"),verbose=TRUE)

print(dim(dmr_gene_file))

dmr_gene_file$region_size <- dmr_gene_file$region_end - dmr_gene_file$region_start +1

#region file information
region_file <- data.table::fread(file.path(PATH_scr11, "dat/region_file_1_chunk.csv"))

library(data.table)

# ensure both are data.table
setDT(dmr_gene_file)
setDT(region_file)

# perform merge by chr, region_start, region_end
dmr_merged <- merge(
    dmr_gene_file,
    region_file[, .(chr, region_start, region_end, n_cpgs)],
    by = c("chr", "region_start", "region_end"),
    all.x = TRUE
)
}



n_cpg_list <- dmr_merged$n_cpgs



data_chunk_list <- c("0014", "0001", "0002", "0003", "0001", "017", "0226", "0227", "0907")
data_chunk_list <- paste0("chunk_", data_chunk_list, ".csv")


n_cpg_mat <- NULL 

for(j_chunk in data_chunk_list){
    print(paste0("Loading data chunk: ", j_chunk))
    chunk1  <-  fread(paste0(PATH_chunk, j_chunk))

    cpg_cor_list <- NULL

for(i in n_cpg_list){



   print( paste0("Processing n_cpgs = ", i) )

i_N_cpg <- n_cpg_list[i]



random.seed <- 12345
set.seed(random.seed)



#randomly select a row 
rand_start_row <- sample(1:(nrow(chunk1) - i_N_cpg + 1), 1)
rand_end_row <- rand_start_row + i_N_cpg -1

#now the random region is:
rand_region <- chunk1[rand_start_row:rand_end_row, ]


#1. match the IDs of AA and control 
# Extract only the methylation proportion columns (those that end with '_meth')
meth_data <- rand_region[, grep("_meth$", colnames(rand_region)), with = FALSE]

# Convert the data frame to a matrix
meth_matrix <- as.matrix(meth_data)

#colect the IDs for meth_matrix
meth_ID <- colnames(meth_matrix)
meth_ID <- sub("^X", "", meth_ID)      # Remove the leading "X"
meth_ID <- sub("_meth$", "", meth_ID)   # Remove the trailing "_meth"
meth_ID <- gsub("\\.", "-", meth_ID)     # Replace the dot with a dash


# columns with all NA
all_na_cols <- apply(meth_matrix, 2, function(x) all(is.na(x)))

# columns with all zero (excluding NA)
all_zero_cols <- apply(meth_matrix, 2, function(x) all(x == 0, na.rm = TRUE))

cols_to_keep <- !(all_na_cols | all_zero_cols)

meth_filtered <- meth_matrix[, cols_to_keep, drop = FALSE]

cor_mat <- cor(t(meth_filtered), use = "pairwise.complete.obs")

upper_vals <- cor_mat[upper.tri(cor_mat)]
median_corr <- median(upper_vals, na.rm = TRUE)
cpg_cor_list <- c(cpg_cor_list, median_corr)

}

n_cpg_mat <- rbind(n_cpg_mat, cpg_cor_list)
}

#for each column, randomly select 1 value to represent that n_cpgs
set.seed(12345)
final_cor_vec <- apply(n_cpg_mat, 2, function(x) max(abs(x)))

dmr_merged$inter_cpg_cor <- final_cor_vec
save(dmr_merged, file = paste0(PATH_wk,"/scr/14_paper/results/5_col_int.RData"))


dmr_merged$cpg_density <- dmr_merged$n_cpgs / dmr_merged$region_size

load( file = paste0(PATH_wk,"/scr/14_paper/results/4_DMR_M1_M2_overlap.RData"),verbose=TRUE)
table(DMR_M2$overlap_CpG)

dmr_merged$overlap_CpG <- DMR_M2$overlap_CpG[match(dmr_merged$region_id, DMR_M2$region_id)]

DMRs_features <- merge(dmr_merged, DMR_M2, by = c("chr", "region_start", "region_end"), all.x = TRUE)


DMRs_features[DMRs_features$overlap_CpG.y==1,]

#remove overlap_CpG.x
DMRs_features <- DMRs_features %>%
  dplyr::select(-overlap_CpG.x) %>%
  dplyr::rename(overlap_CpG = overlap_CpG.y)


save(DMRs_features, file = paste0(PATH_wk,"/scr/14_paper/results/5_DMRs_features.RData"))


#now report the regional concordance:all.x
library(data.table)

vars <- c("region_size", "n_cpgs", "cpg_density", "inter_cpg_cor")

result_table <- lapply(vars, function(v) {
    x1 <- DMRs_features[overlap_CpG == 1][[v]]
    x0 <- DMRs_features[overlap_CpG == 0][[v]]

    data.table(
        feature = v,
        median_overlap = median(x1, na.rm=TRUE),
        median_nonoverlap = median(x0, na.rm=TRUE),
        pvalue = wilcox.test(x1, x0)$p.value
    )
})

result_table <- rbindlist(result_table)
result_table


fit <- glm(overlap_CpG ~ region_size + n_cpgs + cpg_density + inter_cpg_cor,
           data = DMRs_features,
           family = binomial)
summary(fit)

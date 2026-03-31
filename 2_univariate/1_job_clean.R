args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])

base_dir <- getwd()
data_dir <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")

meth_file_loc <- file.path(data_dir, sprintf("chunk_%04d.csv", task_id))
rdata_loc <- file.path(data_dir, "4_1_data.RData")
func_file <- file.path(base_dir, "6_0_func_single_clean.R")

dir.create(results_dir, showWarnings = FALSE)

meth_file <- read.csv(meth_file_loc, sep = "\t")
load(rdata_loc)

meth_data <- meth_file[, grep("_meth$", colnames(meth_file))]
meth_matrix <- as.matrix(meth_data)

meth_id <- colnames(meth_matrix)
meth_id <- sub("^X", "", meth_id)
meth_id <- sub("_meth$", "", meth_id)
meth_id <- gsub("\.", "-", meth_id)

meth_matrix <- meth_matrix[, match(pheno_file$ID, meth_id)]

meth_file$location <- paste(meth_file$chr, paste(meth_file$start, meth_file$end, sep = "-"), sep = ":")

source(func_file)

for (i_cpg in seq_len(nrow(meth_matrix))) {
  meth_response <- asin(sqrt(meth_matrix[i_cpg, ]))
  meth_cpg_info <- meth_file$location[i_cpg]

  pheno_file$meth_response <- meth_response
  pheno_file_feed <- na.omit(pheno_file)

  if (nrow(pheno_file_feed) > 14 && sd(pheno_file_feed$meth_response) != 0) {
    source(func_file)

    write.table(t(output_simple_results),
      file = file.path(results_dir, paste0("result_", task_id, ".txt")),
      append = TRUE, col.names = FALSE, row.names = FALSE, sep = "\t"
    )
  }
}


suppressPackageStartupMessages({
  library(data.table)
})

PATH_wk <- path.expand(Sys.getenv("METH_BASE_DIR", "~/scratch/UQAC/meth/"))
PATH_scr11 <- file.path(PATH_wk, "scr", "11_mgcv")
PATH_save <- file.path(PATH_wk, "results", "15_revision", "5_cv")
dir.create(PATH_save, recursive = TRUE, showWarnings = FALSE)

env <- new.env(parent = emptyenv())
load(file.path(PATH_scr11, "dat", "18_pheno_BMI.RData"), envir = env)
if (!exists("pheno_file", envir = env, inherits = FALSE)) stop("pheno_file is missing.")
pheno_file <- as.data.table(get("pheno_file", envir = env))

all_FIDs <- unique(as.character(pheno_file$FID))
all_FIDs <- all_FIDs[!is.na(all_FIDs) & nzchar(all_FIDs)]
n_FIDs <- length(all_FIDs)
n_train <- 135L
n_split <- 100L

stopifnot(n_FIDs == 169L)
set.seed(20260727)

train_FID_list <- lapply(seq_len(n_split), function(i) {
  sort(sample(all_FIDs, size = n_train, replace = FALSE))
})
train_FID_df <- as.data.table(do.call(rbind, train_FID_list))
setnames(train_FID_df, paste0("train_FID_", seq_len(n_train)))
train_FID_df[, splitID := seq_len(n_split)]
setcolorder(train_FID_df, c("splitID", paste0("train_FID_", seq_len(n_train))))

result_file <- file.path(PATH_save, "100_family_training_splits.csv")
cluster_file <- file.path(PATH_scr11, "dat", "100_family_training_splits.csv")
fwrite(train_FID_df, result_file)
fwrite(train_FID_df, cluster_file)

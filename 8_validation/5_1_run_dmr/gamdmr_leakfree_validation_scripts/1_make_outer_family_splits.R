###############################################################################
# 1_make_outer_family_splits.R
#
# Create the 100 outer family-aware train/test splits BEFORE any
# outcome-dependent methylation feature discovery.
###############################################################################

library(data.table)

PATH_wk <- path.expand("~/scratch/UQAC/meth/")

PATH_scr11 <- file.path(PATH_wk, "scr/11_mgcv")
PATH_out   <- file.path(
  PATH_wk,
  "results/15_revision/5_cv_leakfree"
)

dir.create(PATH_out, recursive = TRUE, showWarnings = FALSE)

pheno_file_path <- file.path(
  PATH_scr11,
  "dat/18_pheno_BMI.RData"
)

###############################################################################
# Load phenotype data
###############################################################################

e <- new.env()

load(
  pheno_file_path,
  envir = e
)

if (!exists("pheno_file", envir = e)) {
  stop("18_pheno_BMI.RData does not contain pheno_file.")
}

pheno <- as.data.table(e$pheno_file)

rm(e)

required_cols <- c(
  "ID",
  "FID",
  "AA_only"
)

missing_cols <- setdiff(
  required_cols,
  names(pheno)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing phenotype columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

pheno[, ID := as.character(ID)]
pheno[, FID := as.character(FID)]

pheno <- pheno[
  !is.na(ID) &
  !is.na(FID) &
  nzchar(ID) &
  nzchar(FID)
]

###############################################################################
# Keep the original FID order so the same seed can reproduce your old splits
###############################################################################

all_FIDs <- unique(pheno$FID)

n_FIDs <- length(all_FIDs)

cat("Number of families:", n_FIDs, "\n")

# Your manuscript cohort currently has 169 families
if (n_FIDs == 169L) {
  n_train <- 135L
} else {
  n_train <- floor(0.80 * n_FIDs)
}

n_split <- 100L

cat("Training families:", n_train, "\n")
cat("Test families:", n_FIDs - n_train, "\n")

###############################################################################
# Generate splits
###############################################################################

set.seed(20260727)

split_list <- lapply(
  seq_len(n_split),
  function(split_id) {

    train_FIDs <- sample(
      all_FIDs,
      size = n_train,
      replace = FALSE
    )

    test_FIDs <- setdiff(
      all_FIDs,
      train_FIDs
    )

    rbind(
      data.table(
        splitID = split_id,
        FID = train_FIDs,
        set = "train"
      ),
      data.table(
        splitID = split_id,
        FID = test_FIDs,
        set = "test"
      )
    )
  }
)

split_long <- rbindlist(split_list)

###############################################################################
# Sanity checks
###############################################################################

for (s in seq_len(n_split)) {

  x <- split_long[
    splitID == s
  ]

  train_ids <- x[
    set == "train",
    FID
  ]

  test_ids <- x[
    set == "test",
    FID
  ]

  stopifnot(
    length(intersect(train_ids, test_ids)) == 0L
  )

  stopifnot(
    length(unique(c(train_ids, test_ids))) == n_FIDs
  )
}

###############################################################################
# Show AA distribution in every split
###############################################################################

subject_split <- merge(
  split_long,
  pheno[, .(ID, FID, AA_only)],
  by = "FID",
  allow.cartesian = TRUE
)

split_summary <- subject_split[
  ,
  .(
    N_subjects = uniqueN(ID),
    N_AA = sum(AA_only == 1, na.rm = TRUE),
    N_control = sum(AA_only == 0, na.rm = TRUE),
    N_families = uniqueN(FID)
  ),
  by = .(
    splitID,
    set
  )
]

print(split_summary)

###############################################################################
# Save
###############################################################################

fwrite(
  split_long,
  file.path(
    PATH_out,
    "outer_family_membership.csv"
  )
)

fwrite(
  split_summary,
  file.path(
    PATH_out,
    "outer_family_split_summary.csv"
  )
)

cat("\nSaved family-aware outer splits.\n")

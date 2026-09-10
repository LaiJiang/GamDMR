#this code collect the GAMDMR regional results with the first selection criteria. to see how to 
#efficiently narrow down regions to run.


#in this file, we collect the DMR regions from the GAM-DMR results for AA.
#and cross-check if these signals overlap with GWAS signals for AA.!!!!


##############################################################################
#first load the DMR regions and location information.

library(data.table)
library(dplyr)
#now extract the genes
PATH_meth <- "C:/Per/LaiJiang/Project/UQAC/meth/"
PATH_scr11 <- "C:/Per/LaiJiang/Project/UQAC/meth/scr/11_mgcv/"

PATH_save <- paste0(PATH_meth,"results/15_revision/5_cv/")


setwd(PATH_save)

#load the region information. 

load( file = paste0(PATH_meth,"/scr/14_paper/results/4_overlap_M2.RData"),verbose=TRUE)
#dmr_gene_file is the object needed here. 

region_file <- fread(paste0(PATH_scr11, "dat/region_file_1_chunk.csv"))
region_file$region_id <- 1:nrow(region_file)
# Merge with results for only intersecting region_ids
dup_cols <- setdiff(intersect(names(dmr_gene_file),
                              names(region_file)),
                    "region_id")

DMRs <- merge(
  dmr_gene_file[, !dup_cols, with = FALSE],
  region_file,
  by = "region_id"
)

load( paste0(PATH_meth,"/scr/12_cpg_DMR/results/5_eval.RData"), verbose=TRUE)


#now collect the CpG M2 pvalues and evaluat them comparing to all CpGs
library(data.table)

# Keep original objects unchanged
DMRs_dt <- as.data.table(copy(DMRs))
M1_dt   <- as.data.table(copy(M1_manhattan))
M2_dt   <- as.data.table(copy(M2_manhattan))


# Extract chromosome and first CpG coordinate
# Example: "1:10563-10564" -> chr = 1, cpg_pos = 10563
prepare_CpG_data <- function(dat) {

  dat[, c("chr", "cpg_pos") := tstrsplit(
    CpG,
    split = "[:-]",
    type.convert = TRUE
  )[1:2]]

  dat[]
}


M1_dt <- prepare_CpG_data(M1_dt)
M2_dt <- prepare_CpG_data(M2_dt)


# Function to calculate the minimum pval_old within every DMR
get_DMR_min_pval <- function(CpG_data, DMR_data, output_name) {

  matched_CpGs <- CpG_data[
    DMR_data,
    on = .(
      chr,
      cpg_pos >= region_start,
      cpg_pos <= region_end
    ),
    nomatch = 0L,
    allow.cartesian = TRUE,
    .(
      region_id = i.region_id,
      pval_old
    )
  ]

  min_pvals <- matched_CpGs[
    ,
    .(
      min_pval = if (all(is.na(pval_old))) {
        NA_real_
      } else {
        min(pval_old, na.rm = TRUE)
      }
    ),
    by = region_id
  ]

  setnames(
    min_pvals,
    old = "min_pval",
    new = output_name
  )

  min_pvals
}


# Calculate minimum p-values separately for M1 and M2
M1_min_pvals <- get_DMR_min_pval(
  CpG_data = M1_dt,
  DMR_data = DMRs_dt,
  output_name = "min_pval_M1"
)

M2_min_pvals <- get_DMR_min_pval(
  CpG_data = M2_dt,
  DMR_data = DMRs_dt,
  output_name = "min_pval_M2"
)


# Remove the old column if it already exists
DMRs_dt[, min_pval_old := NULL]


# Attach M1 minimum p-values
DMRs_dt[
  M1_min_pvals,
  on = .(region_id),
  min_pval_M1 := i.min_pval_M1
]


# Attach M2 minimum p-values
DMRs_dt[
  M2_min_pvals,
  on = .(region_id),
  min_pval_M2 := i.min_pval_M2
]


# Replace the original DMRs object
DMRs <- DMRs_dt

head(DMRs)

intersect(which(is.na(DMRs$min_pval_M1)),
          which(is.na(DMRs$min_pval_M2)))


#replace NA value in DMRs$min_pval_M1 to 1, and NA value in DMRs$min_pval_M2 to 1, to avoid NA values in the next steps.
DMRs$min_pval_M1[is.na(DMRs$min_pval_M1)] <- 0
DMRs$min_pval_M2[is.na(DMRs$min_pval_M2)] <- 0


plot(DMRs$min_pval_M2,DMRs$min_pval_M1)
##########################################################
#now estimate the boundary of p1 and p2s from the DMRs$min_pval_M1 and DMRs$min_pval_M2, to see if we can find a reasonable boundary to select the DMRs regions for further analysis.

library(data.table)

dat <- as.data.table(DMRs)[
  is.finite(min_pval_M2) &
  is.finite(min_pval_M1) &
  min_pval_M2 > 0 &
  min_pval_M2 < 0.05 &
  min_pval_M1 > 0,
  .(
    x = min_pval_M2,
    y = min_pval_M1
  )
]

x_lower <- 0
x_upper <- 0.05
x_mid   <- (x_lower + x_upper) / 2

# For a proposed slope m, the smallest intercept that contains
# every point is max(y - m*x).
objective <- function(m) {

  intercept <- max(dat$y - m * dat$x)

  # Mean height of the line across x = 0 to x = 0.05
  intercept + m * x_mid
}

# Restrict to a downward-sloping upper boundary
fit <- optimize(
  f = objective,
  interval = c(-1000, 0),
  tol = 1e-12
)

slope <- fit$minimum
intercept <- max(dat$y - slope * dat$x)

# y <= intercept + slope*x
# Rearranged: -slope*x + y - intercept <= 0
a <- -slope
b <- 1
c <- -intercept

result <- data.table(
  a = a,
  b = b,
  c = c,
  slope = slope,
  intercept = intercept
)

print(result)
###################################################################
library(data.table)

# Ensure data.table format
M1_dt <- as.data.table(copy(M1_dt))
M2_dt <- as.data.table(copy(M2_dt))

###############################################################################
# 1. Create M12_dt containing M1 and M2 p-values
#
# This uses pval_old from each dataset.
###############################################################################

M1_for_merge <- M1_dt[
  ,
  .(
    CpG,
    chr,
    cpg_pos,
    M1_pval = pval_old
  )
]

M2_for_merge <- M2_dt[
  ,
  .(
    CpG,
    chr,
    cpg_pos,
    M2_pval = pval_old
  )
]

# Full join retains CpGs present in either M1 or M2
M12_dt <- merge(
  M1_for_merge,
  M2_for_merge,
  by = c("CpG", "chr", "cpg_pos"),
  all = TRUE
)

dim(M12_dt)
head(M12_dt)


###############################################################################
# 2. Extract the boundary coefficients previously stored in result
#
# Constraint:
#
# a * M2_pval + b * M1_pval + c <= 0
###############################################################################

a_boundary <- result$a[1]
b_boundary <- result$b[1]
c_boundary <- result$c[1]

cat(
  "Boundary:",
  a_boundary, "* M2_pval +",
  b_boundary, "* M1_pval +",
  c_boundary, "<= 0\n"
)


###############################################################################
# 3. Select CpGs inside the constrained region
###############################################################################

M12_selected <- M12_dt[
  is.finite(M1_pval) &
  is.finite(M2_pval) &
  M1_pval > 0 &
  M2_pval > 0 &
  M2_pval < 0.05 &
  a_boundary * M2_pval +
    b_boundary * M1_pval +
    c_boundary <= 0
]

# Add the constraint value for checking
M12_selected[
  ,
  boundary_value :=
    a_boundary * M2_pval +
    b_boundary * M1_pval +
    c_boundary
]

# More negative values lie farther inside the boundary
setorder(M12_selected, boundary_value)

dim(M12_selected)
head(M12_selected)

###################################################################
#now we know every DMRs region have minimum pvalue < 0.05. 
#now for experimental analysis, we collect all cpgs from M2_manhattan that has pval_old < 0.05, and then extract the regional information.


M12_sig_cpgs<- M12_selected


library(data.table)

# Ensure both objects are data.tables
region_dt <- as.data.table(copy(region_file))
M12_sig_cpgs_dt     <- as.data.table(copy(M12_sig_cpgs))

# Find regions containing at least one CpG from M2_sig_cpgs
overlapping_region_ids <- M12_sig_cpgs_dt[
  region_dt,
  on = .(
    chr,
    cpg_pos >= region_start,
    cpg_pos <= region_end
  ),
  nomatch = 0L,
  allow.cartesian = TRUE,
  .(region_id = i.region_id)
][
  , unique(region_id)
]

# Retain each overlapping region exactly once
region_M12 <- region_dt[
  region_id %in% overlapping_region_ids
]

# Preserve original region_file order
setorder(region_M12, data_chunk_id, region_id)

dim(region_M12)
head(region_M12)
#region_M12: inlucde all rows in region_file that have at least one CpG significant in M12_sig_cpgs.
#######################
#now extrat rows in region_file that have at least one CpG significant in GAM results DMRs.

library(data.table)

# Keep original objects unchanged
region_dt <- as.data.table(copy(region_file))
DMR_dt    <- as.data.table(copy(DMRs))

# Preserve the original order of region_file
region_dt[, region_row_id := .I]

# Prepare DMR intervals
DMR_intervals <- DMR_dt[
  ,
  .(
    DMR_region_id = region_id,
    chr,
    DMR_start = region_start,
    DMR_end = region_end
  )
]

# Set the interval table key required by foverlaps()
setkey(
  DMR_intervals,
  chr,
  DMR_start,
  DMR_end
)

# Identify every region_file row overlapping at least one DMR
region_DMR_overlaps <- foverlaps(
  x = region_dt,
  y = DMR_intervals,
  by.x = c("chr", "region_start", "region_end"),
  by.y = c("chr", "DMR_start", "DMR_end"),
  type = "any",
  nomatch = 0L
)

# Extract each region_file region once
overlapping_rows <- unique(region_DMR_overlaps$region_row_id)

region_GAM <- region_dt[
  region_row_id %in% overlapping_rows
][
  order(region_row_id)
]

# Remove temporary row identifier
region_GAM[, region_row_id := NULL]

dim(region_GAM)
head(region_GAM)

#################
#rbind region_M12 and region_GAM, and remove duplicates.

region_GAM_M12 <- rbind(region_M12, region_GAM)
region_GAM_M12 <- unique(region_GAM_M12)


#now save region_GAM_M12 to a csv file.
fwrite(region_GAM_M12, file = paste0(PATH_save,"region_GAM_M12.csv"), row.names = FALSE)


######################################
#now create 100 train_test split of the sbujects in the pheno_file. By assign 80% of FIDs into training and 20% into test.

PATH_scr11 <- "C:/Per/LaiJiang/Project/UQAC/meth/scr/11_mgcv/"
  load(file = paste0(PATH_scr11, "dat/18_pheno_BMI.RData"), verbose = TRUE)

library(data.table)

###############################################################################
# Generate 100 family-level train/test splits
# Each split contains exactly 135 unique training FIDs
###############################################################################

# Extract unique, non-missing family IDs
all_FIDs <- unique(as.character(pheno_file$FID))
all_FIDs <- all_FIDs[!is.na(all_FIDs)]

n_FIDs  <- length(all_FIDs)
n_train <- 135L
n_split <- 100L

cat("Total unique FIDs:", n_FIDs, "\n")
cat("Training FIDs per split:", n_train, "\n")
cat("Test FIDs per split:", n_FIDs - n_train, "\n")

stopifnot(n_FIDs == 169L)
stopifnot(n_train < n_FIDs)

# Reproducible random sampling
set.seed(20260727)

train_FID_list <- lapply(
  seq_len(n_split),
  function(split_id) {
    sort(sample(
      x = all_FIDs,
      size = n_train,
      replace = FALSE
    ))
  }
)

# Convert the list into a wide matrix:
# one row per split and one column per selected training FID
train_FID_matrix <- do.call(
  rbind,
  train_FID_list
)

# Convert to data.table
train_FID_df <- as.data.table(train_FID_matrix)

# Name columns train_FID_1 through train_FID_135
setnames(
  train_FID_df,
  paste0("train_FID_", seq_len(n_train))
)

# Add split identifier
train_FID_df[
  ,
  splitID := seq_len(n_split)
]

# Put splitID first
setcolorder(
  train_FID_df,
  c("splitID", paste0("train_FID_", seq_len(n_train)))
)

dim(train_FID_df)
head(train_FID_df)



fwrite(train_FID_df, file = paste0(PATH_save,"100_family_training_splits.csv"), row.names = FALSE)



#now upload the 100_family_training_splits.csv to the cluster
#  scp /mnt/c/Per/LaiJiang/Project/UQAC/meth/results/15_revision/5_cv/100_family_training_splits.csv laj773@rorqual.calculquebec.ca:/home/laj773/scratch/UQAC/meth/scr/11_mgcv/dat/100_family_training_splits.csv 

#now upload the region_GAM_M12.csv to the cluster
#  scp /mnt/c/Per/LaiJiang/Project/UQAC/meth/results/15_revision/5_cv/region_GAM_M12.csv laj773@rorqual.calculquebec.ca:/home/laj773/scratch/UQAC/meth/scr/11_mgcv/dat/region_GAM_M12.csv
#in this file, we first prepare the data will be used for DMRcate to run.
#find the region we want to run DMRcate on.   

#we re-design the DMRcate runs for looped runs.

#staeg A:
#select regions of interest -> convert each corresponding data chunk into DMRcate long format once -> split into regions with sominbus default size -> save bundled region data plus a manifest

#stage B:
#Run binomRegMethModel() on batches of regions
#Checkpoint results by region
#Record errors without terminating the job


#(Optional) Stage C:
#Rerun only failed region IDs

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

#first pre-select the regions of interest. 



##############################################################################
#first load the DMR regions and location information.

library(data.table)
library(dplyr)
#now extract the genes
PATH_meth <- "C:/Per/LaiJiang/Project/UQAC/meth/"
PATH_scr11 <- "C:/Per/LaiJiang/Project/UQAC/meth/scr/11_mgcv/"
PATH_scr14 <- paste0(PATH_meth, "scr/14_paper/")

PATH_save <- paste0(PATH_meth,"results/15_revision/5_cv/")


setwd(PATH_save)

#load the region information. 

DMRs <- fread(paste0(PATH_scr14,"/results/7_DMRcate_DMRs.csv"))


region_file <- fread(paste0(PATH_scr11, "dat/region_file_1_chunk.csv"))
region_file$region_id <- 1:nrow(region_file)
# Merge with results for only intersecting region_ids



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
      region_id = region_id,
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
DMRs$min_pval_M1[is.na(DMRs$min_pval_M1)] <- 1
DMRs$min_pval_M2[is.na(DMRs$min_pval_M2)] <- 1


plot(DMRs$min_pval_M2,DMRs$min_pval_M1,xlab="min_pval_M2",ylab="min_pval_M1")

plot(DMRs$min_pval_M2,DMRs$min_pval_M1,xlim=c(0,0.1),ylim=c(0,1),xlab="min_pval_M2",ylab="min_pval_M1")
##########################################################
#now find the 

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

#replace NA values in M12_dt$M1_pval and M12_dt$M2_pval with 0, to avoid NA values in the next steps.
M12_dt$M1_pval[is.na(M12_dt$M1_pval)] <- 1
M12_dt$M2_pval[is.na(M12_dt$M2_pval)] <- 1

#now select 


M12_selected <- M12_dt[
  M1_pval <= 0.1&
  M2_pval <= 0.05 ]




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
  .(region_id = region_id)
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
fwrite(region_GAM_M12, file = paste0(PATH_save,"region_DMRcate_M12.csv"), row.names = FALSE)


#upload the file to the server for DMRcate runs.
# scp region_DMRcate_M12.csv 


#  scp /mnt/c/Per/LaiJiang/Project/UQAC/meth/results/15_revision/5_cv/region_DMRcate_M12.csv laj773@rorqual.calculquebec.ca:/home/laj773/scratch/UQAC/meth/results/15_revision/5_cv/region_DMRcate_M12.csv


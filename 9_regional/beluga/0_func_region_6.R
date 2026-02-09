#this is functional code to create regions for 6_test_chunk.R
library(dplyr)


D_gap <- 10 #span 500 to the left and right of the CpG
#reuqired: overlap_cpgs
#required: meth_file

#output: meth_file

################################################################################
################################################################################

#Convert into a data.frame with ±100 bp regions
# Keep only entries with exactly one colon
split_list <- strsplit(overlap_cpgs, ":", fixed = FALSE)

# Convert to data.frame
split_mat <- do.call(rbind, split_list)
overlap_df <- data.frame(
  chr = split_mat[, 1],
  start = as.integer(split_mat[, 2]),
  stringsAsFactors = FALSE
)

# Add ±500 bp window
overlap_df$region_start <- overlap_df$start - D_gap
overlap_df$region_end   <- overlap_df$start + D_gap


#Step 2: Merge overlapping or adjacent regions
# Ensure sorted by chr and region_start
merged_regions <- overlap_df %>%
  arrange(chr, region_start) %>%
  group_by(chr) %>%
  mutate(
    new_region = cumsum(c(TRUE, diff(region_start) > D_gap*2))
  ) %>%
  group_by(chr, new_region) %>%
  summarize(
    merged_start = min(region_start),
    merged_end = max(region_end),
    .groups = "drop"
  )


#step 3: Filter CpGs in meth_file that fall within these regions
# Create lookup vector for fast matching
meth_file$coord <- paste0(meth_file$chr, ":", meth_file$start)

# Subset meth_file to only rows within any of the merged regions
selected_cpgs <- purrr::map_dfr(seq_len(nrow(merged_regions)), function(i) {
  chr_i <- merged_regions$chr[i]
  start_i <- merged_regions$merged_start[i]
  end_i <- merged_regions$merged_end[i]
  
  meth_file %>%
    filter(chr == chr_i, start >= start_i, start <= end_i)
})

print("number of selected CpGs:")
print(dim(selected_cpgs))
print("number of regions:")
print(nrow(merged_regions))
#approximiately 10 cpgs per region.

#remove coord column
selected_cpgs <- selected_cpgs %>%
  select(-coord)

meth_file <- selected_cpgs


rm(selected_cpgs, overlap_df, merged_regions, split_list, split_mat)
rm(overlap_cpgs)
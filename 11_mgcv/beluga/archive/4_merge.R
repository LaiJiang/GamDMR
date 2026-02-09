#now merger all region data from 3_spacing.R

#merge closely spaced cpg regions into one region.

#then split region with more than 1000 cpgs into smaller regions for computational efficiency.

#note also need to add the chunk ID into the final output file. so we can easily trace back to the original data chunk. for next model running analysis.



# Load necessary libraries
library(dplyr)
library(data.table)
library(stringr)

# Set path
PATH_wk <- "~/scratch/UQAC/meth/"
input_dir <- paste0(PATH_wk, "data/region_spacing/")
max_cpgs <- 100
gap <- 250

# Step 1: Read and merge all region files
region_files <- list.files(input_dir, pattern = "spacing_\\d{4}\\.csv", full.names = TRUE)

merged_regions <- lapply(region_files, function(f) {
  df <- fread(f, header = FALSE, col.names = c("chr", "region_start", "region_end", "n_cpgs"))
  df$data_chunk_id <- as.character(str_extract(basename(f), "\\d{4}"))
  df
}) %>% bind_rows()

# Step 2: Sort regions by chr and start
merged_regions <- merged_regions %>% arrange(chr, region_start)

# Step 3: Merge adjacent regions if gap < 250
merged_list <- list()
current_region <- merged_regions[1, ]

for (i in 2:nrow(merged_regions)) {
  row <- merged_regions[i, ]
  if (row$chr == current_region$chr &&
      row$region_start - current_region$region_end < gap) {
    # ensure both are character before concatenation
    current_region$data_chunk_id <- as.character(current_region$data_chunk_id)
    row$data_chunk_id <- as.character(row$data_chunk_id)
    # merge regions
    current_region$region_end <- row$region_end
    current_region$n_cpgs <- current_region$n_cpgs + row$n_cpgs
    current_region$data_chunk_id <- paste(current_region$data_chunk_id, row$data_chunk_id, sep = ",")
  } else {
    merged_list[[length(merged_list) + 1]] <- current_region
    current_region <- row
  }
}
# Add last region
merged_list[[length(merged_list) + 1]] <- current_region

# Ensure consistent type before binding
merged_list <- lapply(merged_list, function(x) {
  x$data_chunk_id <- as.character(x$data_chunk_id)
  x
})
merged_df <- bind_rows(merged_list)

# Step 4: Split regions with n_cpgs > max_cpgs
final_list <- list()
for (i in 1:nrow(merged_df)) {
  row <- merged_df[i, ]
  if (row$n_cpgs <= max_cpgs) {
    final_list[[length(final_list) + 1]] <- row
  } else {
    # evenly split position range
    n_splits <- ceiling(row$n_cpgs / max_cpgs)
    region_size <- ceiling((row$region_end - row$region_start + 1) / n_splits)
    for (j in 0:(n_splits - 1)) {
      sub_start <- row$region_start + j * region_size
      sub_end <- min(row$region_end, sub_start + region_size - 1)
      sub_n_cpgs <- min(max_cpgs, row$n_cpgs - j * max_cpgs)
      final_list[[length(final_list) + 1]] <- data.frame(
        data_chunk_id = row$data_chunk_id,
        chr = row$chr,
        region_start = sub_start,
        region_end = sub_end,
        n_cpgs = sub_n_cpgs
      )
    }
  }
}

final_df <- bind_rows(final_list)

# Reorder columns
final_df <- final_df %>%
  select(data_chunk_id, chr, region_start, region_end, n_cpgs)

# View a sample

# View a sample
print(head(final_df, 10))

#save the final merged regions
fwrite(final_df, paste0(PATH_wk, "results/11_mgcv/merged_regions.csv"), row.names = FALSE)
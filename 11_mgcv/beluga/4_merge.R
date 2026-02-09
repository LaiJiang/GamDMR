
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
    current_region$data_chunk_id <- as.character(current_region$data_chunk_id)
    row$data_chunk_id <- as.character(row$data_chunk_id)
    current_region$region_end <- row$region_end
    current_region$n_cpgs <- current_region$n_cpgs + row$n_cpgs
    current_region$data_chunk_id <- paste(current_region$data_chunk_id, row$data_chunk_id, sep = ",")
  } else {
    merged_list[[length(merged_list) + 1]] <- current_region
    current_region <- row
  }
}
merged_list[[length(merged_list) + 1]] <- current_region

merged_list <- lapply(merged_list, function(x) {
  x$data_chunk_id <- as.character(x$data_chunk_id)
  x
})
merged_df <- bind_rows(merged_list)

# Step 4: Corrected CpG-aware splitting for regions with > max_cpgs
final_list <- list()
for (i in 1:nrow(merged_df)) {
  row <- merged_df[i, ]

  # Simulate CpG positions evenly if actual positions are not available
  region_cpg_positions <- seq(row$region_start, row$region_end, length.out = row$n_cpgs)

  # Split into chunks of max_cpgs
  split_indices <- split(region_cpg_positions, ceiling(seq_along(region_cpg_positions) / max_cpgs))

  for (sub_pos in split_indices) {
    final_list[[length(final_list) + 1]] <- data.frame(
      data_chunk_id = row$data_chunk_id,
      chr = row$chr,
      region_start = floor(min(sub_pos)),
      region_end   = ceiling(max(sub_pos)),
      n_cpgs = length(sub_pos)
    )
  }
}

final_df <- bind_rows(final_list)

# Reorder columns
final_df <- final_df %>%
  select(data_chunk_id, chr, region_start, region_end, n_cpgs)

# View sample
print(head(final_df, 10))

# Save output
fwrite(final_df, paste0(PATH_wk, "results/11_mgcv/merged_regions.csv"), row.names = FALSE)

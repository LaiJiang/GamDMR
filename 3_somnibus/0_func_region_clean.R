suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

window_size <- 10

if (!exists("overlap_cpgs")) stop("object 'overlap_cpgs' not found")
if (!exists("meth_file")) stop("object 'meth_file' not found")

split_list <- strsplit(overlap_cpgs, ":", fixed = TRUE)
split_lengths <- lengths(split_list)
if (any(split_lengths != 2)) {
  stop("Each entry in overlap_cpgs must have the form 'chr:start'.")
}

split_mat <- do.call(rbind, split_list)
overlap_df <- data.frame(
  chr = split_mat[, 1],
  start = as.integer(split_mat[, 2]),
  stringsAsFactors = FALSE
)

if (anyNA(overlap_df$start)) {
  stop("Failed to parse start positions in overlap_cpgs.")
}

overlap_df$region_start <- overlap_df$start - window_size
overlap_df$region_end <- overlap_df$start + window_size

merged_regions <- overlap_df %>%
  arrange(chr, region_start) %>%
  group_by(chr) %>%
  mutate(new_region = cumsum(c(TRUE, diff(region_start) > window_size * 2))) %>%
  group_by(chr, new_region) %>%
  summarize(
    merged_start = min(region_start),
    merged_end = max(region_end),
    .groups = "drop"
  )

selected_cpgs <- purrr::map_dfr(seq_len(nrow(merged_regions)), function(i) {
  chr_i <- merged_regions$chr[i]
  start_i <- merged_regions$merged_start[i]
  end_i <- merged_regions$merged_end[i]

  meth_file %>%
    filter(chr == chr_i, start >= start_i, start <= end_i)
})

meth_file <- selected_cpgs

rm(selected_cpgs, overlap_df, merged_regions, split_list, split_mat)

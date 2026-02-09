#note max_cpgs = 2000 is important, otherwise the region will be too small to run mgcv.
splitCpGsBySpacing <- function(df, gap = 250, min_cpgs = 10, max_cpgs = 2000) {
  if (!all(c("chr", "position") %in% colnames(df))) {
    stop("Input must have columns 'chr' and 'position'")
  }

  df <- df[order(df$chr, df$position), ]
  result <- list()

  for (ch in unique(df$chr)) {
    chr_df <- df[df$chr == ch, ]
    pos <- chr_df$position
    gaps <- c(Inf, diff(pos))
    region_id <- cumsum(gaps > gap)

    region_list <- split(pos, region_id)

    for (region in region_list) {
      n_cpg <- length(region)

      if (n_cpg < min_cpgs) next  # skip too-small regions

      if (n_cpg <= max_cpgs) {
        result[[length(result) + 1]] <- data.frame(
          chr = ch,
          region_start = min(region),
          region_end = max(region),
          n_cpgs = n_cpg
        )
      } else {
        # Split large region into chunks of max_cpgs
        split_indices <- split(region, ceiling(seq_along(region) / max_cpgs))
        for (sub_region in split_indices) {
          result[[length(result) + 1]] <- data.frame(
            chr = ch,
            region_start = min(sub_region),
            region_end = max(sub_region),
            n_cpgs = length(sub_region)
          )
        }
      }
    }
  }

  do.call(rbind, result)
}

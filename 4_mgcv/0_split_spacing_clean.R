split_cpgs_by_spacing <- function(df, gap = 1000, min_cpgs = 10, max_cpgs = 2000) {
  required_cols <- c("chr", "position")
  if (!all(required_cols %in% names(df))) {
    stop("Input must contain columns: ", paste(required_cols, collapse = ", "))
  }

  df <- df[order(df$chr, df$position), , drop = FALSE]
  result <- list()

  for (ch in unique(df$chr)) {
    chr_df <- df[df$chr == ch, , drop = FALSE]
    pos <- chr_df$position
    gaps <- c(Inf, diff(pos))
    region_id <- cumsum(gaps > gap)
    region_list <- split(pos, region_id)

    for (region in region_list) {
      n_cpg <- length(region)

      if (n_cpg < min_cpgs) next

      if (n_cpg <= max_cpgs) {
        result[[length(result) + 1]] <- data.frame(
          chr = ch,
          region_start = min(region),
          region_end = max(region),
          n_cpgs = n_cpg
        )
      } else {
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

  if (length(result) == 0) {
    return(data.frame(
      chr = character(),
      region_start = numeric(),
      region_end = numeric(),
      n_cpgs = numeric()
    ))
  }

  do.call(rbind, result)
}

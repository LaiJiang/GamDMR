convert_to_somnibus_format <- function(meth_file) {
  position <- meth_file$start
  col_names <- colnames(meth_file)

  meth_cols <- grep("_meth$", col_names)
  tot_cols <- grep("_tot$", col_names)

  if (length(meth_cols) == 0 || length(tot_cols) == 0) {
    stop("Methylation and total-count columns were not found.")
  }

  meth_sample_ids <- sub("_meth$", "", col_names[meth_cols])
  tot_sample_ids <- sub("_tot$", "", col_names[tot_cols])

  if (!identical(sort(meth_sample_ids), sort(tot_sample_ids))) {
    stop("Mismatch between methylation and total-count sample IDs.")
  }

  common_sample_ids <- meth_sample_ids

  long_list <- vector("list", length(common_sample_ids))

  for (i in seq_along(common_sample_ids)) {
    sample_id <- common_sample_ids[i]
    meth_col <- paste0(sample_id, "_meth")
    tot_col <- paste0(sample_id, "_tot")

    long_list[[i]] <- data.frame(
      Meth_Counts = meth_file[[meth_col]],
      Total_Counts = meth_file[[tot_col]],
      Position = position,
      ID = sample_id,
      stringsAsFactors = FALSE
    )
  }

  somnibus_input <- do.call(rbind, long_list)
  somnibus_input$ID <- gsub("\\.", "-", somnibus_input$ID)
  somnibus_input$ID <- gsub("^X", "", somnibus_input$ID)

  somnibus_input
}

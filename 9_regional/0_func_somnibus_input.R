convert_to_somnibus_format <- function(meth_file) {
  # Extract position
  position <- meth_file$start
  
  # Find the indices of meth and total count columns
  col_names <- colnames(meth_file)
  meth_cols <- grep("_meth$", col_names)
  tot_cols <- grep("_tot$", col_names)

  # Check that each meth has a corresponding tot column
  if (length(meth_cols) != length(tot_cols)) {
    stop("Mismatch between number of methylated and total count columns.")
  }
  
  # Extract sample IDs from column names
  sample_ids <- sub("_meth$", "", col_names[meth_cols])

  # Initialize list to store per-sample data frames
  long_list <- vector("list", length(sample_ids))
  
  for (i in seq_along(sample_ids)) {
    sample_id <- sample_ids[i]
    meth_counts <- meth_file[[paste0(sample_id, "_meth")]]
    total_counts <- meth_file[[paste0(sample_id, "_tot")]]
    
    long_list[[i]] <- data.frame(
      Meth_Counts = meth_counts,
      Total_Counts = total_counts,
      Position = position,
      ID = sample_id,
      stringsAsFactors = FALSE
    )
  }

  # Combine all sample data into a single data frame
  somnibus_input <- do.call(rbind, long_list)


  somnibus_input$ID <- gsub("\\.", "-", somnibus_input$ID)  # Replace the dot with a dash
  somnibus_input$ID <- gsub("^X", "", somnibus_input$ID)    # Remove the leading "X"
  
  return(somnibus_input)
}

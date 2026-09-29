###############################################################################
# 3_collect_mgcv.R
###############################################################################

library(data.table)

PATH_wk <- path.expand(
  "~/scratch/UQAC/meth/"
)

PATH_validation <- file.path(
  PATH_wk,
  "results/15_revision/5_cv"
)

PATH_raw <- file.path(
  PATH_validation,
  "1_mgcv_raw"
)

PATH_selected <- file.path(
  PATH_validation,
  "2_mgcv_selected"
)

dir.create(
  PATH_selected,
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################################
###############################################################################

P_THRESHOLD <- 1e-5

FDR_THRESHOLD <- 0.05

EDF_THRESHOLD <- 0.5

R2_THRESHOLD <- 0.50

MEAN_DIFF_THRESHOLD <- 0.05

###############################################################################
# Process every outer split independently
###############################################################################

summary_list <- list()

for (
  split_id in 1:100
) {

  cat(
    "\nProcessing split",
    split_id,
    "\n"
  )

  files <- list.files(
    PATH_raw,
    pattern = sprintf(
      "^gam_split_%03d_job_[0-9]+\\.tsv$",
      split_id
    ),
    full.names = TRUE
  )

  if (
    length(files) == 0L
  ) {

    warning(
      "No files found for split ",
      split_id
    )

    next
  }

  x <- rbindlist(
    lapply(
      files,
      fread
    ),
    fill = TRUE
  )

  ###########################################################################
  # Make sure no region was accidentally duplicated
  ###########################################################################

  if (
    anyDuplicated(
      x$region_row_id
    )
  ) {

    stop(
      "Duplicated regions found for split ",
      split_id
    )
  }

  ###########################################################################
  # Only successfully fitted regions have a valid GAM test
  ###########################################################################

  x[
    ,
    FDR := NA_real_
  ]

  valid <- which(
    x$status == "ok" &
    is.finite(x$p_smooth_AA)
  )

  ###########################################################################
  # BH correction across ALL successfully tested regions in THIS training set
  ###########################################################################

  x[
    valid,
    FDR :=
      p.adjust(
        p_smooth_AA,
        method = "BH"
      )
  ]

  ###########################################################################
  # Strict GAM-DMR definition
  ###########################################################################

  x[
    ,
    GAM_DMR :=
      status == "ok" &

      is.finite(p_smooth_AA) &
      p_smooth_AA < P_THRESHOLD &

      is.finite(FDR) &
      FDR < FDR_THRESHOLD &

      is.finite(edf_smooth_AA) &
      edf_smooth_AA >= EDF_THRESHOLD &

      is.finite(R2) &
      R2 >= R2_THRESHOLD &

      is.finite(mean_diff) &
      mean_diff >= MEAN_DIFF_THRESHOLD
  ]

  ###########################################################################
  # Save complete training-only regional results
  ###########################################################################

  fwrite(
    x,
    file.path(
      PATH_selected,
      sprintf(
        "all_GAM_results_split_%03d.tsv",
        split_id
      )
    ),
    sep = "\t"
  )

  ###########################################################################
  # Save only DMRs selected within this split
  ###########################################################################

  dmr <- x[
    GAM_DMR == TRUE
  ]

  fwrite(
    dmr,
    file.path(
      PATH_selected,
      sprintf(
        "GAM_DMR_split_%03d.tsv",
        split_id
      )
    ),
    sep = "\t"
  )

  ###########################################################################
  # Diagnostics
  ###########################################################################

  summary_list[[split_id]] <-
    data.table(

      splitID =
        split_id,

      N_region_rows =
        nrow(x),

      N_successful_GAM =
        sum(
          x$status == "ok"
        ),

      N_failed_GAM =
        sum(
          x$status == "gam_error"
        ),

      N_DMR =
        nrow(dmr)
    )

  cat(
    "Successfully fitted:",
    sum(x$status == "ok"),
    "\n"
  )

  cat(
    "Selected GAM-DMRs:",
    nrow(dmr),
    "\n"
  )
}

###############################################################################
# Save summary
###############################################################################

summary_table <- rbindlist(
  summary_list,
  fill = TRUE
)

fwrite(
  summary_table,
  file.path(
    PATH_validation,
    "GAM_DMR_split_summary.tsv"
  ),
  sep = "\t"
)

print(
  summary_table
)

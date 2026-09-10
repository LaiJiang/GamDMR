#!/usr/bin/env Rscript

# =============================================================================
# 20_generate_Figure5_AUROC_PRAUC.R
#
# Generates Figure 5 from the already-combined all-method split-level
# performance results. This figure is identical in style to
# Figure1_all_methods_split_performance_distributions, but includes only:
#   1) AUROC
#   2) PR-AUC
#
# Brier score is excluded.
#
# Outputs:
#   Figure5.pdf
#   Figure5.png
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

data.table::setDTthreads(1L)

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

PATH_wk <- path.expand(Sys.getenv("PATH_WK", "~/scratch/UQAC/meth/"))

pred_root <- file.path(
  PATH_wk,
  "results/15_revision/5_cv/5_prediction"
)

out_root <- path.expand(Sys.getenv(
  "ALL_METHOD_COMPARISON_OUTPUT",
  file.path(pred_root, "7_all_method_comparison")
))

dat_dir <- file.path(out_root, "data")
fig_dir <- file.path(out_root, "figures")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

perf_file <- file.path(
  dat_dir,
  "all_methods_primary_performance_all_splits.tsv"
)

if (!file.exists(perf_file)) {
  stop(
    "Combined split-level performance file not found:\n",
    perf_file,
    "\nRun 19_combine_all_primary_prediction_methods.R first."
  )
}

# -------------------------------------------------------------------------
# Method display settings: identical to the original all-method figure
# -------------------------------------------------------------------------

method_order <- c(
  "GAM-DMR",
  "SOMNiBUS",
  "DMRcate",
  "BSmooth",
  "Univariate method M1",
  "Univariate method M2",
  "Univariate method M3",
  "Covariates only"
)

method_palette <- c(
  "GAM-DMR" = "#D55E00",
  "SOMNiBUS" = "#0072B2",
  "DMRcate" = "#009E73",
  "BSmooth" = "#CC79A7",
  "Univariate method M1" = "#E69F00",
  "Univariate method M2" = "#56B4E9",
  "Univariate method M3" = "#000000",
  "Covariates only" = "#777777"
)

# -------------------------------------------------------------------------
# Load split-level performance
# -------------------------------------------------------------------------

perf <- fread(perf_file)

required_cols <- c(
  "splitID",
  "benchmark_method",
  "AUROC_test",
  "PR_AUC_test"
)

missing_cols <- setdiff(required_cols, names(perf))

if (length(missing_cols) > 0L) {
  stop(
    "Missing required columns in performance file: ",
    paste(missing_cols, collapse = ", ")
  )
}

perf[, benchmark_method := factor(
  as.character(benchmark_method),
  levels = method_order
)]

# -------------------------------------------------------------------------
# Long format: AUROC and PR-AUC only
# -------------------------------------------------------------------------

pv <- rbindlist(list(
  perf[, .(
    splitID,
    benchmark_method,
    metric = "AUROC",
    value = AUROC_test
  )],
  perf[, .(
    splitID,
    benchmark_method,
    metric = "PR-AUC",
    value = PR_AUC_test
  )]
))

pv[, metric := factor(
  metric,
  levels = c("AUROC", "PR-AUC")
)]

# Mean for the white diamond used in the original figure
means <- pv[
  ,
  .(value = mean(value, na.rm = TRUE)),
  by = .(benchmark_method, metric)
]

# -------------------------------------------------------------------------
# Figure 5
# -------------------------------------------------------------------------

p5 <- ggplot(
  pv,
  aes(
    x = benchmark_method,
    y = value,
    fill = benchmark_method
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.45,
    linewidth = 0.4
  ) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    linewidth = 0.4
  ) +
  geom_point(
    data = means,
    shape = 23,
    size = 2.4,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    nrow = 1
  ) +
  scale_fill_manual(
    values = method_palette,
    drop = FALSE
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Predictive performance across repeated family-aware test splits",
    x = NULL,
    y = NULL
  )

# -------------------------------------------------------------------------
# Save
# -------------------------------------------------------------------------

pdf_file <- file.path(fig_dir, "Figure5.pdf")
png_file <- file.path(fig_dir, "Figure5.png")

ggsave(
  filename = pdf_file,
  plot = p5,
  width = 9.5,
  height = 5.5
)

ggsave(
  filename = png_file,
  plot = p5,
  width = 9.5,
  height = 5.5,
  dpi = 300
)

cat("\nFigure 5 generated successfully.\n")
cat("PDF: ", pdf_file, "\n", sep = "")
cat("PNG: ", png_file, "\n", sep = "")

cat("\nMean split-level performance shown by white diamonds:\n")
print(
  means[
    order(metric, benchmark_method)
  ]
)

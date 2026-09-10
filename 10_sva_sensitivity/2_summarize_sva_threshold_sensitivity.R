#!/usr/bin/env Rscript

# =============================================================================
# 2_summarize_sva_threshold_sensitivity.R
#
# Compare SVA solutions for vfilter = 10k, 25k, 50k.
#
# Outputs:
#   SVA_nSV_summary.csv
#   SVA_pairwise_stability.csv
#   SVA_primary_SV_match_summary.csv
#   abs_cor_SV_*.csv
#   heatmap_abs_cor_SV_*.png
#   subject_distance_*.png
#   SVA_sensitivity_summary.txt
#
# Stability metrics:
#   1. number of retained SVs
#   2. absolute cross-correlation among individual SVs
#   3. rotation-invariant similarity between entire SV subspaces
#      (singular values / principal-angle cosines)
#   4. Spearman correlation of pairwise subject distances in orthonormal SV space
#   5. cophenetic correlation of hierarchical clustering trees
#
# Metrics 3-5 remain meaningful if SV labels/signs change across thresholds.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

data.table::setDTthreads(1L)

PATH_wk <- path.expand(
  Sys.getenv("METH_BASE_DIR", unset = "~/scratch/UQAC/meth")
)

out_dir <- path.expand(
  Sys.getenv(
    "SVA_SENS_OUT",
    unset = file.path(
      PATH_wk, "results", "15_revision", "7_sva_threshold_sensitivity"
    )
  )
)

vfilters <- c(10000L, 25000L, 50000L)

files <- setNames(
  file.path(out_dir, sprintf("sva_vfilter_%05d.rds", vfilters)),
  as.character(vfilters)
)

missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing SVA results:\n",
    paste(missing_files, collapse = "\n")
  )
}

res <- lapply(files, readRDS)

# --------------------------- Basic checks -----------------------------------

ids <- lapply(res, function(x) as.character(x$subject_ids))

if (!all(vapply(ids[-1], identical, logical(1), ids[[1]]))) {
  stop("Subject order differs among SVA threshold results.")
}

nsv_summary <- rbindlist(lapply(res, function(x) {
  data.table(
    vfilter = x$vfilter,
    n_sv = x$n_sv,
    n_subjects = nrow(x$sv),
    sva_version = x$sva_version,
    seed = x$seed
  )
}))

setorder(nsv_summary, vfilter)
fwrite(nsv_summary, file.path(out_dir, "SVA_nSV_summary.csv"))

cat("Number of SVs:\n")
print(nsv_summary)

# ---------------------------- Helpers ---------------------------------------

orthonormal_basis <- function(X) {

  X <- as.matrix(X)
  X <- scale(X, center = TRUE, scale = FALSE)

  q <- qr(X)
  rank <- q$rank

  if (rank < 1L) stop("SV matrix has rank zero.")

  qr.Q(q)[, seq_len(rank), drop = FALSE]
}

pairwise_stability <- function(A, B, label_a, label_b) {

  A <- as.matrix(A)
  B <- as.matrix(B)

  # Direct individual-SV correlation matrix.
  cor_mat <- abs(cor(A, B, use = "pairwise.complete.obs"))

  # Rotation-invariant subspace similarity.
  QA <- orthonormal_basis(A)
  QB <- orthonormal_basis(B)

  svals <- svd(
    crossprod(QA, QB),
    nu = 0,
    nv = 0
  )$d

  svals <- pmin(pmax(svals, 0), 1)

  # Pairwise subject geometry in SV space.
  dist_A <- as.vector(dist(QA))
  dist_B <- as.vector(dist(QB))

  distance_rho <- suppressWarnings(
    cor(dist_A, dist_B, method = "spearman", use = "complete.obs")
  )

  # Hierarchical clustering stability.
  # Scale columns so an arbitrary SV scale does not dominate.
  scale_safe <- function(X) {
    X <- scale(X)
    X[, apply(X, 2, function(z) all(is.finite(z))), drop = FALSE]
  }

  As <- scale_safe(A)
  Bs <- scale_safe(B)

  hc_A <- hclust(dist(As), method = "ward.D2")
  hc_B <- hclust(dist(Bs), method = "ward.D2")

  cophenetic_rho <- suppressWarnings(
    cor(
      as.vector(cophenetic(hc_A)),
      as.vector(cophenetic(hc_B)),
      method = "spearman",
      use = "complete.obs"
    )
  )

  list(
    summary = data.table(
      comparison = paste0(label_a, " vs ", label_b),
      n_sv_A = ncol(A),
      n_sv_B = ncol(B),
      mean_principal_angle_cosine = mean(svals),
      min_principal_angle_cosine = min(svals),
      median_principal_angle_cosine = median(svals),
      subject_distance_spearman = distance_rho,
      clustering_cophenetic_spearman = cophenetic_rho,
      mean_max_abs_sv_correlation_A_to_B =
        mean(apply(cor_mat, 1L, max, na.rm = TRUE)),
      min_max_abs_sv_correlation_A_to_B =
        min(apply(cor_mat, 1L, max, na.rm = TRUE))
    ),
    cor_mat = cor_mat,
    QA = QA,
    QB = QB
  )
}

# --------------------------- Pairwise comparisons ---------------------------

pairs <- list(
  c("10000", "25000"),
  c("10000", "50000"),
  c("25000", "50000")
)

stability_rows <- list()

for (k in seq_along(pairs)) {

  a <- pairs[[k]][1]
  b <- pairs[[k]][2]

  A <- res[[a]]$sv
  B <- res[[b]]$sv

  comp <- pairwise_stability(
    A, B,
    paste0(as.integer(a) / 1000, "k"),
    paste0(as.integer(b) / 1000, "k")
  )

  stability_rows[[k]] <- comp$summary

  # Save the full absolute correlation matrix.
  cor_dt <- as.data.table(comp$cor_mat, keep.rownames = "SV_A")
  fwrite(
    cor_dt,
    file.path(
      out_dir,
      sprintf("abs_cor_SV_%sk_vs_%sk.csv",
              as.integer(a) / 1000,
              as.integer(b) / 1000)
    )
  )

  # Correlation heatmap.
  png(
    file.path(
      out_dir,
      sprintf("heatmap_abs_cor_SV_%sk_vs_%sk.png",
              as.integer(a) / 1000,
              as.integer(b) / 1000)
    ),
    width = 1800,
    height = 1500,
    res = 220
  )

  par(mar = c(6, 6, 3, 2))

  heatmap(
    comp$cor_mat,
    Rowv = NA,
    Colv = NA,
    scale = "none",
    margins = c(8, 8),
    xlab = paste0(as.integer(b) / 1000, "k CpGs"),
    ylab = paste0(as.integer(a) / 1000, "k CpGs"),
    main = "Absolute correlation between surrogate variables"
  )

  dev.off()

  # Pairwise subject-distance comparison.
  dA <- as.vector(dist(comp$QA))
  dB <- as.vector(dist(comp$QB))

  png(
    file.path(
      out_dir,
      sprintf("subject_distance_%sk_vs_%sk.png",
              as.integer(a) / 1000,
              as.integer(b) / 1000)
    ),
    width = 1600,
    height = 1400,
    res = 220
  )

  plot(
    dA, dB,
    pch = 16,
    cex = 0.35,
    xlab = paste0("Pairwise subject distance: ", as.integer(a) / 1000, "k"),
    ylab = paste0("Pairwise subject distance: ", as.integer(b) / 1000, "k"),
    main = paste0(
      "Latent sample structure (Spearman rho = ",
      sprintf("%.3f", comp$summary$subject_distance_spearman),
      ")"
    )
  )

  dev.off()
}

stability_table <- rbindlist(stability_rows)
fwrite(
  stability_table,
  file.path(out_dir, "SVA_pairwise_stability.csv")
)

# ------------------- Primary 10k SV matching summary ------------------------

sv10 <- as.matrix(res[["10000"]]$sv)

match_rows <- list()

for (b in c("25000", "50000")) {

  B <- as.matrix(res[[b]]$sv)
  C <- abs(cor(sv10, B))

  tmp <- rbindlist(lapply(seq_len(nrow(C)), function(i) {

    j <- which.max(C[i, ])

    data.table(
      primary_sv = rownames(C)[i],
      comparison_vfilter = as.integer(b),
      best_matching_sv = colnames(C)[j],
      abs_correlation = C[i, j]
    )
  }))

  match_rows[[b]] <- tmp
}

match_table <- rbindlist(match_rows)

fwrite(
  match_table,
  file.path(out_dir, "SVA_primary_SV_match_summary.csv")
)

# --------------------------- Human-readable summary -------------------------

summary_txt <- file.path(out_dir, "SVA_sensitivity_summary.txt")

sink(summary_txt)

cat("SVA CpG-selection threshold sensitivity analysis\n")
cat("================================================\n\n")

cat("Number of surrogate variables\n")
print(nsv_summary)
cat("\n")

cat("Pairwise latent-space stability\n")
print(stability_table)
cat("\n")

cat("Best matches for each primary 10k surrogate variable\n")
print(match_table)
cat("\n")

cat("Interpretation guide\n")
cat("--------------------\n")
cat("* Individual SV signs and ordering are arbitrary; use absolute correlations.\n")
cat("* Principal-angle cosine values close to 1 indicate nearly identical SV subspaces.\n")
cat("* Subject-distance Spearman values close to 1 indicate stable sample geometry/clustering.\n")
cat("* Cophenetic Spearman values close to 1 indicate similar hierarchical clustering patterns.\n")
cat("* Downstream GAM-DMR stability should be evaluated separately after replacing the SVs.\n")

sink()

cat("\nSensitivity summaries written to:", out_dir, "\n")
print(stability_table)

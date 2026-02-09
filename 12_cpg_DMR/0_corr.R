## ---- correlogram_pheno.R ----
## Assumes `pheno_file` is a data.frame already in memory
## Excludes ID and FID; encodes binary factors; skips multi-level non-numeric

# Choose correlation method: "pearson" (default) or "spearman"
corr_method <- "pearson"

# 1) Drop ID/FID and prep numeric matrix
df <- pheno_file
drop_cols <- c("ID", "FID")
df <- df[ , setdiff(names(df), drop_cols), drop = FALSE]

# Helper: coerce columns to numeric where safe
prep_numeric <- function(d) {
  out <- list()
  dropped <- c()
  for (nm in names(d)) {
    x <- d[[nm]]
    if (is.numeric(x)) {
      out[[nm]] <- x
    } else if (is.logical(x)) {
      out[[nm]] <- as.numeric(x)
    } else if (is.factor(x) || is.character(x)) {
      ux <- unique(na.omit(x))
      if (length(ux) == 2L) {
        # Binary factor → 0/1 (alphabetical order defines 0/1)
        f <- factor(x, levels = sort(ux))
        out[[nm]] <- as.numeric(f) - 1L
      } else {
        dropped <- c(dropped, nm)
      }
    } else {
      dropped <- c(dropped, nm)
    }
  }
  if (length(dropped)) {
    message("Dropped non-numeric/non-binary columns: ", paste(dropped, collapse = ", "))
  }
  as.data.frame(out, check.names = TRUE)
}

X <- prep_numeric(df)

# Need at least 2 numeric columns
stopifnot(ncol(X) >= 2)

# 2) Correlation and p-values
cor_mat <- cor(X, use = "pairwise.complete.obs", method = corr_method)

# p-value matrix (no extra packages required)
cor_pmat <- function(mat, method = "pearson") {
  m <- ncol(mat)
  p <- matrix(NA_real_, m, m, dimnames = list(colnames(mat), colnames(mat)))
  diag(p) <- 0
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      ct <- suppressWarnings(cor.test(mat[, i], mat[, j], method = method))
      p[i, j] <- p[j, i] <- ct$p.value
    }
  }
  p
}
p_mat <- cor_pmat(X, method = corr_method)

# 3) Save numeric outputs
#write.csv(cor_mat, file = "cor_matrix.csv", row.names = TRUE)
#write.csv(p_mat,  file = "cor_pvalues.csv", row.names = TRUE)

# 4) Plot: use corrplot if present; else base heatmap
plot_file <- file.path(PATH_results, "12_cpg_DMR", "correlogram.jpg")

# JPEG device (pixels). Tweak width/height/res/quality as you like.
jpeg(plot_file, width = 2000, height = 2000, res = 220, quality = 95)

if (requireNamespace("corrplot", quietly = TRUE)) {
  corrplot::corrplot(
    cor_mat,
    method = "color",
    type = "upper",
    order = "hclust",
    addCoef.col = "black",
    tl.col = "black",
    tl.srt = 45,
    p.mat = p_mat,
    sig.level = 0.05,
    insig = "blank",
    diag = FALSE
  )
  mtext(sprintf("Correlogram (%s) — blanked if p ≥ 0.05", corr_method), line = 0.5)
} else {
  message("Package 'corrplot' not found; using base heatmap. For nicer plot, install.packages('corrplot').")
  heatmap(cor_mat, symm = TRUE, Colv = NA, Rowv = NA, scale = "none")
  title(main = sprintf("Correlogram (%s)", corr_method))
}
dev.off()


message("Done. Wrote: cor_matrix.csv, cor_pvalues.csv, and ", plot_file)




## Age vs BMI colored by AA_only → JPEG

# Use your existing results path; fall back to cwd if not set
if (!exists("PATH_results")) PATH_results <- getwd()
out_dir   <- file.path(PATH_results, "12_cpg_DMR")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
plot_file <- file.path(out_dir, "Age_vs_BMI_by_AA.jpg")

# Build a minimal data frame
df <- pheno_file[, c("AgeCalc", "BMI", "AA_only")]
df <- df[is.finite(df$AgeCalc) & is.finite(df$BMI) & !is.na(df$AA_only), , drop = FALSE]

# Ensure AA_only is a factor with sensible labels
if (!is.factor(df$AA_only)) {
  u <- sort(unique(df$AA_only))
  if (length(u) == 2 && all(u %in% c(0, 1))) {
    df$AA_only <- factor(df$AA_only, levels = c(0, 1), labels = c("Control", "AA"))
  } else {
    df$AA_only <- factor(df$AA_only)
  }
}

# Compute Pearson r for subtitle
r_val <- suppressWarnings(cor(df$AgeCalc, df$BMI, use = "complete.obs", method = "pearson"))
n_pair <- sum(stats::complete.cases(df[, c("AgeCalc","BMI")]))
subtitle_txt <- sprintf("Pearson r = %.2f (n = %d)", r_val, n_pair)

# Plot (ggplot2 if available; else base)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  jpeg(plot_file, width = 1800, height = 1400, res = 220, quality = 95)
  p <- ggplot(df, aes(x = AgeCalc, y = BMI, color = AA_only)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = "Age vs BMI", subtitle = subtitle_txt,
         x = "Age (years)", y = "BMI (kg/m²)", color = "AA_only") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")
  print(p)
  dev.off()
} else {
  jpeg(plot_file, width = 1800, height = 1400, res = 220, quality = 95)
  cols <- as.integer(df$AA_only)
  plot(df$AgeCalc, df$BMI, col = cols, pch = 16,
       xlab = "Age (years)", ylab = "BMI (kg/m²)", main = "Age vs BMI")
  abline(lm(BMI ~ AgeCalc, data = df), lwd = 2, lty = 2)
  legend("topright", legend = levels(df$AA_only), col = seq_along(levels(df$AA_only)),
         pch = 16, bty = "n", title = "AA_only")
  mtext(subtitle_txt)
  dev.off()
}

message("Saved: ", plot_file)

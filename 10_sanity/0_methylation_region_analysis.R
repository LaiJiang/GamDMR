
# -------------------------------
# Sliding Window & Clustering Analysis on Methylation Region
# -------------------------------

library(dplyr)
library(ggplot2)
library(mgcv)
library(data.table)
library(tidyr)
library(stats)

# -------------------------------
# Load methylation long-format data
# -------------------------------
load("C:/Per/LaiJiang/Project/UQAC/meth/results/10_sanity/meth_dat_1.RData")  # meth_long assumed loaded

# -------------------------------
# Method 1: Sliding Window Analysis
# -------------------------------
window_size <- 50
step_size <- 10

unique_cpgs <- sort(unique(meth_long$start))
n_cpgs <- length(unique_cpgs)

sliding_results <- data.frame()

for (i in seq(1, n_cpgs - window_size + 1, by = step_size)) {
  window_cpgs <- unique_cpgs[i:(i + window_size - 1)]
  subdata <- meth_long %>% filter(start %in% window_cpgs)

  if (length(unique(subdata$AA_only)) < 2) next  # Skip if only one group

  mod <- tryCatch(
    gam(meth ~ s(start) + AA_only, data = subdata),
    error = function(e) NULL
  )

  if (!is.null(mod)) {
    pval <- summary(mod)$p.table["AA_only", "Pr(>|t|)"]
    sliding_results <- rbind(sliding_results, data.frame(
      window_start = min(window_cpgs),
      window_end = max(window_cpgs),
      pvalue = pval,
      n_cpgs = length(window_cpgs)
    ))
  }
}

# -------------------------------
# Method 2: Clustering Based on Methylation Pattern
# -------------------------------
meth_matrix <- meth_long %>%
  pivot_wider(names_from = sample, values_from = meth) %>%
  column_to_rownames("start") %>%
  as.matrix()

# Remove rows with NA
meth_matrix <- meth_matrix[complete.cases(meth_matrix), ]

# Perform hierarchical clustering
dist_matrix <- dist(meth_matrix)
hc <- hclust(dist_matrix, method = "ward.D2")
cluster_assignments <- cutree(hc, h = 0.3 * max(hc$height))

# Attach cluster labels
cluster_df <- data.frame(start = as.integer(rownames(meth_matrix)), cluster = cluster_assignments)
meth_clustered <- meth_long %>% inner_join(cluster_df, by = "start")

# Analyze each cluster with GAM
cluster_results <- data.frame()

for (cl in unique(meth_clustered$cluster)) {
  subdata <- meth_clustered %>% filter(cluster == cl)
  if (length(unique(subdata$AA_only)) < 2) next

  mod <- tryCatch(
    gam(meth ~ s(start) + AA_only, data = subdata),
    error = function(e) NULL
  )

  if (!is.null(mod)) {
    pval <- summary(mod)$p.table["AA_only", "Pr(>|t|)"]
    cluster_results <- rbind(cluster_results, data.frame(
      cluster = cl,
      start = min(subdata$start),
      end = max(subdata$start),
      pvalue = pval,
      n_cpgs = length(unique(subdata$start))
    ))
  }
}

# Save results
write.csv(sliding_results, "C:/Per/LaiJiang/Project/UQAC/meth/results/10_sanity/sliding_window_results.csv", row.names = FALSE)
write.csv(cluster_results, "C:/Per/LaiJiang/Project/UQAC/meth/results/10_sanity/cluster_results.csv", row.names = FALSE)

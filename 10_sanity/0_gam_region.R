#1. Predict methylation for both AA groups across region


# Create prediction grid
grid_pred <- expand.grid(
  start = seq(min(i_meth_long_cov$start), max(i_meth_long_cov$start), length.out = 1000),
  AA_only = c(0, 1)
)

vars_to_fix <- c("AgeCalc", "Sex", "Non.smoker", "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc", "sv1","sv2","sv3","sv4","sv5")

# Add fixed covariate means
for (v in vars_to_fix) {
  val <- mean(as.numeric(i_meth_long_cov[[v]]), na.rm = TRUE)
  grid_pred[[v]] <- rep(val, nrow(grid_pred))
}

grid_pred$Sex <- as.integer(mean(i_meth_long_cov$Sex, na.rm = TRUE))
grid_pred$Non.smoker <- as.integer(mean(i_meth_long_cov$Non.smoker, na.rm = TRUE))
grid_pred$FID <-   # random effect held at most freqnuent FID in i_meth_long_cov
  factor(names(which.max(table(i_meth_long_cov$FID)))) # most frequent FID


# Predict excluding random effect
grid_pred$fit <- predict(gam_model, newdata = grid_pred, type = "response", exclude = "s(FID)")


#Compute the difference between AA=1 and AA=0 curves



library(tidyr)
diff_df <- grid_pred %>%
  select(start, AA_only, fit) %>%
  pivot_wider(names_from = AA_only, values_from = fit, names_prefix = "AA_") %>%
  mutate(diff = AA_1 - AA_0)

# Identify peaks in differences
# Optional: Smooth difference to reduce noise
diff_df$diff_smooth <- stats::filter(diff_df$diff, rep(1/5, 5), sides = 2)

# Flag sub-regions with high differences
threshold <- quantile(abs(diff_df$diff_smooth), 0.95, na.rm = TRUE)  # top 5%
candidate_regions <- diff_df %>%
  mutate(high_diff = abs(diff_smooth) > threshold) %>%
  filter(high_diff == TRUE)


library(ggplot2)
ggplot_gam_region <- ggplot(diff_df, aes(x = start, y = diff_smooth)) +
  geom_line(color = "darkred") +
  geom_hline(yintercept = c(-threshold, threshold), linetype = "dashed", color = "blue") +
  labs(title = "Smoothed Methylation Difference (AA=1 - AA=0)",
       y = "Smoothed Δ arcsin(sqrt(Methylation))")

#saveplot

#ggsave(filename = paste0(PATH_wk, "results/10_sanity/2_col_gam_region_diff.jpg"), plot=ggplot_gam_region, width = 10, height = 5)



library(ggplot2)
library(mgcv)
library(dplyr)

# Create a new data grid across the observed start positions for both AA groups
new_data <- expand.grid(
  start = seq(min(i_meth_long_cov$start), max(i_meth_long_cov$start), length.out = 200),
  AA_only = c(0, 1)
)

# Fill in mean values for covariates (keep them constant across prediction)
covariate_means <- i_meth_long_cov %>%
  summarise(across(c(AgeCalc, Sex, Non.smoker, EOSINOpc, LYMPHOpc, MONOpc, NEUTROpc, sv1:sv5), mean, na.rm = TRUE))

new_data <- cbind(new_data, covariate_means[rep(1, nrow(new_data)), ])
new_data$FID <-   # random effect held at most freqnuent FID in i_meth_long_cov
  factor(names(which.max(table(i_meth_long_cov$FID)))) # most frequent FID

new_data$Sex <- 1
new_data$Non.smoker<- as.integer(floor(new_data$Non.smoker))
# Predict fitted methylation values
new_data$fit <- predict(gam_model, newdata = new_data, type = "response", exclude = "s(FID)")

new_data$fit[1:5]
# Plot
ggplot_gam <- ggplot(new_data, aes(x = start, y = fit, color = factor(AA_only))) +
  geom_line(size = 1) +
  labs(
    title = "Fitted GAM Methylation Curves by AA Status",
    x = "Genomic Position (start)",
    y = "Fitted arcsin(sqrt(Methylation))",
    color = "AA Status"
  ) +
  scale_color_manual(values = c("blue", "red"), labels = c("AA = 0", "AA = 1")) +
  theme_minimal()
# Save the plot
#ggsave(filename = paste0(PATH_wk, "results/10_sanity/2_col_gam_fitted_curves.jpg"), plot=ggplot_gam,width = 10, height = 5)  


############Possibly identify sub-regions where the smooth difference is sharp, for hypothesis-driven regional testing.


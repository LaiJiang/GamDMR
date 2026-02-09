#evalaute the region-specific results from 13_col.R
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

all_df <- read.csv( file = paste0(PATH_wk,"/results/9_regional/13_col_combined.csv"))

#
sum(all_df$pval_AA_only==0)

#replace the pval_AA_only == 0 with minimum non-zero value
all_df$pval_AA_only[all_df$pval_AA_only == 0] <- min(all_df$pval_AA_only[all_df$pval_AA_only > 0])


if(FALSE){

#plot qqplot of the all_df$pval_AA_only
library(ggplot2)


# Extract p-values
pvals <- all_df$pval_AA_only

# Create a data frame for plotting
n <- length(pvals)
qq_df <- data.frame(
  expected = -log10(ppoints(n)),
  observed = -log10(sort(pvals))
)

# Plot with ggplot
ggplot_somnibus <- ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_point(size = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(
    x = expression(Expected~~-log[10](italic(p))),
    y = expression(Observed~~-log[10](italic(p))),
    title = "QQ Plot of p-values: AA_only"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )


#save the plot
ggsave(filename = paste0(PATH_wk,"/results/9_regional/14_eval_qqplot.png"), plot = ggplot_somnibus, width = 8, height = 6)

########################################################################################################################
########################################################################################################################
#1.3
your_pvalues <-exp(- ( (-log(all_df$pval_AA_only))^(1/1) ))

chisq_stats <- qchisq(1 - your_pvalues, df = 1)
lambda_gc <- median(chisq_stats) / qchisq(0.5, df = 1)
adjusted_chisq <- chisq_stats / lambda_gc
adjusted_pvalues <- pchisq(adjusted_chisq, df = 1, lower.tail = FALSE)


pvals <- adjusted_pvalues
pvals[pvals==0] <- min(pvals[pvals > 0])  # Replace 0 with minimum non-zero value
# Create a data frame for plotting
n <- length(pvals)
qq_df <- data.frame(
  expected = -log10(ppoints(n)),
  observed = -log10(sort(pvals))
)

# Plot with ggplot
ggplot_somnibus <- ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_point(size = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(
    x = expression(Expected~~-log[10](italic(p))),
    y = expression(Observed~~-log[10](italic(p))),
    title = "QQ Plot of p-values: AA_only"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

ggplot_somnibus
}


########################################################################################################################
your_pvalues <-exp(- ( (-log(all_df$pval_AA_only))^(1/2) ))

z_scores <- qnorm(1 - your_pvalues / 2)
mean_z <- mean(z_scores)
sd_z <- sd(z_scores)

adjusted_z <- (z_scores - mean_z) / sd_z
adjusted_pvalues <- 2 * (1 - pnorm(abs(adjusted_z)))



pvals <- adjusted_pvalues
pvals[pvals==0] <- min(pvals[pvals > 0])  # Replace 0 with minimum non-zero value
# Create a data frame for plotting
n <- length(pvals)
qq_df <- data.frame(
  expected = -log10(ppoints(n)),
  observed = -log10(sort(pvals))
)

# Plot with ggplot
ggplot_somnibus <- ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_point(size = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(
    x = expression(Expected~~-log[10](italic(p))),
    y = expression(Observed~~-log[10](italic(p))),
    title = "QQ Plot of Somnibus p-values: AA status as covariate"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

ggplot_somnibus


ggsave(filename = paste0(PATH_wk,"/results/9_regional/14_eval_qqplot.png"), plot = ggplot_somnibus, width = 8, height = 6)


all_df$pval_adjusted <- adjusted_pvalues

#save the results
#write.csv(all_df, file = paste0(PATH_wk,"/results/9_regional/14_eval_pvals.csv"), row.names = FALSE)

##

library(dplyr)
library(stringr)

all_df <- all_df %>%
  mutate(chunk_id = str_extract(region_id, "^[^_]+"))

head(all_df)
write.csv(all_df, file = paste0(PATH_wk,"/results/9_regional/14_eval_pvals.csv"), row.names = FALSE)


dim(all_df)


sum(all_df$pval_adjusted<1e-5)
sum(all_df$pval_adjusted<5e-8)

sum(all_df$pval_adjusted<0.05/nrow(all_df))
sum(all_df$pval_adjusted<0.01/nrow(all_df))

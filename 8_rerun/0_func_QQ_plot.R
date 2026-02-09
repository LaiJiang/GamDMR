# Base‐R version
# ----------------
# 1. Extract and clean p-values



# install.packages("ggplot2")   # if you haven’t already
library(ggplot2)



#save this plot as a jpeg
jpeg(filename = paste0(PATH_wk,"/scr/8_rerun/results/8_QQplot_",plot_name,".jpeg"))



pv <- na.omit(manhattan_dat$pval)

# 2. Sort
pv_sorted <- sort(pv)

# 3. Compute expected quantiles under Uniform(0,1)
n <- length(pv_sorted)
exp_p <- (1:n) / (n + 1)

# 4. Plot observed vs. expected on –log10 scale
plot(
  -log10(exp_p), -log10(pv_sorted),
  pch = 19, cex = 0.4,
  xlab = "Expected -log10(p)",
  ylab = "Observed -log10(p)",
  main = paste0("QQ-plot : ", plot_name)
)
abline(0, 1, col = "red")

dev.off()
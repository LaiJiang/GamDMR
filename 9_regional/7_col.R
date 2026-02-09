#this file collect and evluatae the results of the somnibus analysis
#plot the smoothed covaraite effects from a sample result


library(SOMNiBUS)


PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
PATH_scr9 <- paste0(PATH_wk,"scr/9_regional/")
PATH_results <- paste0(PATH_scr9,"results/")

task_id <- 20
load(paste0(PATH_results,"/AA_",sprintf("%04d",task_id),"_somnbibus.RData"), verbose=TRUE)



(outs$region_6239325_6240956$reg.out)

summary(outs$region_6239325_6240956)

################################################################################################
# Extract result
res <- outs$region_6239325_6240956

# Positions of CpGs
pos <- res$uni.pos


# Number of covariates
ncovs <- res$ncovs

# Smoothed beta(t) effects for AA_only
# Note: est is a matrix with rows = CpG positions, columns = 1+ncovs
# Column 1 is intercept, Column 2 is for AA_only
res$reg.out.gam

beta_AA <- res$Beta.out[,2]

se_AA <- res$SE.out[,2]

jpeg(paste0(PATH_results, "smoothed_effect_example.jpg"), width = 800, height = 600)
# Plot with 95% confidence band
plot(pos, beta_AA, type = "l", lwd = 2, col = "blue",
     ylab = "Smoothed effect of AA_only", xlab = "Genomic position",
     main = "Smoothed effect of AA_only on methylation")

# Add confidence interval
lines(pos, beta_AA + 1.96 * se_AA, col = "blue", lty = 2)
lines(pos, beta_AA - 1.96 * se_AA, col = "blue", lty = 2)
abline(h = 0, col = "gray50", lty = 3)

dev.off()

f<-outs$region_6239325_6240956

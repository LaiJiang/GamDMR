#we collect the results from the previous analysis in 3_2_run.R


PATH_UQAC <- "C:/Per/LaiJiang/Project/UQAC/"

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"



#load the results from 3_2_run.R
models_res <- read.table(file=paste0(PATH_wk,"dat/AA_impute_simple.txt"),header=FALSE)

colnames(models_res) <- c("CpG",
                     "model1_beta","model1_pvalue", "model1_FID_sd",
                        "model2_beta","model2_pvalue", "model2_FID_sd",
                     "model3_beta","model3_pvalue", "model3_FID_sd")

#first plot the pvalues from the model2 and model3

#save plot
jpeg(paste0(PATH_wk, "/results/5_col_Fig1.jpeg"), width = 800, height = 800)
plot(-log10(models_res$model2_pvalue), -log10(models_res$model3_pvalue),
col="blue",xlab="-log10(p-value): model 2",ylab="model 3",main="-log10(p-values) from Model 2 versus Model 3", pch=19)
abline(a=0,b=1)
dev.off()



sum(models_res$model1_beta!=0)

sum(models_res$model2_pvalue<5e-8)

sum(models_res$model3_beta!=0)


# Load necessary package
if (!require(VennDiagram)) {
  install.packages("VennDiagram")
  library(VennDiagram)
}

# Assume models_res is a matrix with dimensions 6 x nCpGs,
# where:
# - Row 1: beta estimates for Model 1 (significant if != 0)
# - Row 3: p-values for Model 2 (significant if < 0.05)
# - Row 5: beta estimates for Model 3 (significant if != 0)

# Extract significant CpG indices for each model:
sig_model1 <- which(models_res$model1_beta!= 0)
sig_model2 <- which(models_res$model2_pvalue < 0.01)
sig_model3 <- which(models_res$model3_beta != 0)

# Optionally, print the number of significant CpGs per model:
cat("Model 1 significant CpGs: ", length(sig_model1), "\n")
cat("Model 2 significant CpGs: ", length(sig_model2), "\n")
cat("Model 3 significant CpGs: ", length(sig_model3), "\n")

# Create a list of the significant sets
sig_list <- list(
  Model1 = sig_model1,
  Model2 = sig_model2,
  Model3 = sig_model3
)

# Plot the Venn diagram.
venn.plot <- venn.diagram(
  x = sig_list,
  filename = NULL,  # if NULL, the diagram is returned as a grid object
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  main = "Overlap of Significant CpG Associations"
)


# Open a JPEG device to save the plot
jpeg(paste0(PATH_wk, "/results/5_col_Fig2.jpeg"), width = 800, height = 800)

# Draw the Venn diagram on the device
grid.draw(venn.plot)

# Close the device to finalize the file
dev.off()



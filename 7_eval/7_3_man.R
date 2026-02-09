#mahattan plot


library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(qqman)

library(VennDiagram)

#first load the results 

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
library(data.table)
AA_simple_results <- fread(paste0(PATH_wk,"/scr/6_beluga/all_simple.txt"), header=FALSE)
colnames(AA_simple_results) <- c("Chunk_index","CpG",
                                  "M1_coef","M1_pval","M1_F_sdv",
                                  "M2_coef","M2_pval","M2_F_sdv",
                                  "M3_coef","M3_pval","M3_F_sdv" ,
                                  "M1_lambda","M3_lambda")                                
                                  

#remove dupliated rows
AA_simple_results <- AA_simple_results[!duplicated(AA_simple_results), ]

#replace the NA in M2_coef as 0 
AA_simple_results$M2_coef[is.na(AA_simple_results$M2_coef)] <- 0
AA_simple_results$M2_pval[is.na(AA_simple_results$M2_pval)] <- 1

#plot M1_pval only for those M1_coef!=0
M1_pval <- AA_simple_results$M1_pval[AA_simple_results$M1_coef!=0]
M3_pval <- AA_simple_results$M3_pval[AA_simple_results$M3_coef!=0]
M2_pval <- AA_simple_results$M2_pval


sum(M1_pval < 0.05/length(M1_pval))
sum(M3_pval < 0.05/length(M3_pval))
sum(M2_pval < 0.05/length(M2_pval))

sum(AA_simple_results$M1_coef !=0)
###########################################################################
#first draw manhattan plot for M1_pval

M1_manhattan <- data.frame(CpG = AA_simple_results$CpG[AA_simple_results$M1_coef!=0],
                          pval = AA_simple_results$M1_pval[AA_simple_results$M1_coef!=0])


M3_manhattan <- data.frame(CpG = AA_simple_results$CpG[AA_simple_results$M3_coef!=0],
                          pval = AA_simple_results$M3_pval[AA_simple_results$M3_coef!=0])



M2_manhattan <- data.frame(CpG = AA_simple_results$CpG[AA_simple_results$M2_pval<0.05],
                          pval = AA_simple_results$M2_pval[AA_simple_results$M2_pval<0.05])


M2_manhattan <- data.frame(CpG = AA_simple_results$CpG,
                          pval = AA_simple_results$M2_pval)

source(paste0(PATH_wk,"scr/7_eval/func_plot_manhattan.R"))

man_plot1 <- man_plot(M1_manhattan, name ="M1_manhattan")
man_plot3 <- man_plot(M3_manhattan, name ="M3_manhattan")
man_plot2 <- man_plot(M2_manhattan, name ="M2_manhattan")

sum(M1_manhattan$pval < 1e-5)
sum(M1_manhattan$pval < 1e-4)

###########################################################################

sum(M1_manhattan$pval < 1e-5)
sum(M2_manhattan$pval < 1e-5)
sum(M3_manhattan$pval < 1e-5)


###########################################################################
# Assume models_res is a matrix with dimensions 6 x nCpGs,
# where:
# - Row 1: beta estimates for Model 1 (significant if != 0)
# - Row 3: p-values for Model 2 (significant if < 0.05)
# - Row 5: beta estimates for Model 3 (significant if != 0)

# Extract significant CpG indices for each model:
sig_model1 <- which(AA_simple_results$M1_pval<1e-5)
sig_model2 <- which(AA_simple_results$M2_pval<1e-5)

# Create a list of the significant sets
sig_list <- list(
  Model1 = sig_model1,
  Model2 = sig_model2)

# Plot the Venn diagram.
venn.plot <- venn.diagram(
  x = sig_list,
  filename = NULL,  # if NULL, the diagram is returned as a grid object
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  main = "Overlap of Significant CpGs from Model 1 and Model 2 (pvalue < 1e-5)"
)


# Open a JPEG device to save the plot
jpeg(paste0(PATH_wk, "/scr/7_eval/venn_12.jpeg"), width = 800, height = 800)
# Draw the Venn diagram on the device
grid.draw(venn.plot)
# Close the device to finalize the file
dev.off()
###########################################################################
# Assume models_res is a matrix with dimensions 6 x nCpGs,
# where:
# - Row 1: beta estimates for Model 1 (significant if != 0)
# - Row 3: p-values for Model 2 (significant if < 0.05)
# - Row 5: beta estimates for Model 3 (significant if != 0)

# Extract significant CpG indices for each model:
sig_model1 <- order(AA_simple_results$M1_pval, decreasing = FALSE)[1:1000]
sig_model2 <- order(AA_simple_results$M2_pval, decreasing = FALSE)[1:1000]
sig_model3 <- order(AA_simple_results$M3_pval, decreasing = FALSE)[1:1000]


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
  main = "Overlap of Top 1000 significant CpGs from three models"
)


# Open a JPEG device to save the plot
jpeg(paste0(PATH_wk, "/scr/7_eval/top1k_123.jpeg"), width = 800, height = 800)

# Draw the Venn diagram on the device
grid.draw(venn.plot)
# Close the device to finalize the file
dev.off()
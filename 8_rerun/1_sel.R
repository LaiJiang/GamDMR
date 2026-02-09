#this file select the CpGs that need to rerun
#and prepare the feature selection set for the rerun


library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)


#first load the results 

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"


#load the id_rerun to a file
id_rerun <-read.table( file = paste0(PATH_wk,"/scr/7_eval/7_1_id_rerun.txt"))

#load the complete resultswhich were extracted from the whole genome results on Beluga

library(data.table)
complete_results <- fread(paste0(PATH_wk,"/scr/7_eval/extracted_complete.txt"), header=FALSE)



# List of features to check
ft_names <- c("Intercept","meth_AA", "AgeCalc", "Sex", "Non.smoker", 
              "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc", 
              "sv1", "sv2", "sv3", "sv4", "sv5")

col_names_complete <- c("batch","CpG", 
paste0("model1_pvals_",ft_names), paste0("model1_coef_",ft_names) , "model1_FID_sdv",
paste0("model2_pvals_",ft_names), paste0("model2_coef_",ft_names) , "model2_FID_sdv",
paste0("model3_pvals_",ft_names), paste0("model3_coef_",ft_names) , "model3_FID_sdv"
)

#the original run did not include the pvalue of the intercept only from model1
#so we need to remove it from the col_names_complete
#remove feature "model1_vals_Intercept" from the col_names_complete
col_names_complete <- col_names_complete[-which(col_names_complete == "model1_pvals_Intercept")]

colnames(complete_results) <- col_names_complete


#sanity check
table(complete_results$model1_coef_meth_AA!=0, complete_results$model3_coef_meth_AA!=0)


#we actually do not need these feature selections, we just do it again in the rerun, would not make much differences 


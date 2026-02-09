#collect the BMI and age data and run association analysis
#check if the age and BMI are assocaited, and how they impact the AA status.


PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

library(data.table)
library(dplyr)

library(stringr)
library(mgcv)
library(qqman)

PATH_scr11 <- paste0(PATH_wk, "scr/11_mgcv/")

PATH_scr12 <- paste0(PATH_wk, "scr/12_cpg_DMR/")
PATH_results <- paste0(PATH_wk, "results/")

load(file = paste0(PATH_scr11, "dat/18_pheno_BMI.RData"), verbose = TRUE)


colnames(pheno_file)

#check the association between age and BMI
model <- lm(BMI ~ AgeCalc, data = pheno_file)
summary(model)


plot(pheno_file$AgeCalc, pheno_file$BMI,  
col= ifelse(pheno_file$AA_only == 1, "red", "blue"),
pch = 19,
     xlab = "Age", ylab = "BMI", main = "Scatter plot of Age vs BMI")

#generate the correlogram between all variables except ID and FID
library(GGally)
pheno_file_subset <- pheno_file %>% select(-ID, -FID)
# Create the correlogram
correlogram <- ggpairs(pheno_file_subset, aes(color = as.factor(AA_only), alpha = 0.5)) +
  theme_minimal() +
  labs(title = "Correlogram of Phenotypic Variables by AA Status", color = "AA Status")
# Display the correlogram
print(correlogram)

source(paste0(PATH_scr12, "0_corr.R"))

#this file rerun the univariate analysis with updates and comments from AM in Februrary




library(lmerTest)
library(glmmLasso)
library(lme4)
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

#load the updated data 

load(file=paste0(PATH_wk,"results/2_3_2_transformed.RData"),verbose=TRUE)



#SAVE RESULTS
load( file=paste0(PATH_wk,"results/2_3_3_svojb.RData"),verbose=TRUE)

#now attach sv variables to the phenotype data
sv <- svobj$sv 
colnames(sv) <- paste0("sv",1:ncol(sv))
pheno_file <- cbind(pheno_file, sv[,1:5])

#centralization and scaling of the covariates
pheno_file$AgeCalc <- scale(pheno_file$AgeCalc)
pheno_file$EOSINOpc <- scale(pheno_file$EOSINOpc)
pheno_file$LYMPHOpc <- scale(pheno_file$LYMPHOpc)
pheno_file$MONOpc <- scale(pheno_file$MONOpc)
pheno_file$NEUTROpc <- scale(pheno_file$NEUTROpc)
pheno_file$sv1 <- scale(pheno_file$sv1)
pheno_file$sv2 <- scale(pheno_file$sv2)
pheno_file$sv3 <- scale(pheno_file$sv3)
pheno_file$sv4 <- scale(pheno_file$sv4)
pheno_file$sv5 <- scale(pheno_file$sv5)
pheno_file$FID <- factor(pheno_file$FID)
#sex and non.smoker are two categorical variables, we need to convert them to factors


beta_asr_impute <- asin(sqrt(meth_matrix_imputed))
beta_asr <- asin(sqrt(meth_matrix))



##############################################################################
##In AM methylation prepration step, all CpGs with at least 120 individuals are retained. (out of 757) samples.
#this means the proportion of missing values for each CpG is capped at 0.84.
#Note that In my AA analysis, we filtered the data from 755 subject to 380 subjects.

#for each row in meth_matrix, how many values are missing?
#this is the number of missing values in each row of the matrix
missing_values <- rowSums(is.na(meth_matrix))/ncol(meth_matrix)

print(max(missing_values))

#we will be comparing the analysis result with original asr and imputed asr.
##############################################################################

save(
pheno_file, beta_asr, beta_asr_impute, meth_matrix, meth_matrix_imputed,
file=paste0(PATH_wk,"results/4_1_data.RData")  # Save the data to a file

)

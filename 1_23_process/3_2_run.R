

library(lmerTest)
library(glmmLasso)
library(lme4)
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

load( file=paste0(PATH_wk,"results/AA_only_final_model.RData"))

#first scale data 
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

#Approach 1: AA Status as the Outcome:

#model 1: binary logisitic regression

i_cpg <- 1

for(i_cpg in 1:nrow(meth_matrix_imputed)){

print(i_cpg)
meth_prop <- meth_matrix_imputed[i_cpg,]



# model 2: Fit a lasso-penalized mixed model, this would resolve the converge issue, but we would
# need to tune the penalty parameter lambda and also it lacks the interpretability of the coefficients
mod_glmmLasso <- glmmLasso(
  fix = AA_only ~ meth_prop + AgeCalc + Sex + Smoker + EOSINOpc + LYMPHOpc + MONOpc + sv1 + sv2 +
  sv3 + sv4 + sv5 ,
  rnd = list(FID = ~ 1),
  data = pheno_file,
  family = binomial(link = "logit"),
  lambda = 10  # adjust the penalty parameter as needed
)

#record the estimates of meth_prop
meth_prop_est_LASSO <- coef(mod_glmmLasso)["meth_prop"]
meth_FID_LASSO <- mod_glmmLasso$StdDev

#########################
#Approach 2: Methylation as the Outcome:


#we first tried the whole model, and met the singularity issue.



#then lets try simpler model
#model 2: simple one




lmer_simple <- lmer(meth_prop ~ AA_only + AgeCalc + Sex + Smoker + 
                   EOSINOpc + LYMPHOpc + MONOpc  + sv1 + sv2 +sv3 +sv4 + sv5 +(1|FID),  
            data = pheno_file)

f<-summary(lmer_simple)
# Summarize the model results
fixef_summary <- summary(lmer_simple)$coefficients
AA_est_lmer <- fixef_summary["AA_only", "Pr(>|t|)"]
#extract random effect estimates
FID_est_lmer <- ranef(lmer_simple)$FID[,1]

vc <- VarCorr(lmer_simple)

#we can also extract the standard deviation of the random effect
fid_std <- attr(vc$FID, "stddev")


#model 3: next, we try glmmLASSO to shrink fixed effect, while keep FID random effect



mod_glmmLasso_AA <- glmmLasso(
  fix = meth_prop  ~ AA_only + AgeCalc + Sex + Smoker + EOSINOpc + LYMPHOpc + MONOpc + sv1 + sv2 +
  sv3 + sv4 + sv5 ,
  rnd = list(FID = ~ 1),
  data = pheno_file,
  family =  gaussian(link = "identity"), 
  lambda = 10  # adjust the penalty parameter as needed
)

#record the estimates of meth_prop
AA_est_LASSO <- coef(mod_glmmLasso_AA)["AA_only"]

AA_FID_LASSO <- mod_glmmLasso_AA$StdDev

output_results <- c(meth_prop_est_LASSO, meth_FID_LASSO, 
             AA_est_lmer, fid_std,
              AA_est_LASSO, AA_FID_LASSO)

#attach this line to a file
write.table(output_results, file=paste0(PATH_wk,"results/AA_only_final_model_results.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")


}
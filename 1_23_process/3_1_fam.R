#collect family information and sva variables to generate the final dataset to be modeled.

impute_missing_meth <- TRUE

PATH_UQAC <- "C:/Per/LaiJiang/Project/UQAC/"

PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"

pheno_file <- read.csv(paste0(PATH_wk, "/results/AA_only.csv"))

#note that , the smoking status 
#("Non.smoker", "Ex.smoker", "Smoker"), 
#which are perfectly collinear (they sum to 1 for each individual).
# therefore we will only include the "Ex.smoker", "Smoker" into the model


#sanity check
#cor(pheno_file[,c("Non.smoker","Ex.smoker", "Smoker")])
#cor(pheno_file[,c("EOSINOpc","LYMPHOpc", "MONOpc", "NEUTROpc")])
#plot(apply(pheno_file[,c("EOSINOpc","LYMPHOpc", "MONOpc", "NEUTROpc")],1,sum))


# Full model: include your known covariates (replace 'Age', 'Sex', etc. with your actual covariates)
mod <- model.matrix(~ AgeCalc + Sex + Ex.smoker + Smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc, data = pheno_file)
# Null model: include only the intercept
mod0 <- model.matrix(~ 1, data = pheno_file)

#################
#now load the methylation profile filtered from whole set with var > 0.05

#meth_file <- read.csv(paste0(PATH_wk, "dat/peek_methylation_sva_0035.csv"),sep="\t")
meth_file <- read.csv(paste0(PATH_wk, "dat/peek_methylation_sva_short.csv"),sep="\t")

#1. match the IDs of AA and control 
# Extract only the methylation proportion columns (those that end with '_meth')
meth_data <- meth_file[, grep("_meth$", colnames(meth_file))]

# Convert the data frame to a matrix
meth_matrix <- as.matrix(meth_data)

#colect the IDs for meth_matrix
meth_ID <- colnames(meth_matrix)
meth_ID <- sub("^X", "", meth_ID)      # Remove the leading "X"
meth_ID <- sub("_meth$", "", meth_ID)   # Remove the trailing "_meth"
meth_ID <- gsub("\\.", "-", meth_ID)     # Replace the dot with a dash

# Check dimensions and a preview of the matrix

meth_matrix <- meth_matrix[,match(pheno_file$ID,meth_ID)]

if(impute_missing_meth==TRUE){
#2. address the missing value
# Impute missing values with the median for each CpG (row-wise)
meth_matrix_imputed <- t(apply(meth_matrix, 1, function(x) {
  med <- median(x, na.rm = TRUE)  # calculate median ignoring NAs
  x[is.na(x)] <- med              # replace NAs with the median value
  return(x)
}))
}

# Check for any remaining missing values
any(is.na(meth_matrix_imputed))



rm(meth_file)
rm(meth_data)
rm(meth_matrix)

#load sv variables
library(sva)
load(file=paste0(PATH_wk,"results/svojb.RData"),verbose=TRUE)

#we will use the top 5 sv variables to explain the methylation data
sv <- svobj$sv

#we now attach the family ID to the phenotype data
pheno_file$FID  <- factor( sub("-.*", "", pheno_file$ID))

#we now plot the histogram of family size to see the distribution of the family size
# Calculate family sizes
family_sizes <- as.numeric(table(pheno_file$FID))

jpeg(paste0(PATH_wk,"results/family_hist.jpg"))
hist(family_sizes,
     breaks = seq(0, max(family_sizes) + 1, by = 1),
     main = "Histogram of Family Sizes in the AA-only cohort (N=482)",
     xlab = "Family Size",
     ylab = "Frequency",
     col = "skyblue",
     border = "white")
dev.off()


###now attach the first 5 sv variables ontto the phenotype data

colnames(sv) <- paste0("sv-",1:ncol(sv))

pheno_file <- cbind(pheno_file, sv[,1:5])

#save the final dataset
write.csv(pheno_file, file = paste0(PATH_wk, "results/AA_only_final.csv"), row.names = FALSE)

#pheno_file <- read.csv( file = paste0(PATH_wk, "results/AA_only_final.csv"))

save(pheno_file, meth_matrix_imputed, file=paste0(PATH_wk,"results/AA_only_final_model.RData"))
#############################model



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

#binary logistic regression with FID as random effect


mod_logit <- glmer(AA_only ~ meth_prop + AgeCalc + Sex + Ex.smoker + Smoker + 
                   EOSINOpc + LYMPHOpc + MONOpc  + sv1 + sv2 + sv3 + sv4 + sv5 +
                   (1 | FID), 
                   data = pheno_file, family = binomial)

#conclusion:  the unidentifiablity issue can be addressed by scaling and removing 
#the highly correlated features. 
#but the model is still not converging. 


mod_logit_simple <- glmer(AA_only ~ meth_prop + AgeCalc + Sex +  Smoker + 
                   EOSINOpc + LYMPHOpc + MONOpc  + sv1 + sv2 +
                   (1 | FID), 
                   data = pheno_file, family = binomial)
#even with simple model. it did not converge.

#next, we try glmmLASSO to shrink fixed effect, while keep FID random effect 




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


# Fit the model
mod <- lmer(meth_prop ~ AA_only + AgeCalc + Sex + Ex.smoker + Smoker + 
                   EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5 +
              (1 | FID),
            data = pheno_file)

# Summarize the model results


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
             AA_est_lmer, FID_est_lmer, fid_std,
              AA_est_LASSO, AA_FID_LASSO)

              #write this line to a file
write.table(output_results, file = paste0(PATH_wk, "results/AA_only_final_model_results.txt"), row.names = FALSE, col.names = FALSE)

}


library(lmerTest)
library(glmmLasso)
library(lme4)
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"


load(file=paste0(PATH_wk,"results/4_1_data.RData"), verbose=TRUE  )

#Approach 1: AA Status as the Outcome:

#model 1: binary logisitic regression

i_cpg <- 1


print(i_cpg)

#the asr transformation
meth_response <- beta_asr_impute[i_cpg,]
#will try original beta_asr later

#####MOdel 1: AA status as outcome
#Fit a lasso-penalized mixed model, this would resolve the converge issue, but we would
# need to tune the penalty parameter lambda and also it lacks the interpretability of the coefficients


##two ways to obtain pvalues for glmmlasso:

#method 1: run glmmlasso first to select features, then run lmer on the selected set again to obtain pvalues.
#drawbacks: Note that this approach does not account for the uncertainty
# introduced during the variable selection process with glmmLASSO. It gives an approximate p‑value as if the model had been specified a priori.


#method 2: likelihood ratio test. run glmmlasso twice. 
#drabacks: A standard likelihood ratio test (LRT) relies on comparing the (unpenalized) likelihoods of two nested models. However, in glmmLASSO the likelihood is modified by a penalty term (the L1 norm on the fixed effects), so its value is not directly comparable to the “pure” likelihood. In other words, 
#the penalized likelihood does not have the same asymptotic distribution as a maximum‐likelihood estimator, and the usual chi‑square test does not apply.

##method 3: bootstrap approach. but it is too computational expensive.


if(TRUE){


# Define a grid of lambda values to test
lambda_grid <- c(5, 10, 20, 40, 50, 60, 80, 100)

# Initialize a vector to store BIC for each lambda
bic_values <- numeric(length(lambda_grid))
aic_values <- numeric(length(lambda_grid))

print("model 1: tune lambda")
for (i in seq_along(lambda_grid)) {
    print(i)
  mod_tmp <- glmmLasso(
    fix = AA_only ~ meth_response + AgeCalc + Sex + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5,
    rnd = list(FID = ~ 1),
    data = pheno_file,
    family = binomial(link = "logit"),
    lambda = lambda_grid[i]
  )
  # Store the BIC (or AIC) value; note that glmmLasso objects often have a 'BIC' component
  bic_values[i] <- mod_tmp$bic
  aic_values[i] <- mod_tmp$aic
}

#we choose BIC beacause it is more robust than AIC when the sample size is enough for model estimation
# Identify the optimal lambda (minimum BIC)
optimal_lambda <- lambda_grid[which.min(bic_values)]
cat("Model 1: Optimal lambda (minimizing BIC):", optimal_lambda, "\n")

#now rerurn with this optimal lambda
  mod1_optimal <- glmmLasso(
    fix = AA_only ~ meth_response + AgeCalc + Sex + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5,
    rnd = list(FID = ~ 1),
    data = pheno_file,
    family = binomial(link = "logit"),
    lambda = optimal_lambda
  )

  
#record the estimates of meth_prop
model1_coefficients <- coef(mod1_optimal)
# Extract the coefficients for the fixed effects
model1_coef_meth <- model1_coefficients["meth_response"]
# Extract the standard deviation of the random effect (FID)
model1_FID_sdv <- mod1_optimal$StdDev


#post selection inference: we can use the selected features to fit a lmer model to obtain pvalues


#step1:  Identify predictors with nonzero coefficients (ignoring the intercept)
selected_features <- names(model1_coefficients)[abs(model1_coefficients) > 0 & names(model1_coefficients) != "(Intercept)"]


# --- Step 2: Build the model formula for refitting ---
# Note: glmmLASSO was fit with fix = AA_only ~ meth_response + AgeCalc + Sex + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5.
# We will refit a model with AA_only as the outcome and the selected predictors as fixed effects,
# plus a random intercept for family (FID).
formula_str <- paste("AA_only ~", paste(selected_features, collapse = " + "), "+ (1|FID)")
#cat("Refit model formula:\n", formula_str, "\n")

# Convert the string to a formula object
formula_obj <- as.formula(formula_str)

# --- Step 3: Refit the model using glmer ---
mod_refit <- glmer(formula_obj, data = pheno_file, family = binomial(link = "logit"))
# --- Step 4: Extract p-values for the fixed effects ---
# The summary from glmer (with lmerTest loaded) will include p-values for fixed effects.
# List of features to check
features <- c("meth_response", "AgeCalc", "Sex", "Non.smoker", 
              "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc", 
              "sv1", "sv2", "sv3", "sv4", "sv5")

# Extract fixed effects summary from the refitted model (using glmer via lmerTest)
fixed_effects_table <- coef(summary(mod_refit))

# Initialize a named vector to store p-values
p_values <- setNames(rep(1, length(features)), features)

# Loop over each feature and update p_values if the feature is present in the fixed effects table
for(feat in features){
  if(feat %in% rownames(fixed_effects_table)){
    # Assuming the p-values are in the column "Pr(>|z|)"
    p_values[feat] <- fixed_effects_table[feat, "Pr(>|z|)"]
  } else {
    p_values[feat] <- 1
  }
}

#the pvalues of all fixed effects in the model
model1_coef_pvals <- p_values 
#model1_coefficients: fixed effects estimates
#model1_FID_sdv: random effect estimates
#there three variables need to be recorded for model1
cat("Model 1: finish ")


}
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################


#########################
#Approach 2: Methylation as the Outcome:


#we first tried the whole model, and met the singularity issue.



#then lets try simpler model
#model 2: simple one



if(TRUE){
lmer_simple <- lmer(meth_response ~ AA_only + AgeCalc + Sex + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5 +(1|FID),  
            data = pheno_file)

model2_summary<-summary(lmer_simple)
model2_coef_pvals <- model2_summary$coefficients[, "Pr(>|t|)"]

# Summarize the model results
fixef_summary <- model2_summary$coefficients
AA_est_lmer <- fixef_summary["AA_only", "Pr(>|t|)"]
#extract random effect estimates
FID_est_lmer <- ranef(lmer_simple)$FID[,1]
vc <- VarCorr(lmer_simple)
#we can also extract the standard deviation of the random effect
model2_FID_sdv <- attr(vc$FID, "stddev")
model2_coefficients <- fixef_summary[, "Estimate"]

cat("Model 2: finish ")

#three variables need to be recorded for model2
#model2_coefficients: fixed effects estimates
#model2_FID_sdv: random effect estimates
#model2_coef_pvals: pvalues of fixed effects
}


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################


#model 3: next, we try glmmLASSO to shrink fixed effect, while keep FID random effect

if(TRUE){


# Define a grid of lambda values to test
lambda_grid <- c(5, 10, 20, 40, 50, 60, 80, 100)

# Initialize a vector to store BIC for each lambda
bic_values <- numeric(length(lambda_grid))
aic_values <- numeric(length(lambda_grid))

print("model 3: tune lambda")
for (i in seq_along(lambda_grid)) {
    print(i)
  mod_tmp <- glmmLasso(
    fix = meth_response ~ AA_only + AgeCalc + Sex + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5,
    rnd = list(FID = ~ 1),
    data = pheno_file,
    family = gaussian(link = "identity"),
    lambda = lambda_grid[i]
  )
  # Store the BIC (or AIC) value; note that glmmLasso objects often have a 'BIC' component
  bic_values[i] <- mod_tmp$bic
  aic_values[i] <- mod_tmp$aic
}

#we choose BIC beacause it is more robust than AIC when the sample size is enough for model estimation
# Identify the optimal lambda (minimum BIC)
optimal_lambda <- lambda_grid[which.min(bic_values)]
cat("Model 3: Optimal lambda (minimizing BIC):", optimal_lambda, "\n")

#now rerurn with this optimal lambda
  mod3_optimal <- glmmLasso(
    fix = meth_response ~ AA_only + AgeCalc + Sex + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5,
    rnd = list(FID = ~ 1),
    data = pheno_file,
    family = gaussian(link = "identity"),
    lambda = optimal_lambda
  )

  
#record the estimates of meth_prop
model3_coefficients <- coef(mod3_optimal)
# Extract the coefficients for the fixed effects
model3_coef_meth <- model3_coefficients["meth_response"]
# Extract the standard deviation of the random effect (FID)
model3_FID_sdv <- mod3_optimal$StdDev

cat("Model 3: turn over")

#post selection inference: we can use the selected features to fit a lmer model to obtain pvalues


#step1:  Identify predictors with nonzero coefficients (ignoring the intercept)
selected_features <- names(model3_coefficients)[abs(model3_coefficients) > 0 & names(model3_coefficients) != "(Intercept)"]


# --- Step 2: Build the model formula for refitting ---
# Note: glmmLASSO was fit with fix = AA_only ~ meth_response + AgeCalc + Sex + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5.
# We will refit a model with AA_only as the outcome and the selected predictors as fixed effects,
# plus a random intercept for family (FID).
formula_str <- paste("meth_response ~", paste(selected_features, collapse = " + "), "+ (1|FID)")
#cat("Refit model formula:\n", formula_str, "\n")

# Convert the string to a formula object
formula_obj <- as.formula(formula_str)

# --- Step 3: Refit the model using glmer ---
mod_refit <- lmer(formula_obj, data = pheno_file )
# --- Step 4: Extract p-values for the fixed effects ---
# The summary from glmer (with lmerTest loaded) will include p-values for fixed effects.
# List of features to check
features <- c("AA_only", "AgeCalc", "Sex", "Non.smoker", 
              "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc", 
              "sv1", "sv2", "sv3", "sv4", "sv5")

# Extract fixed effects summary from the refitted model (using glmer via lmerTest)
fixed_effects_table <- coef(summary(mod_refit))

# Initialize a named vector to store p-values
p_values <- setNames(rep(1, length(features)), features)

# Loop over each feature and update p_values if the feature is present in the fixed effects table
for(feat in features){
  if(feat %in% rownames(fixed_effects_table)){
    # Assuming the p-values are in the column "Pr(>|z|)"
    p_values[feat] <- fixed_effects_table[feat, "Pr(>|t|)"]
  } else {
    p_values[feat] <- 1
  }
}

#the pvalues of all fixed effects in the model
model3_coef_pvals <- p_values 
#model3_coefficients: fixed effects estimates
#model3_FID_sdv: random effect estimates
#there three variables need to be recorded for model1

}

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#now record all results 

if(TRUE){

output_complete_results <- c(model1_coef_pvals, model1_coefficients , model1_FID_sdv, 
              model2_coef_pvals, model2_coefficients , model2_FID_sdv, 
              model3_coef_pvals, model3_coefficients , model3_FID_sdv) 


output_simple_results <- c(model1_coef_meth,
model1_coef_pvals["meth_response"],
model1_FID_sdv,
model2_coefficients["AA_only"],
model2_coef_pvals["AA_only"],
model2_FID_sdv,
model3_coefficients["AA_only"],
model3_coef_pvals["AA_only"],
model3_FID_sdv )

#attach this line to a file
write.table(output_complete_results, file=paste0(PATH_wk,"results/AA_impute_complete.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(output_simple_results, file=paste0(PATH_wk,"results/AA_impute_simple.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")


}

#next: 4_3_roaw.R, use raw value instead of imputed values


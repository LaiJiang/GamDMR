

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
meth_response <- beta_asr[i_cpg,]

pheno_file$meth_response <- meth_response

#remove NA rows

pheno_file2 <- na.omit(pheno_file)
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
print(mod_refit)
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

print(model1_coef_pvals)

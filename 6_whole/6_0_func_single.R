
# Define a grid of lambda values to test
lambda_grid <- c(5, 10, 20, 40, 50, 60, 80, 100)

# Initialize a vector to store BIC for each lambda
bic_values <- numeric(length(lambda_grid))
aic_values <- numeric(length(lambda_grid))

#print("model 1: tune lambda")
for (i in seq_along(lambda_grid)) {
    print(i)
  mod_tmp <- glmmLasso(
    fix = AA_only ~ meth_response + AgeCalc + Sex + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5,
    rnd = list(FID = ~ 1),
    data = pheno_file_feed,
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
optimal_lambda_model1 <- optimal_lambda
#cat("Model 1: Optimal lambda (minimizing BIC):", optimal_lambda, "\n")

#now rerurn with this optimal lambda and get pvalues 
  mod1_optimal <- glmmLasso(
    fix = AA_only ~ meth_response + AgeCalc + Sex + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5,
    rnd = list(FID = ~ 1),
    data = pheno_file_feed,
    family = binomial(link = "logit"),
    lambda = optimal_lambda,
    final.re = TRUE
  )

  
#record the estimates of meth_prop
model1_coefficients <- coef(mod1_optimal)
# Extract the coefficients for the fixed effects
model1_coef_meth <- model1_coefficients["meth_response"]
# Extract the standard deviation of the random effect (FID)
model1_FID_sdv <- mod1_optimal$StdDev



# List of features to check
features <- c("meth_response", "AgeCalc", "Sex", "Non.smoker", 
              "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc", 
              "sv1", "sv2", "sv3", "sv4", "sv5")

# Extract fixed effects summary from the refitted model (using glmer via lmerTest)
fixed_effects_table <- coef(summary(mod1_optimal))

# Initialize a named vector to store p-values
model1_coef_pvals <- fixed_effects_table[match(features,rownames(fixed_effects_table)), "p.value"]

# replace NA with 1
model1_coef_pvals[is.na(model1_coef_pvals)] <- 1
#note it is not necessary to test for the significance of the random effect (FID)

#report: model1_coefficients, model1_FID_sdv, model1_coef_pvals


######################################################
#model 2

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
            data = pheno_file_feed)

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

#cat("Model 2: finish ")

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


# Define a grid of lambda values to test
lambda_grid <- c(5, 10, 20, 40, 50, 60, 80, 100)

# Initialize a vector to store BIC for each lambda
bic_values <- numeric(length(lambda_grid))
aic_values <- numeric(length(lambda_grid))

#print("model 3: tune lambda")
for (i in seq_along(lambda_grid)) {
    #print(i)
  mod_tmp <- glmmLasso(
    fix = meth_response ~ AA_only + AgeCalc + Sex + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5,
    rnd = list(FID = ~ 1),
    data = pheno_file_feed,
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
optimal_lambda_model3 <- optimal_lambda

#cat("Model 3: Optimal lambda (minimizing BIC):", optimal_lambda, "\n")

#now rerurn with this optimal lambda
  mod3_optimal <- glmmLasso(
    fix = meth_response ~ AA_only + AgeCalc + Sex + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5,
    rnd = list(FID = ~ 1),
    data = pheno_file_feed,
    family = gaussian(link = "identity"),
    lambda = optimal_lambda,
    final.re = TRUE
  )

  
#record the estimates of meth_prop
model3_coefficients <- coef(mod3_optimal)
# Extract the coefficients for the fixed effects
model3_coef_meth <- model3_coefficients["meth_response"]
# Extract the standard deviation of the random effect (FID)
model3_FID_sdv <- mod3_optimal$StdDev



# Extract fixed effects summary from the refitted model (using glmer via lmerTest)
fixed_effects_table <- coef(summary(mod3_optimal))

# Initialize a named vector to store p-values
model3_coef_pvals <- fixed_effects_table[, "p.value"]
# replace NA with 1
model3_coef_pvals[is.na(model3_coef_pvals)] <- 1

###########################################################################

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#now record all results 

if(TRUE){

output_complete_results <- c(i_cpg, model1_coef_pvals, model1_coefficients , model1_FID_sdv, 
              model2_coef_pvals, model2_coefficients , model2_FID_sdv, 
              model3_coef_pvals, model3_coefficients , model3_FID_sdv) 


output_simple_results <- c(i_cpg, model1_coefficients['meth_response'],
model1_coef_pvals["meth_response"],
model1_FID_sdv,
model2_coefficients["AA_only"],
model2_coef_pvals["AA_only"],
model2_FID_sdv,
model3_coefficients["AA_only"],
model3_coef_pvals["AA_only"],
model3_FID_sdv, 
optimal_lambda_model1,
optimal_lambda_model3 )

}

#next: 4_3_roaw.R, use raw value instead of imputed values


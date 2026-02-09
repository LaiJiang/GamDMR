
yes_rerun_model1 <- rerun_chunk$V4[rerun_chunk$V2 == i_cpg]
yes_rerun_model3 <- rerun_chunk$V5[rerun_chunk$V2 == i_cpg]

#!!!!!!!!!!!!!!!!!!
#yes_rerun_model1 <- 1
#yes_rerun_model3 <- 1


#!!!
model1_coefficients <- rep(0,15)
model1_FID_sdv <- 0
model1_coef_pvals   <- rep(0,14)#note there is a difference here 

#!!!
model3_coefficients <- rep(0,15)
model3_FID_sdv <- 0
model3_coef_pvals   <- rep(0,15)#note there is a difference here 

features <- c("meth_response", "AgeCalc", "Sex", "Non.smoker", 
              "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc", 
              "sv1", "sv2", "sv3", "sv4", "sv5" , "BMI")

features_model3 <- c("AA_only", "AgeCalc", "Sex", "Non.smoker", 
              "EOSINOpc", "LYMPHOpc", "MONOpc", "NEUTROpc", 
              "sv1", "sv2", "sv3", "sv4", "sv5", "BMI")

#stopped here!!!!!!!!!!!!!!!
names(model1_coefficients) <- c("(Intercept)",features)
names(model3_coefficients) <- c("Intercept",features_model3)
names(model1_coef_pvals) <- features
names(model3_coef_pvals) <- c("Intercept",features_model3)

optimal_lambda_model3 <- 0
optimal_lambda_model1 <- 0

######

if(yes_rerun_model1){
# Define a grid of lambda values to test
lambda_grid <- c(5, 10, 20, 40, 50, 60, 80, 100)

# Initialize a vector to store BIC for each lambda
bic_values <- numeric(length(lambda_grid))
aic_values <- numeric(length(lambda_grid))

#print("model 1: tune lambda")
for (i in seq_along(lambda_grid)) {
    print(i)
  mod_tmp <- glmmLasso(
    fix = AA_only ~ meth_response + AgeCalc + Sex + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5 + BMI,
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
    fix = AA_only ~ meth_response + AgeCalc + Sex + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5 + BMI,
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


# Extract fixed effects summary from the refitted model (using glmer via lmerTest)
fixed_effects_table <- coef(summary(mod1_optimal))

# Initialize a named vector to store p-values
model1_coef_pvals <- fixed_effects_table[match(features,rownames(fixed_effects_table)), "p.value"]

# replace NA with 1
model1_coef_pvals[is.na(model1_coef_pvals)] <- 1
#note it is not necessary to test for the significance of the random effect (FID)

#report: model1_coefficients, model1_FID_sdv, model1_coef_pvals
######################################################
#!!!add the second stage refit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#now select these features are significant in model 1
model1_sel_fts <- names(model1_coefficients[model1_coefficients!=0])
#remove (Intercept) from the model1_sel_fts
model1_sel_fts <- model1_sel_fts[model1_sel_fts != "(Intercept)"]


if(length(model1_sel_fts)>0){
#and rerun model 2 with these features
# Define the formula for model 2
formula_model1_stage2 <- as.formula(paste0("AA_only ~", paste(model1_sel_fts, collapse = " + "), "+ (1|FID)"))

res1 <- tryCatch({


# 3) fit the non‑penalized GLMM
glmod <- glmer(
  formula  = formula_model1_stage2,
  data   = pheno_file_feed,
  family = binomial(link = "logit")
)

# 4) extract coefficients and p‑values
model1_FID_sdv <-   as.numeric(attr(VarCorr(glmod)$FID, "stddev"))

# 4a) extract p‑values via summary()
smry <- summary(glmod)$coefficients
# smry is a matrix with columns: Estimate, Std. Error, z value, Pr(>|z|)
stage2_pvals <- smry[, "Pr(>|z|)"]

#update model1_coef_pvals according to the stage2_pvals
model1_coef_pvals[model1_sel_fts] <- stage2_pvals[model1_sel_fts]
# and update model1_coefficients according to the stage2_pvals
model1_coefficients[model1_sel_fts] <- smry[model1_sel_fts, "Estimate"]


      list(
        coefficients = model1_coefficients,
        FID_sdv      = model1_FID_sdv,
        coef_pvals   = model1_coef_pvals
  )

}, error = function(e) {
  # on error, warn and return defaults
  warning("Model 1 refit: glmer() failed for CpG ", i_cpg, ": ", e$message)
  list(
    coefficients = rep(0, length(model1_coefficients)),
    FID_sdv      = 0,
    coef_pvals   = rep(1, length(model1_coef_pvals))
  )
})


# 3) Unpack into your variables
model1_coefficients <- res1$coefficients
model1_FID_sdv      <- res1$FID_sdv
model1_coef_pvals   <- res1$coef_pvals
}

}
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
if(yes_rerun_model3){

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
    fix = meth_response ~ AA_only + AgeCalc + Sex + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5 + BMI,
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
    fix = meth_response ~ AA_only + AgeCalc + Sex + Non.smoker + EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc + sv1 + sv2 + sv3 + sv4 + sv5 + BMI,
    rnd = list(FID = ~ 1),
    data = pheno_file_feed,
    family = gaussian(link = "identity"),
    lambda = optimal_lambda,
    final.re = TRUE
  )

  
#record the estimates of meth_prop
model3_coefficients <- coef(mod3_optimal)
# Extract the coefficients for the fixed effects !!!!!
model3_coef_meth <- model3_coefficients["AA_only"]
# Extract the standard deviation of the random effect (FID)
model3_FID_sdv <- mod3_optimal$StdDev



# Extract fixed effects summary from the refitted model (using glmer via lmerTest)
fixed_effects_table <- coef(summary(mod3_optimal))

# Initialize a named vector to store p-values
model3_coef_pvals <- fixed_effects_table[, "p.value"]
# replace NA with 1
model3_coef_pvals[is.na(model3_coef_pvals)] <- 1

###########################################################################
#now add the second stage refit also for model 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###############################################################################
# ====== Stage 2 refit for Model 3 (un-penalized LMM via lmer) ===============
###############################################################################

# 1) Identify which fixed effects survived the LASSO
model3_sel_fts <- names(model3_coefficients)[model3_coefficients != 0]
model3_sel_fts <- setdiff(model3_sel_fts, "(Intercept)")

if(length(model3_sel_fts)>0){
# 2) Build the stage-2 formula
#    meth_response is your continuous (logit-transformed) outcome
formula_model3_stage2 <- as.formula(
  paste0("meth_response ~ ", paste(model3_sel_fts, collapse = " + "), " + (1|FID)")
)

# 3) Refit in a tryCatch so you get fall-backs on error
res3 <- tryCatch({


  # refit without penalty
  lmod3 <- lmer(
    formula = formula_model3_stage2,
    data    = pheno_file_feed,
    REML    = FALSE    # or TRUE if you prefer REML
  )

  # extract the new random-effect SD for FID
  FID_sdv3 <- as.numeric(attr(VarCorr(lmod3)$FID, "stddev"))

  # pull the fixed effects table
  cf3 <- summary(lmod3)$coefficients
  # columns are: Estimate, Std. Error, df, t value, Pr(>|t|)
  pvals3    <- cf3[ , "Pr(>|t|)"]
  ests3     <- cf3[ , "Estimate"]

  # update only the selected features
  model3_coef_pvals[model3_sel_fts] <- pvals3[model3_sel_fts]
  model3_coefficients[model3_sel_fts] <- ests3[model3_sel_fts]

  list(
    coefficients = model3_coefficients,
    FID_sdv      = FID_sdv3,
    coef_pvals   = model3_coef_pvals
  )

}, error = function(e) {
  warning("Model 3 stage-2 refit failed for CpG ", i_cpg, ": ", e$message)
  list(
    coefficients = rep(0, length(model3_coefficients)),
    FID_sdv      = 0,
    coef_pvals   = rep(1, length(model3_coef_pvals))
  )
})

# 4) Unpack the results back into your variables
model3_coefficients <- res3$coefficients
model3_FID_sdv      <- res3$FID_sdv
model3_coef_pvals   <- res3$coef_pvals

}

}
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#now record all results 

if(TRUE){

output_complete_results <- c(i_cpg, meth_cpg_info, model1_coef_pvals, model1_coefficients , model1_FID_sdv, 
              #model2_coef_pvals, model2_coefficients , model2_FID_sdv, 
              model3_coef_pvals, model3_coefficients , model3_FID_sdv) 


output_simple_results <- c(i_cpg, meth_cpg_info, model1_coefficients['meth_response'],
model1_coef_pvals["meth_response"],
model1_FID_sdv,
#model2_coefficients["AA_only"],
#model2_coef_pvals["AA_only"],
#model2_FID_sdv,
model3_coefficients["AA_only"],
model3_coef_pvals["AA_only"],
model3_FID_sdv, 
optimal_lambda_model1,
optimal_lambda_model3 )

}

#next: 4_3_roaw.R, use raw value instead of imputed values


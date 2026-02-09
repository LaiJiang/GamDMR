library(data.table)
library(dplyr)

library(tidyr)   # <- needed for drop_na()

#now extract the genes
PATH_wk <- "C:/Per/LaiJiang/Project/UQAC/meth/"
outdir <- paste0(PATH_wk, "results/11_mgcv/")
PATH_scr6 <- paste0(PATH_wk, "scr/6_beluga/")
PATH_scr11 <- paste0(PATH_wk, "scr/11_mgcv/")


load(file = paste0(PATH_scr6, "4_1_data.RData"), verbose = TRUE)

#analyze the phenotypes, any of them associated with the AA?

#read a xlsx file with phenotypes
pheno_file <- paste0(PATH_scr11, "dat/PhenoTable complete_recrutement_2025_04_02.xlsx")
library(readxl)
df <- read_excel(pheno_file, sheet = "PhenoTable_all")


#extract column: "poids", "taille", "IMC" "Plaque"
#Amount of IgE in blood (µg/L) ?
#platelet number in 10^9/L?



#PH_OTH: Personal history of other disease, 2=the individual suffer or have already suffered from other disease, 1=no other disease
#not interesting.

#thought process 1:
#BMI already captures weight adjusted for height.
#Best practice: include BMI only, not all three, to reduce redundancy and instability in model fitting.
#Rationale: Obesity and BMI are well-documented confounders in asthma epidemiology and are associated with systemic methylation changes. Adjusting for BMI makes your AA–methylation associations more specific to asthma rather than obesity-related signals.

# process 2:
# IgE is a core biomarker of allergic sensitization and sits on the causal pathway between atopy and asthma.

#process 3: 
# Platelets are not a classical confounder of asthma but are relevant for systemic inflammation and immune response.
#Including platelets can help reduce false positives driven by general inflammation rather than asthma-specific processes.
#Downside: could reduce power if platelet variation is not strongly confounding your data.

#final strategy: 
#Core covariates (always include): Age, sex, smoking status, cell composition (lymphocytes, eosinophils, etc.), family ID.
#confoudner set: BMI. 
#platelets: only if the resulting model have systemaic inflammation effects. 

#plan: test for weight, height, BMI association with AA in phenoype file.
#run model by adding BMI. then see result if that makes sense. Then decdie if we need to add platelets.



#extract column: "poids", "taille", "IMC" "Plaque"
pheno_cols <- c("ID","poids", "taille", "IMC", "Plaque")
pheno_data <- df %>% 
  select(all_of(pheno_cols)) %>%
  drop_na() %>% rename(ID = ID, weight_kg = poids, height_m = taille, BMI = IMC, platelets = Plaque) %>%
    filter(!grepl("na", weight_kg, ignore.case = TRUE) & 
           !grepl("na", height_m, ignore.case = TRUE) & 
           !grepl("na", BMI, ignore.case = TRUE) & 
           !grepl("na", platelets, ignore.case = TRUE))


#merge pheno_file and pheno_data by ID
pheno_data_merged <- pheno_data %>%
  mutate(ID = as.character(ID)) %>%
  mutate(weight_kg = as.numeric(weight_kg),
         height_m = as.numeric(height_m),
         BMI = as.numeric(BMI),
         platelets = as.numeric(platelets)) %>%
  filter(!is.na(weight_kg) & !is.na(height_m) & !is.na(BMI) & !is.na(platelets)) %>%
    left_join(pheno_file %>% select(ID, AA_only), by = "ID") %>%
    na.omit()



#boxplot 

library(ggplot2)
library(dplyr)
library(tidyr)

# List of covariates you want to plot
covariates <- c("weight_kg", "height_m", "BMI", "platelets")

# Reshape data into long format for easier plotting
pheno_long <- pheno_data_merged %>%
  pivot_longer(cols = all_of(covariates),
               names_to = "Covariate",
               values_to = "Value")

# Boxplot for each covariate by AA_only
boxplots <- ggplot(pheno_long, aes(x = factor(AA_only), y = Value, fill = factor(AA_only))) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.7) +
  facet_wrap(~ Covariate, scales = "free_y") +
  labs(x = "Allergic Asthma (AA_only: 0 = Control, 1 = Case)",
       y = "Value",
       title = "Distribution of Covariates by AA Status") +
  scale_fill_manual(values = c("0" = "#4DAF4A", "1" = "#E41A1C"),
                    labels = c("0" = "Control", "1" = "Asthma")) +
  theme_bw(base_size = 13) +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

#test if the BMI is associated with AA_only
bmi_model <- lm(BMI ~ AA_only, data = pheno_data_merged)
summary(bmi_model)
# Save the boxplot
ggsave(filename = paste0(outdir, "18_pheno_boxplots.jpeg"), plot = boxplots, width = 10, height = 6)

#we found that On average, AA subjects had BMI ~1.15 units lower than controls in your cohort.
#This could be due to Allergic asthma (particularly in children/younger subjects) can occur in leaner individuals with strong atopic predisposition.
# The positive obesity–asthma relationship in the literature mostly refers to non-allergic, adult-onset asthma.


#Specifically, Asthma and obesity are often positively associated: obesity increases the risk of asthma,
#and many studies report higher BMI among asthma patients (especially in adults, women, and non-atopic asthma).
#However, in allergic asthma (AA), the relationship is more nuanced:
#Several pediatric and adolescent cohorts show that atopy-related asthma is not as strongly linked to obesity.
#In some studies of children and young adults, leaner individuals with high IgE/allergic sensitization have higher asthma prevalence compared to obese individuals.
#The “obese asthma” phenotype is usually non-allergic, late-onset asthma (different from your AA cohort).

#now curate data containg BMI

load(file = paste0(PATH_scr6, "4_1_data.RData"), verbose = TRUE)


#extract column: "poids", "taille", "IMC" "Plaque"
pheno_data <- df %>% 
  select("ID","IMC") %>%
  drop_na() %>% rename(BMI = IMC) %>%
    filter( 
           !grepl("na", BMI, ignore.case = TRUE) )


#merge pheno_file and pheno_data by ID
pheno_data_merged <- pheno_data %>%
  mutate(ID = as.character(ID)) %>%
  mutate(      BMI = as.numeric(BMI)) %>%
  filter( !is.na(BMI) ) %>%
    left_join(pheno_file, by = "ID") %>%
    na.omit() %>% #reorder columsn to move BMI to the end
  select(-BMI, everything(), BMI)

subject_selection <- match(pheno_data_merged$ID, pheno_file$ID)

pheno_file <- pheno_data_merged 
beta_asr <- beta_asr[, subject_selection]
beta_asr_impute <- beta_asr_impute[, subject_selection]
meth_matrix <- meth_matrix[, subject_selection]
meth_matrix_imputed <- meth_matrix_imputed[, subject_selection]

save(beta_asr, beta_asr_impute, meth_matrix, meth_matrix_imputed, pheno_file, file = paste0(PATH_scr11, "dat/18_pheno_BMI.RData"))
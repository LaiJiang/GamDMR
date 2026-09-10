PRIMARY SINGLE-CpG TRAIN/TEST PREDICTION
=========================================

Scripts
-------
16_run_primary_single_cpg_prediction_cv.R
16_run_primary_single_cpg_prediction_cv.sh
17_collect_primary_single_cpg_prediction_cv.R
17_collect_primary_single_cpg_prediction_cv.sh

Design
------
One SLURM array task = one of the same 100 outer family train/test splits.

For each split, three separate single-CpG prediction models are evaluated:

  SingleCpG_M1_LASSO
  SingleCpG_M2_LASSO
  SingleCpG_M3_LASSO

The split-specific selected file is:
  .../6_single_cpg/2_collected_significant/selected_by_split/
      single_cpg_selected_split_XXX.tsv

M1 features are rows with finite p_M1.
M2 features are rows with finite p_M2.
M3 features are rows with finite p_M3.

The union is used only to reduce raw methylation I/O. Each prediction model uses
only its own M1/M2/M3 feature set.

Critical zero-CpG behavior
--------------------------
If a method has zero selected CpGs in a split, that split is retained.

Example:
  M3 has zero CpGs in split 037

Then SingleCpG_M3_LASSO for split 037 uses the covariate-only fallback and records:
  model_fallback = TRUE
  fallback_reason = NO_SELECTED_CPGS
  n_CpG_selected_input = 0
  n_CpG_usable = 0

The same fallback is used if CpGs were selected but all fail TRAIN-only
missingness/variance QC.

This is essential for paired comparison across all 100 outer splits.

Prediction pipeline
-------------------
Individual selected CpG beta values are transformed:
  asin(sqrt(beta))

CpG QC is TRAIN only:
  observed fraction >= 0.50
  training median imputation
  same training median applied to test missing values
  training variance > 1e-8

Prediction:
  logistic LASSO
  alpha = 1
  lambda = lambda.1se
  type.measure = auc
  same shared family-grouped 5-fold inner CV used by regional methods
  CpG penalty.factor = 1
  AgeCalc / Sex / Non.smoker / BMI penalty.factor = 0
  Youden threshold from inner OOF training predictions only
  untouched outer test families used only for final performance

Metrics match the regional primary pipelines:
  AUROC, PR-AUC, PR-AUC lift, Brier, log loss,
  calibration intercept/slope, sensitivity, specificity,
  balanced accuracy, PPV, NPV.

Smoke test
----------
First test one outer split and cap at five CpGs per model:

  sbatch \
    --array=1 \
    --export=ALL,MAX_CPGS_PER_MODEL=5,OVERWRITE=1 \
    16_run_primary_single_cpg_prediction_cv.sh

Inspect:
  cat \
  ~/scratch/UQAC/meth/results/15_revision/5_cv/5_prediction/6_single_cpg_primary/performance/single_cpg_primary_performance_split_001.tsv

The performance file should contain four rows:
  SingleCpG_M1
  SingleCpG_M2
  SingleCpG_M3
  Covariates_only

Full run
--------
  sbatch 16_run_primary_single_cpg_prediction_cv.sh

Default:
  #SBATCH --array=1-100%50

Output
------
~/scratch/UQAC/meth/results/15_revision/5_cv/5_prediction/6_single_cpg_primary/

  performance/
  predictions/
  coefficients/
  feature_qc/
  split_qc/

Collect after all 100 jobs
--------------------------
  sbatch 17_collect_primary_single_cpg_prediction_cv.sh

Collected output:
  .../6_single_cpg_primary/collected/
    single_cpg_primary_performance_all_splits.tsv
    single_cpg_primary_predictions_all_splits.tsv
    single_cpg_primary_coefficients_all_splits.tsv
    single_cpg_primary_feature_qc_all_splits.tsv
    single_cpg_primary_split_qc_all_splits.tsv
    single_cpg_primary_performance_summary.tsv
    single_cpg_primary_collection_qc.tsv

Later all-method comparison
---------------------------
The performance schema deliberately matches the core regional-method schema:
  splitID
  DMR_method
  model
  model_type
  alpha
  lambda
  lambda_rule
  threshold
  threshold_source
  inner_folds
  fit_status
  model_fallback
  n_train / n_test
  *_inner_oof metrics
  *_test metrics

Single-CpG adds:
  n_CpG_selected_input
  n_CpG_usable

and preserves:
  n_DMR_selected_input = NA
  n_DMR_usable = NA

Thus the later comparison can combine with:
  GAM-DMR
  SOMNiBUS
  DMRcate
  BSmooth
  SingleCpG_M1
  SingleCpG_M2
  SingleCpG_M3

using rbindlist(..., fill=TRUE).

The covariate-only baseline should be retained only once per split in the final
all-method comparison.

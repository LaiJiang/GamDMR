ALL-METHOD PRIMARY PREDICTION COMPARISON
========================================

Scripts
-------
18_combine_all_primary_prediction_methods.R
18_combine_all_primary_prediction_methods.sh

Purpose
-------
Combine the PRIMARY 100-split train/test prediction results for:

  GAM-DMR
  SOMNiBUS
  DMRcate
  BSmooth
  Single-CpG M1
  Single-CpG M2
  Single-CpG M3
  Covariates only

The covariate-only model is taken ONCE from the GAM-DMR pipeline so duplicate
baseline rows from the other pipelines are not included.

Required collected inputs
-------------------------
5_prediction/1_mgcv_primary/collected/
5_prediction/2_somnibus_primary/collected/
5_prediction/3_dmrcate_primary/collected/
5_prediction/4_bsmooth_primary/collected/
5_prediction/6_single_cpg_primary/collected/

Every method must contain all 100 outer family-aware splits.

Run
---
mkdir -p logs
sbatch 18_combine_all_primary_prediction_methods.sh

Default output
--------------
~/scratch/UQAC/meth/results/15_revision/5_cv/5_prediction/
  7_all_method_comparison/

  data/
    all_methods_primary_performance_all_splits.tsv
    all_methods_primary_predictions_all_splits.tsv
    all_methods_subject_average_test_predictions.tsv
    all_methods_subject_average_ROC_coordinates.tsv
    all_methods_subject_average_PR_coordinates.tsv

  tables/
    Table1_all_methods_primary_performance.tsv
    Table1_all_methods_primary_performance_formatted.tsv
    Table2_GAM_DMR_vs_all_alternatives_corrected_CV.tsv
    Table_S1_all_methods_pairwise_corrected_repeated_CV.tsv
    Table_S2_all_methods_subject_average_performance.tsv

  figures/
    Figure1_all_methods_split_performance_distributions.pdf
    Figure1_all_methods_split_performance_distributions.png
    Figure4A_all_methods_subject_average_ROC.pdf
    Figure4A_all_methods_subject_average_ROC.png
    Figure4B_all_methods_subject_average_PR.pdf
    Figure4B_all_methods_subject_average_PR.png

ROC / PR curves
---------------
As in the previous regional comparison, repeated held-out test probabilities are
first averaged WITHIN SUBJECT and method across all outer splits in which the
subject appeared in the test set.

The resulting subject-averaged predictions are then used to construct secondary
descriptive ensemble ROC and precision-recall curves.

The ROC legend is displayed as, for example:
  GAM-DMR (AUROC=0.74)

The PR legend is displayed as:
  GAM-DMR (AUPRC=0.64)

IMPORTANT:
  AUROC and AUPRC values in the FIGURE LEGENDS are formatted to exactly
  TWO decimal places using sprintf("%.2f", ...).

Underlying TSV files retain full numerical precision.

The PR figure includes a dashed horizontal line at observed AA prevalence.

Corrected repeated-CV comparisons
---------------------------------
Pairwise comparisons use paired split-specific metric differences and the
corrected repeated-CV SE:

  sqrt((1/R + mean(n_test/n_train)) * var(d))

Effects are oriented so POSITIVE always favors method1:
  AUROC: method1 - method2
  PR-AUC: method1 - method2
  Brier: method2 - method1

Holm-adjusted p-values are computed separately within each metric family.

Table2 contains:
  GAM-DMR vs SOMNiBUS
  GAM-DMR vs DMRcate
  GAM-DMR vs BSmooth
  GAM-DMR vs Single-CpG M1
  GAM-DMR vs Single-CpG M2
  GAM-DMR vs Single-CpG M3
  GAM-DMR vs Covariates only

M3 zero-CpG splits
------------------
Single-CpG M3 splits that used the pre-specified covariate-only fallback remain
in the comparison. They are NOT dropped. Therefore the paired 100-split design
is preserved.

Interpretation of subject-average ROC/PR
----------------------------------------
The subject-average ROC/PR curves are SECONDARY descriptive ensemble summaries.
Primary inference remains the distribution of held-out outer-split metrics and
the corrected repeated-CV paired comparisons.

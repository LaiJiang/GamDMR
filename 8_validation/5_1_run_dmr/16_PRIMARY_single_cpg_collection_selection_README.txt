PRIMARY SINGLE-CpG COLLECTION AND FINAL PREDICTOR SELECTION
===========================================================

Stage 1
-------
Run only AFTER all 990 primary association jobs are complete:

  sbatch 14_collect_primary_single_cpg_significant.sh

This launches 99 collector jobs. Each collector reads about 10 raw association
files one at a time, so the full ~1.38-billion-row result set is never loaded
into memory at once.

It creates under:
  ~/scratch/UQAC/meth/results/15_revision/5_cv/6_single_cpg/2_collected_significant/partial/

  association_summary_collect_001.tsv
  association_candidates_collect_001.tsv
  association_collect_001_qc.tsv
  ...
  through collection task 099.

The collector counts p-values at:
  0.05, 0.01, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8

and retains rows with:
  p <= 0.001

The wider p<=0.001 retained tail allows stricter threshold sensitivity without
rereading all raw files.

Check Stage 1:
  find \
  ~/scratch/UQAC/meth/results/15_revision/5_cv/6_single_cpg/2_collected_significant/partial \
  -maxdepth 1 -type f -name "association_collect_*_qc.tsv" | wc -l

Expected: 99


Stage 2
-------
After all 99 collector jobs finish:

  sbatch 15_finalize_primary_single_cpg_selection.sh

Primary default selection rule:
  fit_status == FIT_OK
  AND p_value <= 1e-5
  within the CURRENT outer training split.

Then take the UNION of M1, M2, M3 within that split.

A CpG significant in multiple models is retained ONCE as a prediction feature.
Its split-specific file records which models selected it and the model-specific
p-values/coefficients.

To use a different pre-specified threshold <=0.001:
  sbatch \
    --export=ALL,SIGNIFICANCE_P=1e-6 \
    15_finalize_primary_single_cpg_selection.sh


Main outputs
------------
Root:
  .../6_single_cpg/2_collected_significant/

tables/association_significance_summary_by_split_model.tsv
  Counts for every split x M1/M2/M3 at all recorded p thresholds, plus the
  split/model Bonferroni threshold.

tables/association_significance_summary_overall_model.tsv
  Across-split descriptive summary.

tables/bonferroni_significant_associations.tsv
  Bonferroni-significant associations. This is descriptive unless you choose
  Bonferroni as the feature rule.

tables/primary_significant_associations.tsv
  Long-format associations passing SIGNIFICANCE_P.

tables/primary_significant_counts_by_split_model.tsv
  Significant counts for M1/M2/M3 in each split.

tables/primary_selected_cpg_union_all_splits.tsv
  All split-specific final CpG predictor sets in one table.

tables/primary_selected_cpg_counts_by_split.tsv
  Number of M1/M2/M3 hits and deduplicated final CpGs per split.

tables/primary_selected_cpg_model_overlap_by_split.tsv
  Counts of M1/M2/M3 overlap patterns.

tables/primary_selected_cpg_stability_DESCRIPTIVE_ONLY.tsv
  Across-split selection frequency. DO NOT use this table to choose predictors
  for the outer-CV prediction model.

Prediction input files
----------------------
selected_by_split/single_cpg_selected_split_001.tsv
...
selected_by_split/single_cpg_selected_split_100.tsv

These are the files to use in train/test prediction.

For outer split s, use ONLY:
  single_cpg_selected_split_s.tsv

Do not use across-split stability to select CpGs.

Recommended prediction feature construction
-------------------------------------------
For each outer split:
  1. read its selected CpGs;
  2. pull individual beta values from the raw methylation chunks;
  3. beta -> asin(sqrt(beta));
  4. TRAIN-only observed-fraction filter;
  5. TRAIN-only median imputation;
  6. TRAIN-only zero-variance filter;
  7. use individual CpGs as penalized methylation predictors;
  8. use AgeCalc, Sex, Non.smoker, BMI as unpenalized covariates;
  9. use the same shared family-grouped inner folds and logistic LASSO setup
     used for GAM-DMR/SOMNiBUS/DMRcate/BSmooth;
 10. evaluate on untouched outer-test families.

If a split has zero selected CpGs, retain the split and run the covariate-only
fallback.

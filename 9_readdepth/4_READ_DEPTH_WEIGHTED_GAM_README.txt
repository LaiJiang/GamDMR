READ-DEPTH-WEIGHTED GAM-DMR SENSITIVITY ANALYSIS
================================================

Purpose
-------
This workflow directly addresses the reviewer's suggestion to incorporate
read-depth weights into the GAM.

The weighted model is intentionally identical to the uploaded
20_run_region_BMI.R model except that each subject-CpG observation receives a
relative precision weight proportional to its total read depth:

  depth_weight = read_depth / mean(read_depth within region)

The model remains:

  meth_arcsin ~
    s(start, bs="cs", k=k_cpg) +
    s(start, by=AA_only, bs="cs", k=k_cpg) +
    AA_only + AgeCalc + Sex + Non.smoker +
    EOSINOpc + LYMPHOpc + MONOpc + NEUTROpc +
    sv1 + sv2 + sv3 + sv4 + sv5 + BMI +
    s(FID, bs="re")

with method="REML".

Why weight by read depth?
-------------------------
For a binomial methylation proportion p_hat=Y/n, the arcsine-square-root
transformation has approximate sampling variance proportional to 1/n. Thus,
weights proportional to n are a natural inverse-variance sensitivity analysis.

Weights are divided by the within-region mean depth so their mean is 1.
This keeps only the RELATIVE read-depth precision information and avoids an
arbitrary overall change in the weight scale.

The script additionally enforces read_depth >= 5.

Files
-----

1. 22_run_region_BMI_read_depth_weighted.R
   Runs weighted GAM-DMR region by region.

2. 22_run_region_BMI_read_depth_weighted.sh
   SLURM array runner. Default:
     --array=1-900%60
     N_GAM_JOBS=900

3. 23_collect_compare_read_depth_weighted_GAM.R
   Collects weighted results and compares them with the original unweighted
   files:
     results/11_mgcv/BMI_B1/7_results_job_*.txt

4. 23_collect_compare_read_depth_weighted_GAM.sh
   Runs the collector/comparison.

Run
---

mkdir -p logs_weighted_GAM

sbatch 22_run_region_BMI_read_depth_weighted.sh

After all array jobs finish:

sbatch 23_collect_compare_read_depth_weighted_GAM.sh

Weighted output
---------------

~/scratch/UQAC/meth/results/15_revision/7_read_depth_weighted_GAM/

Comparison output
-----------------

~/scratch/UQAC/meth/results/15_revision/7_read_depth_weighted_GAM/comparison/

Key outputs
-----------

comparison/weighted_GAM_sensitivity_results_summary.txt

comparison/tables/weighted_GAM_all_regions.tsv
comparison/tables/unweighted_GAM_all_regions_reconstructed.tsv
comparison/tables/weighted_vs_unweighted_region_comparison.tsv
comparison/tables/weighted_GAM_concordance_summary.tsv

comparison/figures/Figure_weighted_vs_unweighted_pvalue_concordance.pdf/.png
comparison/figures/Figure_weighted_vs_unweighted_effect_concordance.pdf/.png
comparison/figures/Figure_weighted_vs_unweighted_DMR_overlap.pdf/.png

Primary interpretation
----------------------
The strongest evidence of robustness would be:

- high correlation between weighted and unweighted regional association
  statistics;
- high correlation between weighted and unweighted effect-size estimates;
- substantial retention of the original GAM-DMR calls;
- substantial Jaccard overlap.

The collector uses the following default final-DMR definition:

  BH-FDR < 0.05
  EDF > 0.5
  mean_diff > 0.05

These are configurable in the collector shell script. If the final manuscript
uses different DMR thresholds, edit those environment variables before running
the collector.

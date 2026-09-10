SVA vfilter sensitivity analysis
================================

Purpose
-------
Address the reviewer request to assess whether selecting the top 10,000,
25,000, or 50,000 variable CpGs changes SVA estimation or latent sample
clustering.

The input file is the already-extracted high-variance methylation dataset:
  ~/scratch/UQAC/meth/dat/peek_methylation_sva_0035.csv

The user reported:
  dim = 85,936 x 1,517
so it is already large enough to support vfilter=50,000. No additional
genome-wide CpG extraction is needed.

Files
-----
0_prepare_sva_sensitivity_input.R
    Reads only analytic-subject *_meth columns from the 85,936-CpG file,
    applies ASR transformation, performs CpG-wise median imputation on the
    ASR scale, removes zero-variance rows, and saves one reusable input RDS.

0_prepare_sva_sensitivity_input.sh
    SLURM launcher for the preparation step.

1_run_sva_threshold_sensitivity.R
    Runs one vfilter setting: 10k, 25k, or 50k. The SVA model is kept
    identical to the original code:
      mod  = ~ AgeCalc + Sex + smoking + EOSINOpc + LYMPHOpc +
               MONOpc + NEUTROpc + BMI
      mod0 = ~ 1
      num.sv(..., method="be")
      sva(...)

    AA_only is deliberately NOT added because the goal here is to isolate
    the effect of CpG-selection threshold relative to the primary SVA
    specification.

1_run_sva_threshold_sensitivity.sh
    Three-task SLURM array:
      task 1 = 10,000
      task 2 = 25,000
      task 3 = 50,000

2_summarize_sva_threshold_sensitivity.R
    Compares number of SVs, absolute SV correlations, entire SV subspaces,
    pairwise subject geometry, and hierarchical clustering.

2_summarize_sva_threshold_sensitivity.sh
    SLURM launcher for summary/QC.

Recommended run order
---------------------
1. In the directory containing these scripts:
     mkdir -p logs

2. Prepare the 85,936 x analytic-subject SVA matrix:
     sbatch 0_prepare_sva_sensitivity_input.sh

   Check the log carefully. It should report 349 analytic subjects.

3. After preparation finishes:
     sbatch 1_run_sva_threshold_sensitivity.sh

4. Wait for all 3 array tasks to finish.

5. Summarize:
     sbatch 2_summarize_sva_threshold_sensitivity.sh

Default output directory
------------------------
~/scratch/UQAC/meth/results/15_revision/7_sva_threshold_sensitivity/

Key outputs
-----------
sva_sensitivity_input.rds
sva_input_qc.csv

sva_vfilter_10000.rds
sva_vfilter_25000.rds
sva_vfilter_50000.rds

SV_vfilter_10000.csv
SV_vfilter_25000.csv
SV_vfilter_50000.csv

SVA_nSV_summary.csv
SVA_pairwise_stability.csv
SVA_primary_SV_match_summary.csv

abs_cor_SV_10k_vs_25k.csv
abs_cor_SV_10k_vs_50k.csv
abs_cor_SV_25k_vs_50k.csv

heatmap_abs_cor_SV_*.png
subject_distance_*.png
SVA_sensitivity_summary.txt

Important notes
---------------
1. The reviewer specifically requested 10k, 25k, and 50k, so these are the
   only thresholds tested.

2. Individual SV signs/order can change without meaning that the latent
   structure changed. Therefore the summary includes rotation-invariant
   subspace and subject-distance metrics.

3. The prepared matrix uses ASR-transformed methylation values, matching the
   manuscript description of the primary SVA analysis.

4. sva requires complete methylation matrices. This workflow uses CpG-wise
   median imputation on the ASR scale. If the original primary SVA used a
   different imputation rule, change the preparation script to that exact
   rule before using the sensitivity results in the manuscript.

5. This workflow stops at SVA stability. Once the results are inspected, the
   25k and 50k SVs can be substituted into the existing GAM-DMR pipeline for
   a small downstream DMR-stability sensitivity analysis.

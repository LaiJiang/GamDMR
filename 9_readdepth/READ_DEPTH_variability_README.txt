READ-DEPTH VARIABILITY ANALYSIS
===============================

Purpose
-------
This workflow quantifies read-depth heterogeneity across:
1. genomic CpGs / genomic position;
2. analytic subjects;
3. chromosomes;
4. allergic-asthma groups descriptively.

It uses the same wide methylation chunk files as the GAM-DMR/SOMNiBUS workflow:

  ~/scratch/UQAC/meth/data/meth_split/chunk_XXXX.csv

Expected sample columns:
  <sample>_meth  methylation proportion
  <sample>_tot   total read depth

The *_tot definition is the same one used by the existing
convert_to_somnibus_format() helper.

Primary cohort
--------------
By default, the workflow restricts the data to subjects present in:

  ~/scratch/UQAC/meth/scr/11_mgcv/dat/18_pheno_BMI.RData

object:
  pheno_file

This corresponds to the analytic cohort used by the GAM-DMR revision analyses.

Two read-depth definitions
--------------------------
all_depth:
  Every finite, non-negative *_tot value.

gam_observed:
  Read depth only when total depth is >=5 and the paired *_meth value is finite
  and within [0,1]. This explicitly mirrors the stated low-depth filter and most
  closely reflects subject-CpG observations that can contribute to the GAM-DMR
  methylation-proportion analysis.

Step 1: run chunk-wise summaries
--------------------------------

mkdir -p logs

sbatch 20_run_read_depth_variability_chunks.sh

Default:
  300 array tasks
  maximum 60 concurrent jobs

Each task automatically discovers the chunk files and receives a round-robin
subset. Therefore, the script does not assume exactly 1737 chunks.

Smoke test
----------

sbatch \
  --array=1 \
  --export=ALL,N_READ_DEPTH_JOBS=300,OVERWRITE=1 \
  20_run_read_depth_variability_chunks.sh

Step 2: collect
---------------

After all array jobs finish:

sbatch 21_collect_read_depth_variability.sh

Output root
-----------

~/scratch/UQAC/meth/results/15_revision/6_read_depth_variability/

Important tables
----------------
tables/read_depth_subject_summary.tsv
  Exact subject-level mean/SD/CV, call rate, threshold fractions, and
  histogram-derived median/IQR read depth.

tables/read_depth_subject_variability_overall.tsv
  Distribution of subject-level depth metrics across the analytic cohort.

tables/read_depth_subject_by_chromosome.tsv
  Subject x chromosome depth summaries.

tables/read_depth_cpg_variability_overall.tsv
  Genome-wide distribution of CpG-specific mean/median/SD/CV/call rate.

tables/read_depth_genomic_bins_1Mb.tsv
  1-Mb genomic-bin summaries for plotting spatial coverage heterogeneity.

tables/read_depth_by_chromosome.tsv
  Chromosome-level summaries.

tables/read_depth_observation_histogram.tsv
  Observation-level read-depth distribution.

tables/read_depth_summary_by_AA_status_DESCRIPTIVE.tsv
  Descriptive depth summaries by AA status.

tables/read_depth_cpg_summary_random100k.tsv
  Reproducible random 100,000-CpG sample for manual inspection.

Figures
-------
figures/Figure_read_depth_A_across_subjects.pdf/.png
figures/Figure_read_depth_B_CpG_median_distribution.pdf/.png
figures/Figure_read_depth_C_CpG_CV_distribution.pdf/.png
figures/Figure_read_depth_D_across_genome_1Mb.pdf/.png
figures/Figure_read_depth_E_by_AA_status_DESCRIPTIVE.pdf/.png

Manuscript-ready output
-----------------------
read_depth_results_summary.txt

This file dynamically reports:
- observation-level mean/median/IQR/90th/95th percentile depth;
- fractions with depth 5-9, 10-19, 20-29, 30-49, >=50;
- between-subject variability in mean/median depth and CV;
- between-CpG variability in median depth, CV, and call rate.

Interpretation
--------------
These diagnostics directly characterize read-depth heterogeneity, but they do
not replace a read-depth-weighted GAM or count-based (binomial/beta-binomial)
sensitivity analysis. They help determine how large and structured the
measurement-precision heterogeneity actually is and provide manuscript/rebuttal
statistics for the reviewer's read-depth concern.

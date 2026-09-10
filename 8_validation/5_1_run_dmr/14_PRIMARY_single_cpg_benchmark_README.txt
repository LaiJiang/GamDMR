PRIMARY SINGLE-CpG BENCHMARK ASSOCIATION ANALYSIS
=================================================

Purpose
-------
This version removes the full-data/outcome-based CpG preselection used for the pilot.

The CpG universe is built directly from the methylation chunk files using only:
  chr, start, end, data_chunk_id, cpg_row_in_chunk

No full-data p-values, coefficients, AA labels, final_table rows, or significant-CpG
lists are used to construct the primary universe.

Every CpG is then analyzed separately within each outer TRAINING-family split.

Models
------
M1:
  AA_only ~ ASR methylation + covariates + (1|FID)
  first-stage glmmLasso, lambda from minimum BIC over:
    5, 10, 20, 40, 50, 60, 80, 100
  then unpenalized glmer refit of selected fixed effects.

M2:
  ASR methylation ~ AA_only + covariates + (1|FID)
  direct lmerTest::lmer, REML=FALSE.

M3:
  ASR methylation ~ AA_only + covariates + (1|FID)
  first-stage Gaussian glmmLasso, lambda from minimum BIC over the same grid,
  then unpenalized lmer refit of selected fixed effects.

Unselected target terms in M1/M3 are returned as coef=0, p=1, matching pilot v3.

Files
-----
12_prepare_primary_single_cpg_manifest.R
12_prepare_primary_single_cpg_manifest.sh
13_run_primary_single_cpg_training_associations.R
13_run_primary_single_cpg_training_associations.sh

Step 1: build outcome-independent manifest
------------------------------------------
mkdir -p logs

sbatch \
  --export=ALL,OVERWRITE_MANIFEST=1 \
  12_prepare_primary_single_cpg_manifest.sh

After it completes:

cat \
~/scratch/UQAC/meth/results/15_revision/5_cv/6_single_cpg/0_primary_all_cpg_manifest/primary_manifest_qc.tsv

The manuscript reference count is 4,609,564 QC-passing CpGs. The preparation script
reports any difference but does not apply outcome-based filtering.

Step 2: smoke test the PRIMARY runner
-------------------------------------
This tests 2 raw CpGs x 2 splits x M1/M2/M3 = 12 output rows and OVERWRITES job 0001
in the pilot output directory:

sbatch \
  --array=1 \
  --export=ALL,MAX_CPGS=2,MAX_SPLITS=2,OVERWRITE=1 \
  13_run_primary_single_cpg_training_associations.sh

Check:

head -20 \
~/scratch/UQAC/meth/results/15_revision/5_cv/6_single_cpg/1_training_associations/job_results/single_cpg_training_assoc_job_0001.tsv

cat \
~/scratch/UQAC/meth/results/15_revision/5_cv/6_single_cpg/1_training_associations/job_qc/single_cpg_training_assoc_job_0001_qc.tsv

Step 3: full run
----------------
The shell is configured as:

  #SBATCH --array=1-990%100
  MODELS=M1,M2,M3
  OVERWRITE=1

Submit with:

sbatch 13_run_primary_single_cpg_training_associations.sh

IMPORTANT COMPUTE WARNING
-------------------------
This exact primary analysis is far larger than the pilot.

If the primary universe has ~4.61 million CpGs:
  4.61M CpGs x 100 splits x 3 models
  = ~1.383 BILLION CpG-model-split association analyses.

With 990 jobs, each job has roughly:
  ~4,656 CpGs x 100 splits x 3 models
  = ~1.397 million model analyses/job.

Therefore the 24-hour walltime in the supplied shell is only a placeholder for a
representative timing test; it should NOT be assumed sufficient.

Before launching all 990 jobs, run a representative timing benchmark, for example:

sbatch \
  --array=1 \
  --export=ALL,MAX_CPGS=20,MAX_SPLITS=100,OVERWRITE=1 \
  13_run_primary_single_cpg_training_associations.sh

Then inspect elapsed time and MaxRSS:

sacct -j <JOBID> \
  --format=JobID,Elapsed,MaxRSS,ReqMem,State,ExitCode

Output compatibility
--------------------
Primary outputs intentionally use the SAME result directory, filenames, and result
columns as the pilot:

.../6_single_cpg/1_training_associations/job_results/
  single_cpg_training_assoc_job_XXXX.tsv

.../6_single_cpg/1_training_associations/job_qc/
  single_cpg_training_assoc_job_XXXX_qc.tsv

The columns:
  original_full_data_pval
  original_full_data_coef
  genes_entrez
  genes_symbol
are deliberately NA in the primary analysis, because importing them would reintroduce
information from the full-data prefiltered analysis.

The compatibility column final_table_row_id is now a deterministic logical ID:
  CpG1:M1 = 1
  CpG1:M2 = 2
  CpG1:M3 = 3
  CpG2:M1 = 4
  ...
It is NOT derived from the old final_table.

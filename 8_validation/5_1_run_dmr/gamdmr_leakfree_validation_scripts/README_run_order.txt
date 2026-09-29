Leakage-free GAM-DMR validation scripts
=======================================

Files
-----
00_make_outer_family_splits.R
01_run_mgcv_leakfree.R
02_collect_select_mgcv_leakfree.R
run_mgcv_leakfree.sh

Recommended run order on Rorqual
--------------------------------

1. Create outer family-aware splits:

   Rscript 00_make_outer_family_splits.R

2. Create a logs directory:

   mkdir -p logs

3. Smoke-test split 1:

   sbatch --export=ALL,SPLIT_ID=1 run_mgcv_leakfree.sh

4. After all split-1 jobs finish, collect/select the results.

   The current collector loops over all 100 splits and will warn for missing ones:

   Rscript 02_collect_select_mgcv_leakfree.R

5. If split 1 looks correct, submit all 100 splits:

   for SPLIT in $(seq 1 100)
   do
       sbatch --export=ALL,SPLIT_ID=${SPLIT} run_mgcv_leakfree.sh
   done

6. After every job is finished:

   Rscript 02_collect_select_mgcv_leakfree.R

Important
---------
The validation discovery pipeline deliberately does NOT load:

  region_GAM_M12.csv
  full-cohort GAM DMRs
  M1_manhattan
  M2_manhattan
  24_dmrs_STRICT.tsv

Each outer split performs GAM fitting and DMR selection using only that
split's training families.

Threshold note
--------------
02_collect_select_mgcv_leakfree.R currently uses R2 >= 0.50, matching the
V10 manuscript text. Your Supplementary Table S1 currently reports R2 >= 0.30.
Resolve this discrepancy before the final analysis.

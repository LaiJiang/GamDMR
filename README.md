# Repository Overview

This repository contains a modular methylation-analysis workflow for single-CpG association testing, regional DMR discovery, sensitivity analyses, benchmarking, and downstream reporting. The project covers several complementary methods: univariate models, SOMNiBUS, GAM-DMR, BSmooth, and DMRcate, together with additional validation and robustness analyses.

Only core scripts required to reproduce the workflow are included. Cohort-specific data and temporary output files are intentionally omitted for privacy reasons. Users must provide their own methylation datasets and harmonize input paths before running the pipeline.

---

# Structure

## 1. Data preparation and preprocessing (1_process)
- `1_process/`
  - `1_pheno_clean.ipynb`: prepares phenotype and covariate tables from the study cohort.
  - `2_var_clean.py`: collects and filters CpG-level variables for exploratory analyses and downstream model inputs.

---

## 2. Single-CpG association analysis (2_univariate)
- `2_univariate/`
  - `0_func_single_combined_clean.R`: core function for single-CpG association testing.
  - `1_job_clean.R`: driver script for running the univariate analysis in batch or cluster jobs.
  - `1_job_clean.sh`: HPC submission script for the single-CpG workflow.

---

## 3. Somnibus regional analysis (3_somnibus)
- `3_somnibus/`
  - `0_func_region_clean.R`: main regional association function for the SOMNiBUS model.
  - `0_func_somnibus_input_clean.R`: prepares the model input objects required for regional testing.
  - `1_test_chunk_clean.R`: runs a single-region chunk/test job for parallel execution.

---

## 4. GAM-DMR regional analysis (4_mgcv)
- `4_mgcv/`
  - `0_split_spacing_clean.R`: defines and organizes CpG regions for GAM-DMR analysis.
  - `0_load_results_updated_clean.R`: standardizes and loads batch outputs from GAM-DMR jobs.
  - `1_run_region_clean.R`: executes the GAM-DMR model for one region at a time.
  - `1_run_region_clean.sh`: batch submission script for parallel regional runs on a cluster.

---

## 5. BSmooth regional analysis (5_BSmooth)
- `5_BSmooth/`
  - `bsmooth_region_clean.R`: runs the BSmooth-based regional methylation analysis and prepares region-level summaries.

---

## 6. DMRcate regional analysis (6_DMRcate)
- `6_DMRcate/`
  - `dmrcate_region_clean.R`: runs the DMRcate-based regional differential methylation analysis.

---

## 7. Post-processing, overlap, and manuscript reporting (7_paper)
- `7_paper/`
  - `3_overlap_clean.R`: evaluates overlap between CpG-level signals and DMR-level findings for the M1 analysis.
  - `4_overlap_M2_clean.R`: extends the overlap analysis to the M2 model and compares concordance across levels.
  - `6_col_BSmooth_clean.R`: aggregates BSmooth results across jobs and retains regions with significant DMR calls.
  - `7_col_DMRcate_clean.R`: combines DMRcate outputs across jobs and filters for significant regions.
  - `8_overlap_4methods_clean.R`: summarizes and visualizes DMR overlap across multiple methods.
  - `9_supplement_clean.R`: generates supplementary tables and annotated gene-level summaries for the manuscript.
  - `12_main_figures_clean.R`: produces the main figures used in the paper, including examples and pathway-related visualizations.

---

## 8. Sensitivity analyses for latent-variable specification (10_sva_sensitivity)
- `10_sva_sensitivity/`
  - `0_prepare_sva_sensitivity_input.R`: prepares the SVA input matrix under alternative CpG-selection thresholds.
  - `0_prepare_sva_sensitivity_input.sh`: job launcher for the input-preparation step.
  - `1_run_sva_threshold_sensitivity.R`: runs SVA estimation across alternative vfilter settings (e.g., 10k, 25k, 50k).
  - `1_run_sva_threshold_sensitivity.sh`: array job submission script for the threshold sensitivity analysis.
  - `2_summarize_sva_threshold_sensitivity.R`: compares SV stability, clustering, and subspace consistency across settings.
  - `2_summarize_sva_threshold_sensitivity.sh`: launcher for the summary/QC step.
  - `3_run_sva_fixed5.R` and `4_summarize_sva_fixed5.R`: additional fixed-5 SVA workflow for targeted sensitivity checks.
  - `README(1).txt`: usage notes and interpretation guidance for the SVA sensitivity workflow.

This folder addresses reviewer-driven robustness checks for whether CpG selection thresholds materially change the latent-variable structure used in downstream methylation modeling.

---

## 9. Epigenomic overlap and methylation QTL-related analyses (11_mqtl)
- `11_mqtl/`
  - `0_prepapre_DMR.R`: prepares DMR inputs or derives region-level features used in downstream overlap analysis.
  - `1_col_filnames.sh`: collects file names for batch processing.
  - `2_epigen_dmr_overlap.py`: performs overlap analysis between epigenomic features and DMR calls.
  - `3_run_epigen.sh`: launches the overlap pipeline on HPC infrastructure.
  - `4_summarize_epigen.py`: collates and summarizes the overlap results across regions or features.

This folder focuses on integrating DMR findings with external epigenomic annotations and related molecular overlap summaries.

---

## 10. DMR report and summary generation (12_dmr_report)
- `12_dmr_report/`
  - `1_prepare_DMR.R`: prepares the DMR feature table and metadata for downstream reporting.
  - `2_prepare_DMR_chunk_manifest.R`: creates the chunk manifest used to distribute DMR extraction across array jobs.
  - `2_submit_DMR_pipeline.sh`: submits the DMR-characteristics pipeline in parallel.
  - `3_extract_DMR_statistics_array.R`: extracts statistics for each chunk in array form.
  - `3_extract_DMR_statistics_array.sh`: shell launcher for the chunked statistics extraction.
  - `4_summarize_DMR_statistics.R`: combines chunk-level outputs into final DMR summary tables.
  - `4_summarize_DMR_statistics.sh`: launcher for the final summarization step.
  - `README_V2.txt`: notes on the chunk-bundled DMR reporting workflow and SLURM job design.

This folder produces the finalized DMR summaries and characteristic tables used for interpretation and reporting.

---

## 11. Validation and method benchmarking (8_validation)
- `8_validation/`
  - `5_1_run_dmr/`: scripts for running DMR validation or benchmark jobs across training splits.
  - `5_2_collect/`: collection and aggregation scripts for combining split-based DMR results.
  - `5_3_build/`: scripts for final validation build, performance summary, and comparison outputs.

The validation folder contains a larger benchmarking workflow for comparing regional methods under alternative training splits, collection steps, and final prediction/validation summaries.

Examples of the scripts in this directory include:
- `1_test.R`, `2_run.sh`, `2_run_mgcv_train_splits.R`
- `3_run_reversal.R`, `4_somnibus_prep.R`, `5_prepare_somnibus_bundles.R`, `6_combine_somnibus_manifests.R`
- `7_run_somnibus_training_splits.R`, `8_prep_DMRcate.R`, `9_DMRcate_family_training_splits.R`
- `10_prep_BSmooth.R`, `12_prepare_cpg.R`, and the primary prediction/CV scripts (`1_run_mgcv_primary_prediction_cv.R`, `3_run_somnibus_primary_prediction_cv.R`, etc.)

This component is used to benchmark DMR methods in prediction and validation settings rather than only in the primary discovery workflow.

---

## 12. Read-depth sensitivity analysis (9_readdepth)
- `9_readdepth/`
  - `1_run_read_depth_variability_chunks.R` and `.sh`: evaluate read-depth variability across methylation chunks.
  - `2_collect_read_depth_variability.R` and `.sh`: aggregate the chunk-level variability summaries.
  - `3_run_region_BMI_read_depth_weighted.R` and `.sh`: fit the GAM-DMR model with read-depth-based weights.
  - `4_collect_compare_read_depth_weighted_GAM.R` and `.sh`: compare weighted vs unweighted GAM-DMR results and summarize concordance.
  - `5_collect_compare_read_depth_GAM_gene.R`: compares gene-level findings from the read-depth weighted sensitivity analysis.
  - `READ_DEPTH_variability_README.txt` and `4_READ_DEPTH_WEIGHTED_GAM_README.txt`: detailed notes on the analysis and interpretation.

This folder tests whether read-depth weighting changes the main GAM-DMR conclusions, which is important for assessing robustness to measurement precision.

---

# Notes
- SLSJ cohort data, raw methylation matrices, and temporary pipeline outputs are not included in the repository for privacy and size reasons; users need to use their own data, or obtain data permission from Catherine Lab and SLSJ cohort authorities. 
- The scripts are intentionally modular so they can be adapted to similar methylation datasets and cluster environments.
- The repository combines discovery, sensitivity, benchmarking, and reporting steps into one reproducible analysis framework.
- For questions about specific scripts or run order, please contact the repository maintainer Lai Jiang or leave comments..
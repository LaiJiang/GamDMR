# Repository Overview

This repository contains scripts for methylation association analyses, including single-CpG analysis, regional (DMR) analysis with GAM-DMR, and other regional based methods such as DMRcate and BSmooth. The scripts are organized into sections corresponding to different stages of the analysis workflow.

Only core scripts required to reproduce the main results are included. SLSJ cohort data input and temporary data files are not included for privacy reasons. Users must obtain and prepare their own methylation datasets as input.

---

# Structure

## 1. Data Preparation (1_process)
- `1_process/`  

  1_pheno.ipynb:    Prepare phenotype and covariates from SLSJ cohort data.
  2_var.py:    Example script to collect sample CpGs for preliminary data exploration and analysis.

---


## 2. Whole-Genome Analysis with Univaraite Models (2_univariate)
- `2_univariate/`  
  0_func_single.R: the main function for single-CpG association analysis.
  1_job.R: the script to run the analysis in parallel on HPC clusters.
  1_job.sh: the batch script to submit jobs on HPC clusters.

---

## 3. Somnibus (DMR) Analysis
- `3_somnibus/`  
  Regional association analysis with somnibus. 
  0_func_region.R: the main function for regional association analysis with Somnibus.
  0_func_somnibus_input.R: the function to prepare input for somnibus.
  1_test_chunk.R: the script to run the Somnibus analysis on a single region in parallel on HPC clusters.

---

## 4. GAM-DMR Analysis

- `4_mgcv/`  
  Regional association analysis with GAM-DMR method as in our manuscript.

 0_split_spacing_clean.R — Utinity function to define CpG Regions
 0_load_results_updated.R - Utinity function to standardize and load model outputs from batch jobs.
 1_run_region_clean.R - Run the core GAM-DMR model for a single region.
 1_run_region_batch_clean.sh - Batch script to run GAM-DMR across multiple regions in parallel on HPC clusters.

---

## 5-6. Other methods
- `5_BSmooth/`  
   Regional association analysis with BSmooth method.

- `6_DMRcate/`  
  Regional association analysis with DMRcate method.
---


## 7. Post analysis and visualization

- `3_overlap_clean.R`
Computes overlap between significant CpG sites from M1 model and DMR regions, and annotates region-level results.

- `4_overlap_M2_clean.R`
Extends overlap analysis for the M2 model and compares concordance between CpG-level and region-level signals.

- `6_col_BSmooth_clean.R`
Aggregates BSmooth per-region results across jobs and extracts regions with at least one detected DMR.

- `7_col_DMRcate_clean.R`
Aggregates DMRcate per-region results across jobs and filters regions with significant DMR findings.

- `8_overlap_4methods_clean.R`
Quantifies and visualizes overlap of DMRs across multiple methods (SOMNiBUS, GAM-DMR, BSmooth, DMRcate).

- `9_supplement_clean.R`
Functional scripts to generates supplementary tables (as in manuscript) of significant CpGs across models and annotates gene-level overlaps.

- `12_main_figures_clean.R`
Functional scripts to generate main manuscript figures (as in manuscript) including pathway enrichment plots and example regional methylation visualizations.

---

# Notes
- SLSJ cohort data and temporary data files are not included for privacy reasons; users must obtain and prepare their own methylation datasets as input.
- All scripts are designed to be modular to allow users to adapt them to their own methylation datasets.
- For questions on scripts, please contact the repository maintainer LJ.
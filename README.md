# Scripts Dictionary

## 123_process

The data processing scripts are organized in the following way:


### 1_pheno.ipynb

curate phenotype AD and AA, as well asd covariates, for Aim 1.

jupyter nbconvert --to html --no-input 1_pheno.ipynb


### 2_sva.Rmd

2_1_sva.r:
is the experimental scripts to try sva.

calculate surrogate variables to add to phenotype. use sva.

Rscript -e "rmarkdown::render('2_4_sva.Rmd', output_format = 'html_document')"

2_2_var.py:
 is the experiment code to filter methylation probes with high variance. 
 we finally run the script on the whole dataset on Beluga in folder /home/laj773/scratch/UQAC/meth/scr/: sva_meth_job.sh with sva_meth.py

2_3_1_dat.r: collect the pheno and methylation data to create the matched data for analysis.

2_3_2_transform.r: perform sensitivity analysis to determine which transformation to use .

2_3_3_SVA.r: run sva on the data after transformation. 

archive/2_3_sva_whole.R:
the old file to  calculate the surrogate variables with the methylation data extracted from sva_meth.py.

2_4_sva.Rmd: generate html report based on pictures and results. 


### 3_4_models.Rmd

3_fam.R:
collect family information and sva variables to generate the final dataset to be modeled.

Rscript -e "rmarkdown::render('3_4_models.Rmd', output_format = 'html_document')"

run univaraite models with 3 models. generate report.


## 4_univariate_pilot.Rmd

update the analysis with data corrections and comments from the group. 

4_1_dat.R:

aggregate all updated data for the analysis. 

4_2_single.R:
Experimental: we try a single CpG analysis with the new data.

4_3_func_single.R: the function to run a single CpG analysis.

4_4_run.R: load data, determine impute/raw data to run and run 4_3_func_single.R.


## 5_collect analysis

we collect analysis from the results run on Beluga and feed that into the 3_4_models.Rmd report.


## 6_whole genome association analysis


6_whole: 
the develop repo for whole genome association analysis.

6_beluga:
the production repo scripts for whole genome association analysis.

on beluga: 
6_4_eval.sh: collect the number of CpGs procesed in each data chunk to see which chunk is not finished. 
6_5_file_uncomplete.txt: extracted the data chunk ID that is not finished.
6_5_rerun_Job.sh: reurn the job only for the uncomplete data chunk. 


we run the whole genome association analysis with the new data.

to do: 
1. we can collect the pvalues with the final.re=TRUE option on the optimal model. fix code: 30 minutes. how to collect qvalues. FDR adjusted pvalues?
2. litearture review of how other methylation analysis is done with logistic. check with GPT which model (1 or 3) is suitable for this research proposal. 
2 hours.
3. chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.jacionline.org/action/showPdf?pii=S0091-6749%2824%2902278-4. check this abstract. how their qvalue, 
pvalue, FDR are used.
2. stopped here: we fix the error and warning issues on previous runs.  3 hours.
3. need to write script to collect 1 cpg, and convert to data for estimation.
3. figure out how to run the analysis on beluga parllelly maybe for 500 jobs. 8 hours.
4. run the analysis on beluga. 4 days.

1737 files to analyze.

5 minutes for 100 rows.
each file have 3000 rows.

#SBATCH --array=1-501

#SBATCH --array=501-1498

tomorrow we will run the analysis on Beluga.: 
#SBATCH --array=1499-1737


count number of files in a directory:
find . -type f | wc -l


error messages:
array:1116 error:  CpG 2778 Error in eval_f(x, ...) : Downdated VtV is not positive definite
Calls: source ... <Anonymous> -> nloptr -> <Anonymous> -> eval_f -> .Call

array: 1167 error: 2055 same error .
array: 1180, cpg: 1241
array 693 cpg:1496
array 934, cpg: 2627

conclusion: all these error were due to most methy_respoonse = 0. it will cause problem for lmer model.
we fix it by:  use try to skip the problematic cpg only for lmer. 

collect the file sizes of all output, they should all have 3000 rows therefore similar size. we just need to check the first frew rows
to see how they fail from *.out file.

stat -c "%n,%s" *simple.txt > ~/scratch/UQAC/meth/file_sizes.csv
sort -t, -k2,2n ~/scratch/UQAC/meth/file_sizes.csv > ~/scratch/UQAC/meth/sorted_file_sizes.csv


After running #SBATCH --array=1499-1737:

we write script to collect the maximum iteration of each *.simple.txt file. and save all the arrayIDs that did not complte into a txt file. 6_4_eval.R
Then run the anlaysis with  these arrayIDs in 3hours first, to resolve these arrayIDs which is due to lmer failing.
Then collect unfinished arrayIDs, and run longer hours (12hours) for those time-consumping jobs.


## 7_eval evaluation of the whole genome results


cat ~/scratch/UQAC/meth/results/whole_genome/*simple.txt > ~/scratch/UQAC/meth/results/all_simple.txt

cat ~/scratch/UQAC/meth/results/whole_genome/*complete.txt > ~/scratch/UQAC/meth/results/all_complete.txt

cd ~/scratch/UQAC/meth/results/whole_genome



7_eval.ipynb: the report summarizing results.


## 8_rerun rerun these CpGs with LASSO selection


we first collect the corresponding chunk ID and also cpg id to rerun.
 combine all *simple.txt files, prefixing each line with its filename



cd ~/scratch/UQAC/meth/results/whole_genome

for f in *simple.txt; do
  awk -v FN="$f" 'BEGIN{OFS="\t"} {print FN, $0}' "$f"
done > ~/scratch/UQAC/meth/results/rerun_results_simple.txt



1_sel.R: experimental script to select the CpGs to rerun.
2_col.R: the actual script to collect the cpgs which need to be rerun. 

3_job.R: the local script to rerun the invalid jobs on Beluga again.

4_eval.R: collet the rerun results and evaluate them. manhattanplot, and venn digaram of overlapping cpgs.

5_test.R: test script. obsolete.

6_GO.R: GO analysis of the rerun results.
7_gsea.R: GSEA analysis of the rerun results. 
###########################################################################################

March 18, 2025

The report is updated:
1_pheno.ipynb: we use the smoker variable in two levels: non smokers vs (smokers + ex-smokers). The reason is that the definition of ex-smokers is ambiguous, some will define ex-smokers as someone who stop smoking since 3 months, or 1 year, etc. Moreover, ex-smokers may have effects on respiratory functions for long time.

we will also update the phenotype defintion of AD and AA.

we also updated the cell-type proporstions with updates from AM.

2_1_sva.r: we will update the script to include the new phenotype definition.


To Do later:
a sensitive analysis to prove the choice of 10,000 probes for sva calculation is enough. maybe a figure for the choice of No_probes versus performance.


sensitivty analysis to justfy the use of arcsin-square root transofmration rather than M value.

mutliple testing. 

################################

## 9_regional regional analysis

1_test.R: test run Somnibus on a single region.


2_beluga.R: run the analysis on Beluga.

3_others.R: other methds including SMSC and GlobalTest.

4_1_sel_region.R: select cpgs based on the previous results for sominbus analysis.


5_sel_cov.R: to check other covriates effects from previous results, to narrow down for the somnibus analysis.

6_test_chunk.R: run the somnibus analysis on a single chunk, by using the sleected cpgs and selcted covariates from 4_1_sel_region.R and 5_sel_cov.R.


9_regional/beluga/*: the final version of scripts to be run on Beluga.

after running on the analysis on Beluga, we will download all somnibus results from Beluga and evaluate them.
9_regional/13_eval.R: collect all somnibus results from Beluga and analyze them.


9_regional/14_eval.R: collect all pvalues and make multiple testing correction.
9_regional/16_eval.R: collect the chr infomration and make a final report.

in the report:

A figure explaining the regional anlaysis, how its done, and how it is different from the previous whole genome association analysis.

summary statsitics; number of regions, number of cpgs, number of significant regions, number of significant cpgs, number of significant genes, number of significant pathways.


table: how many cpgs were overlapping with the previous 3 methods. 

table: the gene of the significant regions, anyone interesting of asthma association? 

GS pathway analysis: and the pathways of the significant regions.



## 10_sanity sanity check analysis

1_test.R: test run the sanity check analysis on a single region.


## 11_mgcv MGCV regional Analysis

1_correct_EMXOS.R： correct the EMXOS analysis with updated gene region fro ENSEMBL.
2_spacing.R: test split the data into regions for MGCV analysis. and also merge regions afterward is they are close enough, but located in different chunks.
then split regions into subregions if the contain too many cpgs.

beluga/7_run_region.sh: experiement with first 10 cores.
beluga/7_2_run_region_batch2.sh: run the analysis on Beluga with 300 cores. each will take about 3 hours to run.

Tomorrow: 
collect the results from Beluga, and evalute how many regions are missing.
set up a pipeline to collect the missing regions and run them again. set up scripts. 
after that, evaluate the existing results.

important:
Look for edf > 0.5 or 1 before trusting the s(start):AA_only p-value.


17_gene.R: we colelct the top genes with highest edf for s(start):AA_only, and check what phenotypes they are related to.

18_eval_pheno.R: we check the comprehenstive phenotype data from AM to see if we need to adjust fro them in the model.


## 12_cpg_DMR

This folder contains codes to re-run single CpG analysis with BMI adjustment.

Then compare CpG results and DMR.

DMR definition evaluations.

1_rerun.R: re-run the DMR definition with BMI adjustment.

2_rerun_rorqual.R: re-run it on rorqual cluster.

3_BMI_age.R: evaluate the BMI and age effect on AA status. Decide whether to include BMI and age in the model.

4_M2.R: run model 2 seprately with BMI adjustment.

5_0_comb.R: combine the M2 results with BMI adjustment, alongside with M1 and M3 without BMI adjustment.


5_1_eval.R: evaluate M1, M2, M3 results.

6_overlap.R: test the overlapping betweem M1 and DMR results.



cd ~/scratch/UQAC/meth/results/12_cpg_DMR_M2/

for f in *simple.txt; do
  awk -v FN="$f" 'BEGIN{OFS="\t"} {print FN, $0}' "$f"
done > ~/scratch/UQAC/meth/results/12_cpg_DMR_M2_simple.txt


scp laj773@rorqual.calculquebec.ca:/home/laj773/scratch/UQAC/meth/results/13_dmrcate/stringent/*.tsv /mnt/c/Per/LaiJiang/Project/UQAC/meth/results/13_dmrcate/stringent/

scp /mnt/c/Per/LaiJiang/Project/UQAC/meth/scr/8_rerun/results/2_col_rerun_M2.txt laj773@rorqual.calculquebec.ca:/home/laj773/scratch/UQAC/meth/scr/8_rerun/dat/2_col_rerun_M2.txt 


# 13_BSmooth

This folder contains codes to run BSmooth regional analysis and DMRcate.

## 14_paper
scripts to collect results and generate contents for the paper.
1_concor.R, 2_mgcv.R: the scripts for the Concordance Analysis I (Single CpG vs GAM model)
3_overlap.R: the script for Concordance Analysis II (Systematic quantification across all tested regions)
9_supplement.R: collect the supplementary tables and figures for the paper.


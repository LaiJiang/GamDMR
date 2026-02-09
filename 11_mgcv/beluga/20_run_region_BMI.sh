#!/bin/bash
#SBATCH --time=02:59:00
#SBATCH --mem=7G
#SBATCH --account=def-claprise
#SBATCH --array=1-300
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err



echo 'Start 20_run_region_BMI.R all regions for AA!'

date  # Print the current system time after starting


module purge
module load StdEnv/2023
module load gcc/12.3
module load r/4.4.0
module load r-bundle-bioconductor/3.20


Rscript  20_run_region_BMI.R $SLURM_ARRAY_TASK_ID

echo 'Finished 20_run_region_BMI.R all regions for AA!'

date  # Print the current system time after starting

sleep 30
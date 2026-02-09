#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem=7G
#SBATCH --account=def-claprise
#SBATCH --output=chr_logs/slurm-%A_%a.out
#SBATCH --error=chr_logs/slurm-%A_%a.err

echo "Start 15_extract_chr.R for to collect the chr information for the selected Somnibus regions"

date

source /home/laj773/Venv/MethEnv/bin/activate
module load r/4.4.0

export R_LIBS_USER=~/scratch/R/library

Rscript  15_extract_chr.R 

echo 'Finished 15_extract_chr.R for AA!'

date


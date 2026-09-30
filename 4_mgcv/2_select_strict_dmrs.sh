#!/bin/bash
#SBATCH --account=def-claprise
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

module load StdEnv/2023
module load gcc/12.3
module load r/4.4.0

Rscript 2_select_strict_dmrs.R


#!/bin/bash
#SBATCH --account=def-claprise
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

module load StdEnv/2023
module load gcc/12.3
module load r/4.4.0

Rscript 4_collect_mgcv_family_split_dmrs.R

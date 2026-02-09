#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --mem=24G
#SBATCH --account=def-claprise
#SBATCH --array=1-455
#SBATCH --output=StageIIIlogs/slurm-%A_%a.out
#SBATCH --error=StageIIIlogs/slurm-%A_%a.err

echo "Start 12_stage3.R for AA with region IDs (not data chunk ID) thats failed in stage 1 and stage 2"

date

source /home/laj773/Venv/MethEnv/bin/activate
module load r/4.4.0

export R_LIBS_USER=~/scratch/R/library

Rscript  12_stage3.R $SLURM_ARRAY_TASK_ID

echo 'Finished 12_stage3.R for AA!'

date


#!/bin/bash
#SBATCH --time=02:50:00
#SBATCH --mem=10G
#SBATCH --array=1-500

module load r

Rscript 6_3_job_clean.R ${SLURM_ARRAY_TASK_ID}

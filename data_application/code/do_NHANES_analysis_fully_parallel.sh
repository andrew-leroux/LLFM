#!/bin/sh
#$ -R y
#$ -t 1-1440
module load conda_R
Rscript LLFM_NHANES_manuscript_cluster_fully_parallel.R $SGE_TASK_ID $b 1001

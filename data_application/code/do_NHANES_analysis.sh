#!/bin/sh
#$ -R y
#$ -cwd
#$ -pe local 20
#$ -l mem_free=1G,h_vmem=1G
#$ -t 1-1000
module load conda_R
Rscript LLFM_NHANES_manuscript_cluster.R $SGE_TASK_ID $1000 $20

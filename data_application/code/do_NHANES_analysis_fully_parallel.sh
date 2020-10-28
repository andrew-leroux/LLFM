#!/bin/sh
#$ -R y
module load conda_R
Rscript LLFM_NHANES_manuscript_cluster_fully_parallel.R $b $s

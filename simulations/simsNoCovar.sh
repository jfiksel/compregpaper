#!/bin/bash
#$ -l mem_free=1G
#$ -l h_vmem=1G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y
#$ -R y
#$ -t 1-420
module load conda_R
Rscript simsNoCovar.R $SGE_TASK_ID
#!/bin/bash

#$ -l h_vmem=80G
#$ -pe smp 4
#$ -binding linear:4
#$ -l h_rt=12:30:00
#$ -o

source /broad/software/scripts/useuse
use .bcftools-1.9
use R-3.5


SCRIPT_PATH=/humgen/diabetes2/users/shsu/prs_project
PROJECT_PATH=/humgen/diabetes2/users/shsu/prs_project
CLUSTER_PATH=/humgen/diabetes2/users/shsu/prs_project/score_info
VCF_PATH=/humgen/diabetes2/users/shsu/prs_project/vcfexample.vcf.gz
PROJECT_NAME=project_name

source $SCRIPT_PATH/prs_pipeline.sh -s $SCRIPT_PATH -p $PROJECT_PATH -c $CLUSTER_PATH -v $VCF_PATH -n $PROJECT_NAME

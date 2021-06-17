#!/bin/bash

#$ -l h_vmem=80G
#$ -l h_rt=12:30:00
#$ -o /humgen/diabetes2/users/shsu/prs_project
#$ -e /humgen/diabetes2/users/shsu/prs_project

# Make sure you have bcftools and R modules loaded!
source /broad/software/scripts/useuse
use R-3.5
use .bcftools-1.9

SCRIPT_PATH=/prs_project/polygenic-risk-scores-pipeline
PROJECT_PATH=/prs_project
CLUSTER_PATH=/prs_project/score_info
VCF_PATH=vcf.dose.vcf.gz
PROJECT_NAME=prs_project

source $SCRIPT_PATH/prs_pipeline.sh -s $SCRIPT_PATH -p $PROJECT_PATH -c $CLUSTER_PATH -f $VCF_PATH -m $PROJECT_NAME

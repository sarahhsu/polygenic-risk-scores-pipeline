#!/bin/bash

# Make sure you have bcftools and R modules loaded!

SCRIPT_PATH=/humgen/diabetes2/users/shsu/prs_project
PROJECT_PATH=/humgen/diabetes2/users/shsu/prs_project
CLUSTER_PATH=/humgen/diabetes2/users/shsu/prs_project/score_info
VCF_PATH=/humgen/diabetes2/users/shsu/prs_project/vcf_file.vcf.gz
PROJECT_NAME=project_name

source $SCRIPT_PATH/prs_pipeline.sh -s $SCRIPT_PATH -p $PROJECT_PATH -c $CLUSTER_PATH -v $VCF_PATH -n $PROJECT_NAME

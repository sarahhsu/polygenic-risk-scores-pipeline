#!/bin/bash

# Make sure you have bcftools, tabix, and R modules loaded!

SCRIPT_PATH=/prs_project/polygenic-risk-scores-pipeline
PROJECT_PATH=/prs_project
CLUSTER_PATH=/prs_project/score_info
VCF_PATH=/prs_project/vcf_file.vcf.gz
PROJECT_NAME=project_name

source $SCRIPT_PATH/prs_pipeline.sh -s $SCRIPT_PATH -p $PROJECT_PATH -c $CLUSTER_PATH -f $VCF_PATH -m $PROJECT_NAME

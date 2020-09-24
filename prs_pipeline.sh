#!/bin/bash

#Updated September 11 2020 - Sarah Hsu
START_TIME=`date +%s`
# Default values of arguments
SCRIPT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PROJECT_PATH=$SCRIPT_PATH
SCORE_PATH=$SCRIPT_PATH
VCF_FILE_PATH="vcf"
VCF_PATH="path"
VCF_ONLY=0
PROJECT_NAME=project
SCRIPT=`basename ${BASH_SOURCE[0]}`
DATE=$(date +%F)


usage () {
  local l_MSG=$1
  echo
  echo "##########################################################################################################"
  echo
  echo "Usage Error: $l_MSG"
  echo
  echo "Usage: $SCRIPT -s </path/to/script/files> -p </path/to/project/files> -c </path/to/score/files> -f <path/to/vcf/file.vcf.gz>"
  echo "Recognized optional command line arguments:"
  echo "-n  -- Project name"
  echo "-v  -- Extract VCF files only, no PRS calculation"
  echo "-h  -- Help message"
  echo
  echo "##########################################################################################################"
  echo
  echo "END TIME: $(date +"%T")."
  echo "RUN TIME: $(expr `date +%s` - $START_TIME) seconds."
  echo
  exit 1
}

help () {
  echo
  echo "#########################################################################################################################################"
  echo
  echo "Usage: $SCRIPT -s </path/to/script/files> -p </path/to/project/files> -c </path/to/score/files> -f <path/to/vcf/file.vcf.gz>"
  echo
  echo "Recognized optional command line arguments:"
  echo "-n  -- Project name"
  echo "-v  -- Extract VCF files only, no PRS calculation"
  echo "-h  -- Help message"
  echo "------------------------------------------------------------------------------------------------------------------------------------------"
  echo "Put all 3 scripts in /path/to/script/files including (1) prs_pipeline_broad_2020-05-21.sh, (2) vcf_pipeline_broad_2020-05-21.sh, (3) prs_2020-03-09.R."
  echo
  echo "Make sure you are using FULL ABSOLUTE PATHS!"
  echo
  echo "All project files should be in /path/to/project/files including a folder called txt_files that has"
  echo "snps.txt (chr<tab>pos<tab>ref<tab>alt ie.3<tab>137844645<tab>T<tab>C)."
  echo
  echo "All files to generate scores should be in /path/to/score/files with each score in a separate csv
  file with headers Chr, Pos, Ref, Alt, RSID, Effect_Allele, Weight."
  echo
  echo "-n option is for when you want to add your project name to your output file folder and score files (DO NOT INCLUDE ANY SPACES,
  ie. polygenic_scores)."
  echo
  echo "-v option is if you don't want to generate scores, and you just want to exctract a VCF file with
  all of your listed variants. You must include your variants in a file /path/to/project/files/snps.txt"
  echo "(chr<tab>pos<tab>ref<tab>alt ie.3<tab>137844645<tab>T<tab>C)"
  echo
  echo "##########################################################################################################################################"
  echo
  echo "END TIME: $(date +"%T")"
  echo "RUN TIME: $(expr `date +%s` - $START_TIME) seconds"
  echo
  exit 1
}

echo "Polygenic Scores Pipeline: Updated September 11 2020"
echo "DATE: $DATE"
echo "START TIME: $(date +"%T")"
echo "Parsing command line arguments. Here is what your arguments were:"
for i; do
   echo $i
done

echo
### check number of command line arguments
NUMARGS=$#
if [ $NUMARGS -le 6 ]; then
  if [ "$1" == "-h" ]; then
  help
  exit 1
  fi
  usage 'Not enough command line arguments specified!'
fi
# Loop through arguments and process them
# Loop through arguments and process them
while [[ $# -gt 0 ]]
  do
    key="$1"
    case $key in
      -s) # set option "s" specifying the path to scripts
      SCRIPT_PATH="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
      -p) # set option "p" specifying the path to project
      PROJECT_PATH="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
      -c) # set option "c" specifying the path to scores
      SCORE_PATH="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
      -f) #set option "f" specifying the vcf file name
      VCF_FILE_PATH="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
      -z) #set option "z" specifying the vcf folder name
      VCF_PATH="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
      -n) #set option "n" specifying the project name
      PROJECT_NAME="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
      -v) #set option "v" for extract vcf only
      VCF_ONLY=1
      shift # Remove argument name from processing
      ;;
      *) usage "Invalid command line argument $1!"
    esac
done
wait


if [[ "$VCF_PATH" != "path" && "$VCF_FILE_PATH" != "vcf" ]].; then
  echo "Please enter path to a folder where the VCF files are (-z) OR a single VCF file (-f)."
  exit 1
fi

if [ "$PROJECT_NAME" != "project" ]; then
  PROJECT_NAME="${PROJECT_NAME}_"
  PROJECT="${PROJECT_NAME}"
fi

if [ "$PROJECT_NAME" == "project" ]; then
  PROJECT=""
fi


echo "Extracting variants from VCF file."
source $SCRIPT_PATH/vcf_pipeline.sh -p $PROJECT_PATH -n $PROJECT_NAME -f $VCF_FILE_PATH -z $VCF_PATH -s $SCORE_PATH -v $VCF_ONLY || exit 1
echo "RUN TIME: $(expr `date +%s` - $START_TIME) seconds"
if [ $VCF_ONLY -eq 1 ]; then
    echo "END TIME: $(date +"%T")"
    echo "RUN TIME: $(expr `date +%s` - $START_TIME) seconds"
    return
fi

wait
echo
echo "Calculating polygenic risk scores."
Rscript --vanilla $SCRIPT_PATH/prs.R $PROJECT_PATH $SCORE_PATH $PROJECT_PATH/${PROJECT}final_VCF_${USER}_${DATE} $USER $DATE $PROJECT_NAME || exit 1

echo "END TIME: $(date +"%T")."
echo "RUN TIME: $(expr `date +%s` - $START_TIME) seconds."

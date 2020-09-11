#!/bin/bash

#Updated September 11 2020 - Sarah Hsu

#Variables
PROJECT_FOLDER="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DATE=$(date +%F)

# Loop through arguments and process them
for arg in "$@"
do
    case $arg in
      -p|--project)
      PROJECT_FOLDER="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
      -n|--name)
      PROJECT="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
      -f|--file)
      VCF_PATH="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
      -s|--score)
      SCORE_PATH="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
      -v|--vcf)
      SCORE_PATH=$2
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
      *)
    esac
done

if [ "$PROJECT" == "project" ]; then
  PROJECT=""
fi

FINAL_VCF_FOLDER=$PROJECT_FOLDER/${PROJECT}final_VCF_${USER}_${DATE}
RESULTS_FOLDER=$PROJECT_FOLDER/${PROJECT}results_${USER}_${DATE}

echo
echo "Extracting SNPs from VCF file ."
if [ ! -f $PROJECT_FOLDER/snps.txt ]; then
    if [ $VCF_ONLY -eq 1 ]; then
      echo "$PROJECT_FOLDER/snps.txt not found! Exiting."
      exit 1
    fi
    echo "$PROJECT_FOLDER/snps.txt not found! Creating from score files."
    awk 'FNR==1 && NR!=1{next;}{print}' $SCORE_PATH/*.csv > $PROJECT_FOLDER/temp_snps.txt
    echo "$(tail -n +2 $PROJECT_FOLDER/temp_snps.txt)" > $PROJECT_FOLDER/temp_snps.txt
    cat $PROJECT_FOLDER/temp_snps.txt | tr "," "\\t" > snps.txt
    mv snps.txt $PROJECT_FOLDER/temp_snps.txt
    cut -f 1-4 $PROJECT_FOLDER/temp_snps.txt > $PROJECT_FOLDER/snps.txt
    mv $PROJECT_FOLDER/snps.txt $PROJECT_FOLDER/temp_snps.txt
    sort $PROJECT_FOLDER/temp_snps.txt | uniq -u > $PROJECT_FOLDER/snps.txt
    rm $PROJECT_FOLDER/temp_snps.txt
fi

mkdir $FINAL_VCF_FOLDER
cd $FINAL_VCF_FOLDER
filename=$(basename $VCF_PATH .vcf.gz)
bcftools view -O v -R $PROJECT_FOLDER/snps.txt $VCF_PATH  | bcftools annotate --output-type z --output ${filename}_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'


echo
zgrep -v "^##" $FINAL_VCF_FOLDER/${filename}_snps_${USER}_${DATE}.recode.vcf.gz | cut -f1-5 > $RESULTS_FOLDER/${PROJECT}extracted_snps_${USER}_${DATE}.txt
COUNT=$(wc -l $RESULTS_FOLDER/${PROJECT}extracted_snps_${USER}_${DATE}.txt)
echo "SNPs successfully extracted located at $RESULTS_FOLDER/${PROJECT}extracted_snps_${USER}_${DATE}.txt ."
echo "Line count of SNP file (with header): $COUNT ."
echo

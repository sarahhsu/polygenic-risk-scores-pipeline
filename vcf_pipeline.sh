#!/bin/bash

#Updated November 13 2020 - Sarah Hsu

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
      -m|--name)
      PROJECT="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
      -f|--file)
      VCF_FILE_PATH="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
      -z|--folder)
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
      VCF_ONLY=$2
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
      *)
    esac
done

if [ "$PROJECT" == "project" ]; then
  PROJECT=""
fi

MERGED_VCF_FOLDER=$PROJECT_FOLDER/merged_VCF
FINAL_VCF_FOLDER=$PROJECT_FOLDER/${PROJECT}final_VCF_${USER}_${DATE}
RESULTS_FOLDER=$PROJECT_FOLDER/${PROJECT}results_${USER}_${DATE}

echo
if [ ! -f $PROJECT_FOLDER/snps.txt ]; then
    if [ $VCF_ONLY -eq 1 ]; then
      echo "$PROJECT_FOLDER/snps.txt not found! Exiting."
      exit 1
    fi
    echo "$PROJECT_FOLDER/snps.txt not found! Creating from score files."
    awk 'FNR==1 && NR!=1{next;}{print}' $SCORE_PATH/*.csv > $PROJECT_FOLDER/temp_snps.txt
    echo "$(tail -n +2 $PROJECT_FOLDER/temp_snps.txt)" > $PROJECT_FOLDER/temp_snps.txt
    cat $PROJECT_FOLDER/temp_snps.txt | tr "," "\\t" > $PROJECT_FOLDER/snps.txt
    mv $PROJECT_FOLDER/snps.txt $PROJECT_FOLDER/temp_snps.txt
    cut -f 1-4 $PROJECT_FOLDER/temp_snps.txt > $PROJECT_FOLDER/snps.txt
    mv $PROJECT_FOLDER/snps.txt $PROJECT_FOLDER/temp_snps.txt
    cat -n $PROJECT_FOLDER/temp_snps.txt  | sort -uk2 | sort -nk1 | cut -f2- > $PROJECT_FOLDER/snps.txt
    rm $PROJECT_FOLDER/temp_snps.txt
fi

if [ "$VCF_PATH" != "path" ] # folder
then
  #Extract SNPs from vcf files
  echo
  echo "Extracting SNPs from VCF files"
  mkdir $MERGED_VCF_FOLDER
  cd $MERGED_VCF_FOLDER
  for file in $VCF_PATH/*.vcf.gz; do
    filename=$(basename $file .vcf.gz)
    gunzip -c $file | awk '{gsub(/^chr/,""); print}' > temp.vcf
    bgzip temp.vcf
    tabix -p vcf temp.vcf.gz
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt temp.vcf.gz | bcftools annotate --output-type z --output ${filename}_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
    tabix -p ${filename}_snps_${USER}_${DATE}.recode.vcf.gz
    rm *temp.vcf.gz*
  done

  wait

  mkdir $FINAL_VCF_FOLDER
  cd $FINAL_VCF_FOLDER
  echo
  echo "Merging files saving into file $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz"
  bcftools concat -a $MERGED_VCF_FOLDER/*.recode.vcf.gz -Oz -o $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz
  tabix -p vcf $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz

  cd
  #rm -rf $PROJECT_FOLDER/merged_VCF

else # one file
  mkdir $FINAL_VCF_FOLDER
  cd $FINAL_VCF_FOLDER
  echo
  echo "Extracting SNPS and saving into file $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz"
  tabix -p vcf $VCF_FILE_PATH
  bcftools view -O v -R $PROJECT_FOLDER/snps.txt $VCF_FILE_PATH  | bcftools annotate --output-type z --output ${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  cd
fi



echo
mkdir $RESULTS_FOLDER
zgrep -v "^##" $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz | cut -f1-5 > $RESULTS_FOLDER/${PROJECT}extracted_snps_${USER}_${DATE}.txt
COUNT=$(wc -l $RESULTS_FOLDER/${PROJECT}extracted_snps_${USER}_${DATE}.txt)
echo "SNPs successfully extracted located at $RESULTS_FOLDER/${PROJECT}extracted_snps_${USER}_${DATE}.txt ."
echo "Line count of SNP file (with header): $COUNT ."
echo

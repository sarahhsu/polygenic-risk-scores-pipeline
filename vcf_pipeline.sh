#!/bin/bash

#Updated February 15 2021 - Sarah Hsu

#Variables
PROJECT_FOLDER="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DATE=$(date +%F)
BROAD_FILES=/humgen/diabetes/users/josep/PartnersBIOBANK/merge_all_35K_datasets/merge_35K_genotypes/
ERISONE_FILES=/data/dgag/raw_data/Biobank_akm41_20190625_115053_x31582_Imputed/biobank_imputed_subset/
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
      -b|--biobank)
      VCF_SOURCE="$2"
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
    sort $PROJECT_FOLDER/temp_snps.txt | uniq -u > $PROJECT_FOLDER/snps.txt
    rm $PROJECT_FOLDER/temp_snps.txt
fi

if [ "$VCF_PATH" != "path" ] # folder
then
  if [ ! -d $VCF_PATH ]; then
        echo "VCF files not found! Are you sure the path is correct?"
        exit 1
  fi
  #Extract SNPs from vcf files
  echo
  echo "Extracting SNPs from VCF files"
  mkdir $MERGED_VCF_FOLDER
  cd $MERGED_VCF_FOLDER
  for file in $VCF_PATH/*.vcf.gz; do
    filename=$(basename $file .vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output ${filename}_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  done

  wait

  mkdir $FINAL_VCF_FOLDER
  cd $FINAL_VCF_FOLDER
  echo
  echo "Merging files saving into file $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz"
  bcftools concat $MERGED_VCF_FOLDER/*.recode.vcf.gz -Oz -o $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz

  cd
  rm -rf $MERGED_VCF_FOLDER

elif [ "$VCF_FILE_PATH" != "vcf" ] # one file
then
  if [ ! -f $VCF_FILE_PATH ]; then
        echo "VCF file not found! Are you sure the path is correct?"
        exit 1
  fi
  mkdir $FINAL_VCF_FOLDER
  cd $FINAL_VCF_FOLDER
  echo
  echo "Extracting SNPS and saving into file $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz"
  bcftools view -O v -R $PROJECT_FOLDER/snps.txt $VCF_FILE_PATH  | bcftools annotate --output-type z --output ${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  cd

elif [ "$VCF_SOURCE" != "broad-eu" ] # broad europeans
then
  if [ ! -f $BROAD_FILES/EU/imputation_TOPMED_R2_RELATED_EU/subset_1/chr1.dose.vcf.gz ]; then
        echo "Europeans files not found! Are you sure you're on the Broad server?"
        exit 1
  fi
  #Extract SNPs from vcf files
  echo
  echo "Extracting SNPs from MGB Biobank files on the Broad server"
  mkdir $MERGED_VCF_FOLDER
  cd $MERGED_VCF_FOLDER
  for file in $BROAD_FILES/EU/imputation_TOPMED_R2_RELATED_EU/subset_1/chr*.dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output ${filename}_1_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  done
  for file in $BROAD_FILES/EU/imputation_TOPMED_R2_RELATED_EU/subset_2/chr*.dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output ${filename}_2_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
    bcftools merge $MERGED_VCF_FOLDER/${filename}_1_snps_${USER}_${DATE}.recode.vcf.gz $MERGED_VCF_FOLDER/${filename}_2_snps_${USER}_${DATE}.recode.vcf.gz -Oz -o $MERGED_VCF_FOLDER/${filename}_snps_${USER}_${DATE}.recode.vcf.gz
    rm $MERGED_VCF_FOLDER/${filename}_1_snps_${USER}_${DATE}.recode.vcf.gz
    rm $MERGED_VCF_FOLDER/${filename}_2_snps_${USER}_${DATE}.recode.vcf.gz
  done

  bcftools view -O v -R $PROJECT_FOLDER/snps.txt $BROAD_FILES/EU/imputation_TOPMED_R2_RELATED_EU/X_females/chrX.dose.vcf.gz | bcftools annotate --output-type z --output female_chrX_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  bcftools view -O v -R $PROJECT_FOLDER/snps.txt $BROAD_FILES/EU/imputation_TOPMED_R2_RELATED_EU/X_males/chrX.dose.vcf.gz | bcftools annotate --output-type z --output male_chrX_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'

  bcftools merge $MERGED_VCF_FOLDER/male_chrX_snps_${USER}_${DATE}.recode.vcf.gz $MERGED_VCF_FOLDER/female_chrX_snps_${USER}_${DATE}.recode.vcf.gz -Oz -o $MERGED_VCF_FOLDER/chrX_snps_${USER}_${DATE}.recode.vcf.gz

  rm $MERGED_VCF_FOLDER/male_chrX_snps_${USER}_${DATE}.recode.vcf.gz
  rm $MERGED_VCF_FOLDER/female_chrX_snps_${USER}_${DATE}.recode.vcf.gz
  wait

  mkdir $FINAL_VCF_FOLDER
  cd $FINAL_VCF_FOLDER
  echo
  echo "Merging files saving into file $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz"
  bcftools concat $MERGED_VCF_FOLDER/*.recode.vcf.gz -Oz -o $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz

  cd
  rm -rf $MERGED_VCF_FOLDER

elif [ "$VCF_SOURCE" != "broad-afr" ] # broad africans
then
  if [ ! -f $BROAD_FILES/AFR/imputation_TOPMED_R2_RELATED_AFR/chr1.dose.vcf.gz ]; then
        echo "Africans files not found! Are you sure you're on the Broad server?"
        exit 1
  fi
  #Extract SNPs from vcf files
  echo
  echo "Extracting SNPs from MGB Biobank files on the Broad server"
  mkdir $MERGED_VCF_FOLDER
  cd $MERGED_VCF_FOLDER
  for file in $BROAD_FILES/AFR/imputation_TOPMED_R2_RELATED_AFR/chr*.dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output ${filename}_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  done

  bcftools view -O v -R $PROJECT_FOLDER/snps.txt $BROAD_FILES/AFR/imputation_TOPMED_R2_RELATED_AFR/X_females_related/chrX.dose.vcf.gz | bcftools annotate --output-type z --output female_chrX_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  bcftools view -O v -R $PROJECT_FOLDER/snps.txt $BROAD_FILES/AFR/imputation_TOPMED_R2_RELATED_AFR/X_males_related/chrX.dose.vcf.gz | bcftools annotate --output-type z --output male_chrX_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'

  bcftools merge $MERGED_VCF_FOLDER/male_chrX_snps_${USER}_${DATE}.recode.vcf.gz $MERGED_VCF_FOLDER/female_chrX_snps_${USER}_${DATE}.recode.vcf.gz -Oz -o $MERGED_VCF_FOLDER/chrX_snps_${USER}_${DATE}.recode.vcf.gz

  rm $MERGED_VCF_FOLDER/male_chrX_snps_${USER}_${DATE}.recode.vcf.gz
  rm $MERGED_VCF_FOLDER/female_chrX_snps_${USER}_${DATE}.recode.vcf.gz
  wait

  mkdir $FINAL_VCF_FOLDER
  cd $FINAL_VCF_FOLDER
  echo
  echo "Merging files saving into file $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz"
  bcftools concat $MERGED_VCF_FOLDER/*.recode.vcf.gz -Oz -o $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz

  cd
  rm -rf $MERGED_VCF_FOLDER
elif [ "$VCF_SOURCE" != "broad-his" ] # broad hispanics
then
  if [ ! -f $BROAD_FILES/HIS/imputation_TOPMED_R2_RELATED_HIS/chr1.dose.vcf.gz ]; then
        echo "Hispanics files not found! Are you sure you're on the Broad server?"
        exit 1
  fi
  #Extract SNPs from vcf files
  echo
  echo "Extracting SNPs from MGB Biobank files on the Broad server"
  mkdir $MERGED_VCF_FOLDER
  cd $MERGED_VCF_FOLDER
  for file in $BROAD_FILES/HIS/imputation_TOPMED_R2_RELATED_HIS/chr*.dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output ${filename}_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  done

  bcftools view -O v -R $PROJECT_FOLDER/snps.txt $BROAD_FILES/HIS/imputation_TOPMED_R2_RELATED_HIS/X_females_related/chrX.dose.vcf.gz | bcftools annotate --output-type z --output female_chrX_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  bcftools view -O v -R $PROJECT_FOLDER/snps.txt $BROAD_FILES/HIS/imputation_TOPMED_R2_RELATED_HIS/X_males_related/chrX.dose.vcf.gz | bcftools annotate --output-type z --output male_chrX_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'

  bcftools merge $MERGED_VCF_FOLDER/male_chrX_snps_${USER}_${DATE}.recode.vcf.gz $MERGED_VCF_FOLDER/female_chrX_snps_${USER}_${DATE}.recode.vcf.gz -Oz -o $MERGED_VCF_FOLDER/chrX_snps_${USER}_${DATE}.recode.vcf.gz

  rm $MERGED_VCF_FOLDER/male_chrX_snps_${USER}_${DATE}.recode.vcf.gz
  rm $MERGED_VCF_FOLDER/female_chrX_snps_${USER}_${DATE}.recode.vcf.gz
  wait

  mkdir $FINAL_VCF_FOLDER
  cd $FINAL_VCF_FOLDER
  echo
  echo "Merging files saving into file $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz"
  bcftools concat $MERGED_VCF_FOLDER/*.recode.vcf.gz -Oz -o $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz

  cd
  rm -rf $MERGED_VCF_FOLDER
elif [ "$VCF_SOURCE" != "broad-all" ] # all broad files
then
  if [ ! -f $BROAD_FILES/EU/imputation_TOPMED_R2_RELATED_EU/subset_1/chr1.dose.vcf.gz ]; then
        echo "Biobank files not found! Are you sure you're on the Broad server?"
        exit 1
  fi
  echo
  echo "Extracting SNPs from MGB Biobank files on the Broad server"
  mkdir $MERGED_VCF_FOLDER
  cd $MERGED_VCF_FOLDER

  #Extract SNPs from EU files
  for file in $BROAD_FILES/EU/imputation_TOPMED_R2_RELATED_EU/subset_1/chr*.dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output EU_${filename}_1_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  done
  for file in $BROAD_FILES/EU/imputation_TOPMED_R2_RELATED_EU/subset_2/chr*.dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output EU_${filename}_2_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
    bcftools merge $MERGED_VCF_FOLDER/EU_${filename}_1_snps_${USER}_${DATE}.recode.vcf.gz $MERGED_VCF_FOLDER/EU_${filename}_2_snps_${USER}_${DATE}.recode.vcf.gz -Oz -o $MERGED_VCF_FOLDER/$EU_{filename}_snps_${USER}_${DATE}.recode.vcf.gz
    rm $MERGED_VCF_FOLDER/${filename}_1_snps_${USER}_${DATE}.recode.vcf.gz
    rm $MERGED_VCF_FOLDER/${filename}_2_snps_${USER}_${DATE}.recode.vcf.gz
  done

  bcftools view -O v -R $PROJECT_FOLDER/snps.txt $BROAD_FILES/EU/imputation_TOPMED_R2_RELATED_EU/X_females/chrX.dose.vcf.gz | bcftools annotate --output-type z --output EU_female_chrX_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  bcftools view -O v -R $PROJECT_FOLDER/snps.txt $BROAD_FILES/EU/imputation_TOPMED_R2_RELATED_EU/X_males/chrX.dose.vcf.gz | bcftools annotate --output-type z --output EU_male_chrX_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'

  bcftools merge $MERGED_VCF_FOLDER/EU_male_chrX_snps_${USER}_${DATE}.recode.vcf.gz $MERGED_VCF_FOLDER/EU_female_chrX_snps_${USER}_${DATE}.recode.vcf.gz -Oz -o $MERGED_VCF_FOLDER/EU_chrX_snps_${USER}_${DATE}.recode.vcf.gz

  rm $MERGED_VCF_FOLDER/EU_male_chrX_snps_${USER}_${DATE}.recode.vcf.gz
  rm $MERGED_VCF_FOLDER/EU_female_chrX_snps_${USER}_${DATE}.recode.vcf.gz

  #Extract SNPs from AFR files
  for file in $BROAD_FILES/AFR/imputation_TOPMED_R2_RELATED_AFR/chr*.dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output AFR_${filename}_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  done

  bcftools view -O v -R $PROJECT_FOLDER/snps.txt $BROAD_FILES/AFR/imputation_TOPMED_R2_RELATED_AFR/X_females_related/chrX.dose.vcf.gz | bcftools annotate --output-type z --output AFR_female_chrX_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  bcftools view -O v -R $PROJECT_FOLDER/snps.txt $BROAD_FILES/AFR/imputation_TOPMED_R2_RELATED_AFR/X_males_related/chrX.dose.vcf.gz | bcftools annotate --output-type z --output AFR_male_chrX_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'

  bcftools merge $MERGED_VCF_FOLDER/AFR_male_chrX_snps_${USER}_${DATE}.recode.vcf.gz $MERGED_VCF_FOLDER/AFR_female_chrX_snps_${USER}_${DATE}.recode.vcf.gz -Oz -o $MERGED_VCF_FOLDER/AFR_chrX_snps_${USER}_${DATE}.recode.vcf.gz

  rm $MERGED_VCF_FOLDER/AFR_male_chrX_snps_${USER}_${DATE}.recode.vcf.gz
  rm $MERGED_VCF_FOLDER/AFR_female_chrX_snps_${USER}_${DATE}.recode.vcf.gz
  wait


#Extract SNPs from HIS files
  for file in $BROAD_FILES/HIS/imputation_TOPMED_R2_RELATED_HIS/chr*.dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output HIS_${filename}_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  done

  bcftools view -O v -R $PROJECT_FOLDER/snps.txt $BROAD_FILES/HIS/imputation_TOPMED_R2_RELATED_HIS/X_females_related/chrX.dose.vcf.gz | bcftools annotate --output-type z --output HIS_female_chrX_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  bcftools view -O v -R $PROJECT_FOLDER/snps.txt $BROAD_FILES/HIS/imputation_TOPMED_R2_RELATED_HIS/X_males_related/chrX.dose.vcf.gz | bcftools annotate --output-type z --output HIS_male_chrX_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'

  bcftools merge $MERGED_VCF_FOLDER/HIS_male_chrX_snps_${USER}_${DATE}.recode.vcf.gz $MERGED_VCF_FOLDER/HIS_female_chrX_snps_${USER}_${DATE}.recode.vcf.gz -Oz -o $MERGED_VCF_FOLDER/HIS_chrX_snps_${USER}_${DATE}.recode.vcf.gz

  rm $MERGED_VCF_FOLDER/HIS_male_chrX_snps_${USER}_${DATE}.recode.vcf.gz
  rm $MERGED_VCF_FOLDER/HIS_female_chrX_snps_${USER}_${DATE}.recode.vcf.gz
  wait

  for i in {1..22}; do
    bcftools merge --force-samples $MERGED_VCF_FOLDER/*chr${i}_snps_${USER}_${DATE}.recode.vcf.gz -Oz -o $MERGED_VCF_FOLDER/final_chr${i}_snps_${USER}_${DATE}.recode.vcf.gz
  done

  bcftools merge --force-samples $MERGED_VCF_FOLDER/*chrX_snps_${USER}_${DATE}.recode.vcf.gz -Oz -o $MERGED_VCF_FOLDER/final_chrX_snps_${USER}_${DATE}.recode.vcf.gz

  mkdir $FINAL_VCF_FOLDER
  cd $FINAL_VCF_FOLDER
  echo
  echo "Merging files saving into file $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz"
  bcftools concat $MERGED_VCF_FOLDER/final*.recode.vcf.gz -Oz -o $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz

  cd
  rm -rf $MERGED_VCF_FOLDER
else
  if [ ! -f $ERISONE_FILES/MEGA/chr1.dose.vcf.gz ]; then
        echo "Biobank files not found! Are you sure you're on the Broad server?"
        exit 1
  fi

  echo
  echo "Extracting SNPs from MGB Biobank files on ErisOne"
  mkdir $MERGED_VCF_FOLDER
  cd $MERGED_VCF_FOLDER

  for file in $ERISONE_FILES/MEGA/*dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output MEGA_${filename}_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  done &

  for file in $ERISONE_FILES/MEG_A1_A/*dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output MEG_A1_A_${filename}_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  done &

  for file in $ERISONE_FILES/MEG_A1_B/*dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output MEG_A1_B_${filename}_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  done &

  for file in $ERISONE_FILES/MEG_C/*dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output MEG_C_${filename}_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  done &

  for file in $ERISONE_FILES/MEG_D*dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output MEG_D_${filename}_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  done &

  for file in $ERISONE_FILES/MEGAEX/*dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output MEGAEX_${filename}_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  done &

  for file in $ERISONE_FILES/PERLIS/*dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output PERLIS_${filename}_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  done &

  for file in /data/dgag/raw_data/biobank_imputed_subset/data/MEG_E/*dose.vcf.gz; do
    filename=$(basename $file .dose.vcf.gz)
    bcftools view -O v -R $PROJECT_FOLDER/snps.txt $file  | bcftools annotate --output-type z --output MEG_E_${filename}_snps_${USER}_${DATE}.recode.vcf.gz  -I '%CHROM:%POS:%REF:%ALT'
  done &

  wait

  mkdir $FINAL_VCF_FOLDER
  cd $FINAL_VCF_FOLDER
  echo
  echo "Merging files saving into file $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz"
  bcftools concat $MERGED_VCF_FOLDER/*.recode.vcf.gz -Oz -o $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz

  cd
  rm -rf $MERGED_VCF_FOLDER

fi


echo
mkdir $RESULTS_FOLDER
zgrep -v "^##" $FINAL_VCF_FOLDER/${PROJECT}ALL_snps_${USER}_${DATE}.recode.vcf.gz | cut -f1-5 > $RESULTS_FOLDER/${PROJECT}extracted_snps_${USER}_${DATE}.txt
COUNT=$(wc -l $RESULTS_FOLDER/${PROJECT}extracted_snps_${USER}_${DATE}.txt)
echo "SNPs successfully extracted located at $RESULTS_FOLDER/${PROJECT}extracted_snps_${USER}_${DATE}.txt ."
echo "Line count of SNP file (with header): $COUNT ."
echo

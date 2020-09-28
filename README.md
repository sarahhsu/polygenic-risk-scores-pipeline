# Polygenic Risk Score Calculation from VCF Files
Updated September 24, 2020

**MAKE SURE ALL VARIANT COORDINATES AND REF/ALT ARE CORRECT IN YOUR SCORE FILES!**


## RUNNING THE PIPELINE:

  * To run pipeline:
  
    1. Clone the gitlab repository by typing these commands after logging into your server:
    
        ```
        mkdir your_prs_scripts_folder
        cd your_prs_scripts_folder
        git clone https://github.com/sarahhsu/polygenic-risk-scores-pipeline.git
        ```
        
    2. Make sure all files are present in the *polygenic-risk-scores-pipeline* folder (check FILES section below for list of files).

    3. Create a copy of the prs_pipeline_submission.sh script in your project folder, rename to prs_pipeline_submission_{username}_{date}.sh and modify parameters.

    4. Submit prs_pipeline_submission_{username}_{date}.sh to your job scheduler.
    
  * To update your version of the pipeline:
   
    1. Go to your version of the *polygenic-risk-scores-pipeline* folder 
    
    2. Type in ```git pull```
    
    3. After updating, your version will now match the most up-to-date version in Github
  
--------------------------------------------------------------------------------------------------------------------------------------
 
### *Usage: source prs_pipeline.sh -s </path/to/scripts/polygenic-risk-scores-pipeline> -p </path/to/project/files> -c </path/to/score/files> [options]*

### Required command line arguments:
    -s -- Path to where the scripts are (should be in the folder named polygenic-risk-scores-pipeline)
    -p -- Path to the main project folder, where the results folders and files will be written to 
    -c -- Path to scoring files, typically a subfolder within the project folder
    -f -- Path to a gzipped VCF file (extension *.vcf.gz), either -f or -z must be used but not both
    -z -- Path to a folder of gzipped VCF files, either -f or -z must be used but not both

### Recognized optional command line arguments:
    -n  -- Add project name to the output files (no spaces)
    -v  -- Extract VCF file only, no PRS calculation
    -h  -- Help message
 
--------------------------------------------------------------------------------------------------------------------------------------
## FOR HELP:

  To print help message, run *./prs_pipeline.sh -h* in the *polygenic-risk-scores-pipeline* folder.
  
## FILES:
  
  * Files included in *polygenic-risk-scores-pipeline*:
    * prs_pipeline_submission.sh - calls the pipeline (needs modification to be able to use on your server)
    * prs_pipeline.sh - main wrapper for all scripts
    * vcf_pipeline.sh - extracts the vcf file based on variants
    * prs.R - calulates the PRS
    * example/snps.txt - variants, one per line separated by tabs, chr <tab> pos <tab> ref <tab> alt (ex. 3 <tab> 137844645 <tab> T <tab> C)
    * example/prs_pipeline_submission_shsu_2020_09-11.sh - example of a submission file for Broad server
    * example/prs_project_final_VCF_shsu_2020-09-11
      * chr22.dose_snps_shsu_2020-09-11.recode.vcf.gz - example of final extracted VCF file
    * example/prs_project_results_shsu_2020-09-11
      * prs_project_extracted_snps_shsu_2020-09-11.txt - example of list of extracted SNPs
      * prs_project_scores_shsu_2020-09-11.csv - example of scores table results (typically used for downstream analysis)
      * prs_project_scaled_scores_shsu_2020-09-11.csv - example of scaled scores table results
      * prs_project_scores_stats_shsu_2020-09-11.csv - example of scores table results with mean, std, max, min for each score
    * example/score_info:
      * score_file.csv - examples of what the risk score information csv should look like, with headers Chr, Pos, Ref, Alt, RSID, Effect_Allele, Weight.
      
## MORE DETAIL:

* Dosages are used for score calculation.

* All 3 scripts should be in /path/to/script/files including (1) prs_pipeline.sh, (2) vcf_pipeline.sh, (3) prs.R

* Make sure you are using **FULL ABSOLUTE PATHS**! 

* All project files should be in /path/to/project/files
  
* All files to generate scores should be in /path/to/score/files with each score in a separate csv file with headers Chr, Pos, Ref, Alt, RSID, Effect_Allele, Weight

* -n option is for when you want to add your project name to your output file folder and score files (DO NOT INCLUDE ANY SPACES, ie. polygenic_scores)

* -v option is if you don't want to generate scores, and you just want to exctract a VCF file with all of your listed variants. You must include your variants in a file /path/to/project/files/snps.txt
  * snps.txt is tab delimited: chr <tab> pos <tab> ref <tab> alt 
  * ie. 3 <tab> 137844645 <tab> T <tab> C
  
 --------------------------------------------------------------------------------------------------------------------------------------



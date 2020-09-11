# Polygenic Risk Score Calculation from a VCF File
Updated September 11, 2020

**MAKE SURE ALL VARIANT COORDINATES AND REF/ALT ARE CORRECT IN YOUR SCORING FILES!**


## RUNNING THE SCRIPTS:

  * To run pipeline:
  
    1. Clone the gitlab repository by typing these commands after logging into your server:
    
        ```
        mkdir your_prs_scripts_folder
        cd your_prs_scripts_folder
        git clone https://gitlab.partners.org/so454/polygenic-risk-scores-pipeline.git
        ```
        
    2. Make sure all files are present (check FILES section below for list of files).

    3. Create a copy of the prs_pipeline_submission.sh script, rename to prs_pipeline_submission_{username}_{date}.sh and modify parameters.

    4. Submit prs_pipeline_submission.sh to your job scheduler.

## FOR HELP:

  To print help message, run *./prs_pipeline.sh -h* in your scripts folder.

## FILES:
  
  * Files included:
    * prs_pipeline_submission.sh - calls the pipeline (needs modification to be able to use on your server)
    * prs_pipeline.sh - main wrapper for all scripts
    * vcf_pipeline.sh - extracts the vcf file based on variants
    * prs.R - calulates the PRS
    * example/snps.txt - variants, one per line separated by tabs, chr <tab> pos <tab> ref <tab> alt (ex. 3 <tab> 137844645 <tab> T <tab> C)
    * example/prs_pipeline_submission_shsu_2020_05_21.sh: example of a submission file for Broad server
    * example/score_info:
      * cluster1.csv - examples of what the risk score information csv should look like, with headers Chr, Pos, Ref, Alt, RSID, Effect_Allele, Weight.
      

 --------------------------------------------------------------------------------------------------------------------------------------
 
### *Usage: source prs_pipeline.sh -s </path/to/script/files> -p </path/to/project/files> -c </path/to/score/files> -f <path/to/vcf/file.vcf.gz>*

### Recognized optional command line arguments:
    -n  -- Add project name to the output files (no spaces)
    -v  -- Extract VCF file only, no PRS calculation
    -h  -- Help message
 
 --------------------------------------------------------------------------------------------------------------------------------------
## MORE DETAIL:

* All 3 scripts should be in /path/to/script/files including (1) prs_pipeline.sh, (2) vcf_pipeline.sh, (3) prs.R.

* Make sure you are using **FULL ABSOLUTE PATHS**! 

* All project files should be in /path/to/project/files. 
  
* All files to generate scores should be in /path/to/score/files with each score in a separate csv file with headers Chr, Pos, Ref, Alt, RSID, Effect_Allele, Weight.

* -n option is for when you want to add your project name to your output file folder and score files (DO NOT INCLUDE ANY SPACES, ie. polygenic_scores)

* -v option is if you don't want to generate scores, and you just want to exctract a VCF file with all of your listed variants. You must include your variants in a file /path/to/project/files/txt_files/snps.txt
  * snps.txt is tab delimited: chr <tab> pos <tab> ref <tab> alt 
  * ie. 3 <tab> 137844645 <tab> T <tab> C
  
 --------------------------------------------------------------------------------------------------------------------------------------



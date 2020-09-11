#!/usr/bin/env Rscript
library(vcfR)
#Updated September 11 2020 - Sarah Hsu

# Input and Output Description --------------------------------------------------------

#Input files:
### 1. GZIPPED VCF with extracted SNPs and individuals
### 2. Score info files

#Output files:
### 1. scores_stats.csv: UNSCALED scores for each person with mean, max, min, stdev added as rows at the bottom
### 2. scores.csv: UNSCALED scores for each person
### 3. scaled_scores.csv: SCALED scores for each person (standardized)

# Variables to be added ---------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

#Add main path where files are located and will save to
project_path <- args[1]
#Add path where cluster CSV files
cluster_path <- args[2]
#Add name of gzipped vcf file path
vcf_file_path <- args[3]
username <- args[4]
cur_date <- args[5]
project_name <- args[6]

if(project_name == "project"){
  project_name <- ""
}

# Functions ---------------------------------------------------------------

convert <- function(r, c, cur_snps, cluster, vcf_info){
  snp <- colnames(cur_snps)[c]
  ids <- paste(cluster$Chr, cluster$Pos, cluster$Ref, cluster$Alt, sep=":")
  if(snp == "NA") {
    print(paste("ERROR:", ids[c], "was not extracted."))
  }
  letter <- cluster[which(ids == snp), "Effect_Allele"]
  dose <- cur_snps[r,c]
  ref_alt <- vcf_info[which(vcf_info[,3]== snp), c("REF", "ALT")]
  if(as.character(ref_alt$REF) != as.character(cluster[which(ids == snp), c("Ref")])|
     as.character(ref_alt$ALT) != as.character(cluster[which(ids == snp), c("Alt")])){
    return(0)
  }
  if(as.character(ref_alt$ALT)==as.character(letter)){
    return(as.numeric(dose))
  } else {
    return(2-as.numeric(dose))
  }
}

# Read in files -----------------------------------------------------------

#GZIPPED VCFs with extracted SNPs and individuals

print("Reading in VCF files.")

vcf_path <- paste0("_snps_", username, "_", cur_date, ".recode.vcf.gz$")
vcf_file_names <- list.files(path=vcf_file_path, pattern=vcf_path, full.names = TRUE)


vcf_names <- c()

for (i in 1:length(vcf_file_names)) {
  name <- gsub(".vcf.gz", "", vcf_file_names[i])
  assign(name, read.vcfR(vcf_file_names[i]))
  vcf_names <- c(vcf_names, name)
}

print("Reading in score info files.")
cluster_files <- list.files(path = cluster_path, pattern="*.csv")
if (length(cluster_files) == 0){
  stop("No input score info files. Score info files must be specified in csv format with header
  Chr, Pos, Ref, Alt, RSID, Effect_Allele, Weight.")
}

cluster_names <- c()
snps <- data.frame(matrix(, nrow = 0, ncol = 4))
for (file in cluster_files){
  cluster_name <- gsub(".csv", "", file)
  assign(cluster_name, read.csv(paste(cluster_path, file, sep="/"), 
                                colClasses = c("numeric", "numeric", "character", "character", "character", "character", "numeric"), 
                                header=TRUE))
  cluster_names <- c(cluster_names, cluster_name)
  snps <- rbind(snps, get(cluster_name)[,c(1:4)])
}



# Filter out any unwanted SNPs --------------------------------------------
idxs_to_remove <- c()
for (i in 1:length(vcf_names)){
  assign("vcf_info", as.data.frame(getFIX(get(vcf_names[i]))))
  idx <- c()
  for (j in 1:nrow(vcf_info)){
    if(length(which(vcf_info[j, "CHROM"] == snps[,1] & vcf_info[j, "POS"] == snps[,2] & vcf_info[j, "REF"] == snps[,3] & vcf_info[j,"ALT"] == snps[,4]))==0){
      idx <- c(idx, j)
    }
  }
  assign(paste0(vcf_names[i],"_remove"), idx)
  idxs_to_remove <- c(idxs_to_remove, paste0(vcf_names[i],"_remove"))
}
# Create weights table ----------------------------------------------------

print("Calculating polygenic risk scores.")
weights <- data.frame(matrix(, nrow = 0, ncol = length(cluster_names) ))

##Loop through vcf_names
for (i in 1:length(vcf_names)){
  assign("dosage_info", as.data.frame(extract.gt(get(vcf_names[i]), element='DS')))
  assign("vcf_info", as.data.frame(getFIX(get(vcf_names[i]))))
  if(length(get(idxs_to_remove[i])) > 0) {
    remove <- get(idxs_to_remove[i])
    dosage_info <- dosage_info[-remove,]
    vcf_info <- vcf_info[-remove,]
  }
  cur_weights <- data.frame(matrix(0, nrow = ncol(dosage_info), ncol = length(cluster_names) ))
  rownames(cur_weights) <- colnames(dosage_info)
  colnames(cur_weights) <- cluster_names

  ##Loop through cluster_names
  for (j in 1:length(cluster_names)){
    cluster <- get(cluster_names[j])
    cur_snps <- t(dosage_info[match(c(paste(cluster$Chr, cluster$Pos, cluster$Ref, cluster$Alt, sep = ":")), rownames(dosage_info)),])
    cur_snps1 <- cur_snps
    for (r in 1:nrow(cur_snps)){
      for (c in 1:ncol(cur_snps)){
        cur_snps1[r, c] <- convert(r, c, cur_snps, cluster, vcf_info)
      }
    }
    if (ncol(cur_snps1) == 1) {
      cur_weights[,j] <- apply(cur_snps1, 1,as.numeric)*as.vector(cluster$Weight)
    } else {
      cur_weights[,j] <- rowSums(t(apply(cur_snps1, 1,as.numeric)*as.vector(cluster$Weight)))
    }
  }

  if (nrow(weights) < 1){
    weights <- cur_weights
  } else {
    weights <- rbind(weights, cur_weights)
  }
}

standardized_weights <- as.data.frame(apply(weights, 2, scale))
rownames(standardized_weights) <- rownames(weights)

print("Calculating summary statistics.")
stdevs <- apply(weights, 2, sd)
means <- apply(weights, 2, mean)
maxs <- apply(weights, 2, max)
mins <- apply(weights, 2, min)

# Write files -------------------------------------------------------------

print("Writing PRS files.")
dir.create(file.path(project_path, paste0(project_name, "results_", username, "_", cur_date)), showWarnings = FALSE)
results_path <- paste(project_path, paste0(project_name, "results_", username, "_", cur_date), sep = "/")

write.csv(rbind(weights, STDEV=stdevs, MEAN=means, MAX=maxs, MIN=mins), paste(results_path, paste0(project_name, "scores_stats_", username, "_", cur_date, ".csv"), sep="/"))
write.csv(standardized_weights, paste(results_path, paste0(project_name, "scaled_scores_", username, "_", cur_date, ".csv"), sep="/"))
write.csv(weights, paste(results_path, paste0(project_name, "scores_", username, "_", cur_date, ".csv"), sep="/"))

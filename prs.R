#!/usr/bin/env Rscript
library(vcfR)
#Updated September 24 2020 - Sarah Hsu

# Input and Output Description --------------------------------------------------------

#Input files:
### 1. GZIPPED VCF with extracted SNPs
### 2. Polygenic risk score info file(s)

#Output files:
### 1. scores_stats.csv: UNSCALED scores for each person with mean, max, min, stdev added as rows at the bottom
### 2. scores.csv: UNSCALED scores for each person 
### 3. scaled_scores.csv: SCALED scores for each person (standardized)

# Variables to be added ---------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

#Add main path where files are located and will save to
project_path <- args[1]
#Add path where cluster CSV files
score_path <- args[2]
#Add path to folder where the gzipped vcf file is
vcf_file_path <- args[3]
username <- args[4]
cur_date <- args[5]
project_name <- args[6]

if(project_name == "project"){
  project_name <- ""
}

# Functions ---------------------------------------------------------------

# Convert given dosage to 0, 1, 2 based on the effect allele.
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

# GZIPPED VCFs with extracted SNPs

print("Reading in VCF files.")

vcf_path <- paste0("snps_", username, "_", cur_date, ".recode.vcf.gz$")
vcf_file_names <- list.files(path=vcf_file_path, pattern=vcf_path, full.names = TRUE)

vcf_names <- c()
for (i in 1:length(vcf_file_names)) {
  name <- gsub(".vcf.gz", "", vcf_file_names[i])
  assign(name, read.vcfR(vcf_file_names[i]))
  vcf_names <- c(vcf_names, name)
}


print("Reading in score info files.")
cluster_files <- list.files(path = score_path, pattern="*.csv")
if (length(cluster_files) == 0){
  stop("No input score info files. Score info files must be specified in csv format with header
  Chr, Pos, Ref, Alt, RSID, Effect_Allele, Weight.")
}

cluster_names <- c()
snps <- data.frame(matrix(, nrow = 0, ncol = 4))
for (file in cluster_files){
  cluster_name <- gsub(".csv", "", file)
  assign(cluster_name, read.csv(paste(score_path, file, sep="/"), 
                                colClasses = c("numeric", "numeric", "character", "character", "character", "character", "numeric"), 
                                header=TRUE,
                                col.names=c("Chr", "Pos", "Ref", "Alt", "RSID", "Effect_Allele", "Weight")))
  cluster_names <- c(cluster_names, cluster_name)
  snps <- rbind(snps, get(cluster_name)[,c("Chr", "Pos", "Ref", "Alt")])
}



# Filter out any unwanted SNPs --------------------------------------------

idxs_to_remove <- c()
for (i in 1:length(vcf_names)){
  assign("vcf_info", as.data.frame(getFIX(get(vcf_names[i]))))
  idx <- c()
  for (j in 1:nrow(vcf_info)){
    if(length(which(vcf_info[j, "CHROM"] == snps[,"Chr"] & 
                    vcf_info[j, "POS"] == snps[,"Pos"] & 
                    vcf_info[j, "REF"] == snps[,"Ref"] & 
                    vcf_info[j,"ALT"] == snps[,"Alt"])) == 0){
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
  cur_weights <- data.frame(matrix(0, nrow = ncol(dosage_info), ncol = length(cluster_names)))
  dimnames(cur_weights) <- list(colnames(dosage_info), cluster_names)

  ##Loop through cluster_names
  for (j in 1:length(cluster_names)){
    cluster <- get(cluster_names[j])
    cur_snps <- t(dosage_info[match(c(paste(cluster$Chr, cluster$Pos, cluster$Ref, cluster$Alt, sep = ":")), rownames(dosage_info)),])
    cur_snps1 <- t(sapply(1:nrow(cur_snps),
           function(r) t(sapply(1:ncol(cur_snps), convert, r=r, cur_snps=cur_snps, cluster=cluster, vcf_info=vcf_info))))
    dimnames(cur_snps1) <- list(rownames(cur_snps), colnames(cur_snps))
    if (ncol(cur_snps1) == 1) {
      cur_weights[,j] <- apply(cur_snps1, 1,as.numeric)*as.vector(cluster$Weight)
    } else {
      cur_weights[,j] <- rowSums(t(apply(cur_snps1, 1,as.numeric)*as.vector(cluster$Weight)))
    }
  }
  if (nrow(weights) == 0){
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

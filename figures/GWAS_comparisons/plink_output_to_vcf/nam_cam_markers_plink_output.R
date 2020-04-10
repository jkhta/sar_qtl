#this file attempt to generate a PLINK formatted file from the BIMBAM formatted file that Marc sent
library(data.table)
library(dplyr)
library(gdata)

rm(list = ls())

#first reading in the marker information that marc sent and then i subsequently changed into
#a data frame
setwd("/home/jkta/projects/col_sha/nam_cam_markers_stated_creation")

#for 
nam_markers <- as.data.frame(fread("nam_cam_markers_stated_pruned_subset.csv", 
                                   sep = ",", 
                                   stringsAsFactors = FALSE))
nam_markers$geno <- sapply(strsplit(nam_markers$geno, split = "_"), function(x) x[2])

#writing the bd_geno file for the 
nam_genotypes <- nam_markers$geno

#grabbing the bd geno data
nam_pheno <- read.csv("nam_blups_combined_univariate.csv", header = TRUE, stringsAsFactors = FALSE)
nam_pheno_subset <- subset(nam_pheno, geno %in% nam_genotypes)

#for now i am going to use the data table i generated from the manual mapping working script
#going to erase the first column and the first two markers because they are useless
nam_markers_subset <- subset(nam_markers, geno %in% nam_pheno_subset$geno)
nam_genotypes_subset <- nam_markers_subset$geno

#first trying to generate the pedigree file from the marker data
#grabbing only the markers from the genotype information data frame
nam_marker_genos_only <- nam_markers_subset[, 2:ncol(nam_markers_subset)]

#doubling the number of snps by just taking the indexes of the original marker data frame,
#and then just sampling from each marker column twice
nam_marker_doubled_indices <- rep(1:ncol(nam_marker_genos_only), each = 2)
nam_marker_genos_only_doubled <- nam_marker_genos_only[, nam_marker_doubled_indices]

#i tried to use an apply function to double the columns for each marker,
#but this didn't work out, so now i'm trying to double the number of markers using a for loop
#using bind_col function in order to efficiently bind the columns together

#now trying to replace all of the non-base information with a 0
ped_missing_replacer <- function(column) {
  #trying to replace anything that isn't a base by a 0
  column[!grepl("A$|C$|G$|T$", column)] <- "0"
  # column[is.na(column)] <- "0"
  # column["ALT"] <- "0"
  # column["REF"] <- "0"
  return(column)
}

nam_marker_genos_only_doubled_replaced <- apply(nam_marker_genos_only_doubled, 2, function(x) ped_missing_replacer(x))
nam_marker_genos_only_doubled_replaced[nam_marker_genos_only_doubled_replaced == "ALT"] <- "0"

#generating the pedigree information
nam_marker_geno_names <- nam_markers_subset$geno

#
nam_ped <- data.frame(family = paste(sapply(strsplit(nam_marker_geno_names, "RV"), function(x) x[1]), "RV", sep = ""),
                      sample = nam_marker_geno_names,
                      paternal = 0,
                      maternal = 0,
                      sex = 1,
                      affection = 0)

nam_ped_with_geno <- cbind(nam_ped, nam_marker_genos_only_doubled_replaced)

#writing the ped file with all of the snp information
fwrite(nam_ped_with_geno, 
       "nam_pruned_subset.ped", 
       sep = "\t", 
       col.names = FALSE, 
       row.names = FALSE,
       nThread = 8)

#writing the map file
nam_genotypes_subset <- factor(nam_genotypes_subset, levels = nam_genotypes_subset)
nam_pheno_subset$geno <- factor(nam_pheno_subset$geno, levels = nam_genotypes_subset)
nam_pheno_subset <- nam_pheno_subset[match(nam_genotypes_subset, nam_pheno_subset$geno), ]
nam_family_information <- sapply(strsplit(as.character(nam_pheno_subset$geno), split = "RV"), function(x) x[1])
nam_family_information <- paste(nam_family_information, "RV", sep = "")
head(nam_pheno_subset)

#now need to write phenotype file; it is in the PLINK format with no header
nam_pheno_plink <- data.frame(family_id = nam_family_information,
                              individual_id = nam_pheno_subset$geno,
                              bolting_days = nam_pheno_subset$bd_geno,
                              stringsAsFactors = FALSE)

fwrite(nam_pheno_plink, file = "nam_bd_plink.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

#now trying to generate the map file
nam_marker_names <- colnames(nam_marker_genos_only)
nam_chr <- sapply(strsplit(nam_marker_names, "_"), function(x) x[1])
nam_pos <- sapply(strsplit(nam_marker_names, "_"), function(x) x[2])

nam_map <- data.frame(chr = nam_chr, 
                      marker_names = nam_marker_names, 
                      gen_dist = 0, 
                      phys_dist = nam_pos)

fwrite(nam_map, 
       "nam_pruned_subset.map", 
       sep = "\t", 
       col.names = FALSE, 
       row.names = FALSE,
       nThread = 4)

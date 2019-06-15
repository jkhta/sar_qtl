#this script will format the bd data to find QTL with the QTLiciMapping program
library(data.table)
library(plyr)
library(stringr)

#first doing the easy thing and converting the phenotype data into the proper format
#i will only be converting bolting days first to test out the program
rm(list = ls())

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/GWAS_comparisons/")

try(dir.create("QTL_icimapping"))
setwd("QTL_icimapping")

#generating the general info part
nam_ici_general <- data.frame(c(4, 2, 2, 1, 5, 7, 1))
fwrite(nam_ici_general, 
       "nam_ici_general.csv", 
       sep = ",", 
       col.names = FALSE, 
       row.names = FALSE)

#generating the chromosome information
setwd("/Users/jkhta/Desktop/nam_cam_fixing/27 - last_data/input/marc_data/")
nam_gmap <- fread("merged_NAM_lines_all_chr_11_NRP_joint_linkage_map_CSVR.csv", 
                  sep = ",",
                  header = TRUE, 
                  stringsAsFactors = FALSE)
nam_gmap <- nam_gmap[-c(1:22), ]
nam_ici_chromosome <- count(factor(nam_gmap$V2, levels=unique(nam_gmap$V2)))
nam_ici_chromosome$x <- paste("'-CH", nam_ici_chromosome$x, sep = "")

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/GWAS_comparisons/QTL_icimapping/")
fwrite(nam_ici_chromosome, 
       "nam_ici_chromosome.csv", 
       sep = ",", 
       col.names = FALSE, 
       row.names = FALSE)

#generating the linkage map information
nam_gmap_subset <- nam_gmap[, 1:3]
fwrite(nam_gmap_subset, 
       "nam_ici_linkagemap.csv", 
       sep = ",",
       col.names = FALSE, 
       row.names = FALSE)

#generating the genotype information and reading in the phenotype information
#need to match ids together because some rils have been genotyped, and some have been phenotyped 
nam_geno <- data.table(t(nam_gmap), keep.rownames = TRUE)
colnames(nam_geno) <- as.character(unlist(nam_geno[1, ]))
nam_geno <- nam_geno[-c(1:3), ]
colnames(nam_geno)[1] <- "id"

#renaming genotype column
nam_geno$id <- gsub("X", "", nam_geno$id)
geno_pop <- paste(sapply(strsplit(nam_geno$id, split = "RV"), function(x) x[1]), "RV", sep = "")
geno_num <- str_pad(sapply(strsplit(nam_geno$id, split = "RV"), function(x) x[2]), width = 3, pad = 0)
nam_geno$id <- paste(geno_pop, "_", geno_pop, geno_num, sep = "")

#relabeling letters to numbers 
nam_geno[nam_geno == "AA"] <- 2
nam_geno[nam_geno == "BB"] <- 0
nam_geno[nam_geno == "-"] <- -1
nam_geno[nam_geno == NA] <- -1

#grabbing the genotypes vector to match with the phenotypes vector 
nam_genotypes <- nam_geno$id

#writing the bd_geno file for the 
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/3_lmm_blups/output/")
nam_pheno <- fread("nam_blups_combined_univariate.csv",
                   sep = ",",
                   header = TRUE,
                   stringsAsFactors = FALSE)

nam_pheno_pop <- paste(sapply(strsplit(nam_pheno$geno, split = "RV"), function(x) x[1]), "RV", sep = "")
nam_pheno$geno <- paste(nam_pheno_pop, nam_pheno$geno, sep = "_")

#grabbing just the genos to get the lines more in sync with my data
nam_pheno_geno <- subset(nam_pheno, select = geno)

#now removing all the NA's and then subsetting the phenotype data based on the genotypes
#that have markers
nam_pheno_subset <- subset(nam_pheno, geno %in% nam_genotypes)

#matching genotypes based on phenotype data
nam_geno_subset <- subset(nam_geno, id %in% nam_pheno_subset$geno)
nam_geno_subset
nam_geno_t <- t(nam_geno_subset)
nam_geno_t <- as.data.frame(data.table(nam_geno_t, keep.rownames = TRUE))
nam_geno_t <- nam_geno_t[2:nrow(nam_geno_t), ]

#setting factors and matching up
nam_geno_subset$id <- factor(nam_geno_subset$id, levels = nam_geno_subset$id)
nam_pheno_subset$id <- factor(nam_pheno_subset$geno, levels = nam_geno_subset$id)

#generating the family information
nam_genotypes <- data.frame(id = nam_geno$id, stringsAsFactors = FALSE)

#writing the population information and getting numbers on each
nam_genotypes$pop <- sapply(strsplit(nam_genotypes$id, split = "RV"), function(x) x[1])
nam_genotypes$pop <- paste(nam_genotypes$pop, "RV", sep = "")
nam_genotypes_subset <- subset(nam_genotypes, id %in% nam_pheno_subset$geno)
nam_ici_family <- count(factor(nam_genotypes_subset$pop, levels=unique(nam_genotypes_subset$pop)))

#writing the population information
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/GWAS_comparisons/QTL_icimapping/")
fwrite(nam_ici_family, 
       "nam_ici_family.csv", 
       sep = ",", 
       col.names = FALSE, row.names = FALSE)

#matching
nam_pheno_subset <- nam_pheno_subset[match(nam_geno_subset$id, nam_pheno_subset$geno), ]

identical(nam_pheno_subset$geno, as.character(nam_geno_subset$id))

nam_pheno_subset$geno <- NULL
nam_pheno_t <- t(nam_pheno_subset)

nam_pheno_bd_geno <- as.data.frame(matrix(nam_pheno_t[1, ], ncol = ncol(nam_pheno_t)))
rownames(nam_pheno_bd_geno) <- "Bdgeno"
nam_pheno_bd_geno <- as.data.frame(data.table(nam_pheno_bd_geno, keep.rownames = TRUE))

nam_pheno_bd_gxe <- as.data.frame(matrix(nam_pheno_t[2, ], ncol = ncol(nam_pheno_t)))
rownames(nam_pheno_bd_gxe) <- "Bdgxe"
nam_pheno_bd_gxe <- as.data.frame(data.table(nam_pheno_bd_gxe, keep.rownames = TRUE))


fwrite(nam_pheno_bd_geno, "nam_ici_bd_pheno.csv", col.names = FALSE, sep = ",")
fwrite(nam_pheno_bd_gxe, "nam_ici_bd_gxe_pheno.csv", col.names = FALSE, sep = ",")

fwrite(nam_geno_t, "nam_ici_genotype_bd_geno.csv", col.names = FALSE, row.names = FALSE, sep = ",")


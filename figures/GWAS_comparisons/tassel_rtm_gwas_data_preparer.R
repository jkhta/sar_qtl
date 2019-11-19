#making output for tassel and rtm-gwas
library(data.table)
library(plyr)

rm(list = ls())

#reading in blups
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/6_push_button/data/")

nam_blups <- fread("nam_blups_combined_univariate.csv",
                   sep = ",",
                   header = TRUE,
                   stringsAsFactors = FALSE)

#changing genotype name to fit with genotype data
nam_family <- paste(sapply(strsplit(nam_blups$geno, split = "RV"), function(x) x[1]), "RV", sep = "")
nam_blups$geno <- paste(nam_family, nam_blups$geno, sep = "_")
nam_blups$family <- nam_family
family_replacement <- data.frame(family = unique(nam_family), family2 = LETTERS[1:length(unique(nam_family))], stringsAsFactors = FALSE)
nam_blups$family <- mapvalues(nam_blups$family, from = family_replacement$family, to = family_replacement$family2)
  
#grabbing only bolting time
nam_bd_geno_subset <- subset(nam_blups, select = c(geno, bd_geno, family))

#saving files
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/GWAS_comparisons/")

#saving file for tassel
colnames(nam_bd_geno_subset) <- c("<Trait>", "bd")
fwrite(nam_bd_geno_subset, "bd_geno_tassel_with_family.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

nam_bd_gxe_subset <- subset(nam_blups, select = c(geno, bd_gxe, family))
#nam_bd_gxe_subset <- subset(nam_blups, select = c(geno, bd_gxe))
colnames(nam_bd_gxe_subset) <- c("<Trait>", "bd", "factor")
fwrite(nam_bd_gxe_subset, "bd_gxe_tassel_with_family.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#saving file for rtm-gwas
colnames(nam_bd_geno_subset) <- c("Indiv", "bd")
fwrite(nam_bd_geno_subset, "bd_geno_rtm_gwas.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/GWAS_comparisons/rtm_gwas/data/")
nam_bd_gxe_subset <- subset(nam_blups, select = c(geno, bd_gxe))
colnames(nam_bd_gxe_subset) <- c("Indiv", "bd")
fwrite(nam_bd_gxe_subset, "bd_gxe_rtm_gwas.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

nam_bd_gxe_covariates <- subset(nam_blups, select = c(geno, family))
colnames(nam_bd_gxe_covariates) <- c("Indiv", "family")
fwrite(nam_bd_gxe_covariates, "bd_gxe_rtm_gwas_covar.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

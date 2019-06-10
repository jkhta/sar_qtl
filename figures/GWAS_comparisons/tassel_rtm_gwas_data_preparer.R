#making output for tassel and rtm-gwas
library(data.table)

rm(list = ls())

#reading in blups
setwd("/Users/James/Documents/GitHub/sar_qtl/figures/path_analysis_GridLMM/")

nam_blups <- fread("nam_blups_combined_univariate.csv",
                   sep = ",",
                   header = TRUE,
                   stringsAsFactors = FALSE)

#changing genotype name to fit with genotype data
nam_family <- paste(sapply(strsplit(nam_blups$geno, split = "RV"), function(x) x[1]), "RV", sep = "")
nam_blups$geno <- paste(nam_family, nam_blups$geno, sep = "_")

#grabbing only bolting time
nam_blups_subset <- subset(nam_blups, select = c(geno, bd_geno))

#saving files
setwd("/Users/James/Documents/GitHub/sar_qtl/figures/GWAS_comparisons/")

#saving file for tassel
colnames(nam_blups_subset) <- c("<Trait>", "bd")
fwrite(nam_blups_subset, "bd_geno_tassel.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#saving file for rtm-gwas
colnames(nam_blups_subset) <- c("Indiv", "bd")
fwrite(nam_blups_subset, "bd_geno_rtm_gwas.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

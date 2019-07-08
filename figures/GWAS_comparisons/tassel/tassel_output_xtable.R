#this script will format the tassel table
library(data.table)
library(xtable)

#reading in the tassel output for bd_geno
setwd("/Users/James/Documents/GitHub/sar_qtl/figures/GWAS_comparisons/tassel/")

tassel_bd_geno <- fread("tassel_bd_geno_output.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

tassel_bd_geno_subset <- subset(tassel_bd_geno, select = c(1, 2, 3, 4, 9))

tassel_bd_geno_subset <- tassel_bd_geno_subset[-c(1, nrow(tassel_bd_geno_subset)), ]

print(xtable(tassel_bd_geno_subset), include.rownames = FALSE)

#xtable for bd_gxe
tassel_bd_gxe <- fread("tassel_bd_gxe_output.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

tassel_bd_gxe_subset <- subset(tassel_bd_gxe, select = c(1, 2, 3, 4, 9))

tassel_bd_gxe_subset <- tassel_bd_gxe_subset[-c(1, nrow(tassel_bd_gxe_subset)), ]

print(xtable(tassel_bd_gxe_subset), include.rownames = FALSE)

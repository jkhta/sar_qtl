#this script will format the tassel table
library(data.table)
library(xtable)

#reading in the tassel output for bd_geno
setwd("/Users/James/Documents/GitHub/sar_qtl/figures/GWAS_comparisons/tassel/output/")

tassel_bd_geno <- fread("tassel_bd_geno_output.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

tassel_bd_geno_subset <- subset(tassel_bd_geno, select = c(1, 2, 3, 4, 9))

tassel_bd_geno_subset <- tassel_bd_geno_subset[-c(1, nrow(tassel_bd_geno_subset)), ]

print(xtable(tassel_bd_geno_subset), include.rownames = FALSE)

#xtable for bd_gxe
#tassel_bd_gxe <- fread("tassel_bd_gxe_output.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
tassel_bd_gxe <- fread("tassel_bd_gxe_output_revised_with_family.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

tassel_bd_gxe_subset <- subset(tassel_bd_gxe, select = c(1, 2, 3, 4, 9))

tassel_bd_gxe_subset <- tassel_bd_gxe_subset[-c(1, nrow(tassel_bd_gxe_subset)), ]

setwd("../figure/")
fwrite(tassel_bd_gxe_subset, "table_s5_tassel_bd_gxe_results.csv", sep = ",", row.names = FALSE)
print(xtable(tassel_bd_gxe_subset), include.rownames = FALSE)

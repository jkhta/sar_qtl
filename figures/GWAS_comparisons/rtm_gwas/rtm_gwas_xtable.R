#this script will format the tassel table
library(data.table)
library(xtable)

rm(list = ls())

#reading in the tassel output for bd_geno
setwd("/Users/James/Documents/GitHub/sar_qtl/figures/GWAS_comparisons/rtm_gwas/")

rtm_bd_geno <- fread("rtm-gwas-assoc.out.190610_103359485.aov1", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rtm_bd_geno_subset <- rtm_bd_geno[!grepl("EV", rtm_bd_geno$Source), ]

print(xtable(rtm_bd_geno_subset), include.rownames = FALSE)

rtm_bd_gxe <- fread("rtm-gwas-assoc.out.190708_135201489.aov1", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rtm_bd_gxe_subset <- rtm_bd_gxe[!grepl("EV", rtm_bd_gxe$Source), ]

print(xtable(rtm_bd_gxe_subset), include.rownames = FALSE)

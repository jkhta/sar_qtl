#this script will generate a table with the accessions used in the experiments
library(data.table)
library(xtable)

rm(list = ls())

setwd("/Users/James/Documents/GitHub/sar_qtl/1_phenotype_standardization/output/")

nam_data <- fread("nam_cam_data_combined_tformed_std.csv",
                  sep = ",",
                  header = TRUE,
                  stringsAsFactors = FALSE)

head(nam_data)

nam_data_acc <- subset(nam_data, !grepl("RV", geno))
nam_data_acc$geno2 <- as.factor(nam_data_acc$geno2)

nam_acc <- data.frame(Accession = levels(nam_data_acc$geno2)[1:27],
                      Accession = levels(nam_data_acc$geno2)[28:54],
                      Accession = levels(nam_data_acc$geno2)[55:81],
                      Accession = c(levels(nam_data_acc$geno2)[82:107], "NA"))

setwd("/Users/James/Documents/GitHub/sar_qtl/figures/accession_table/figure/")
fwrite(nam_acc, "table_s1_accessions.csv", sep = ",", row.names = FALSE)
print(xtable(nam_acc), include.rownames = FALSE, NA.string = "")

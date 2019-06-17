#this script will generate a LaTeX output of different tables
library(data.table)
library(xtable)

rm(list = ls())

#table output for trait effects and heritabilities
setwd("/Users/James/Documents/GitHub/sar_qtl/figures/h2_table/")

h2_table <- fread("trait_effects_and_h2_no_ci.csv", 
                  sep = ",", 
                  header = TRUE, 
                  stringsAsFactors = FALSE)

#changing column names and removing underscores
colnames(h2_table)[1] <- NA
rownames(h2_table) <- NULL
colnames(h2_table) <- c("Trait", "Shelf Fixef", "Treatment Fixef", "Geno Var", "GxE Var", "Res", "H2", "GxE PVE")

#generating xtable
print(xtable(h2_table, caption = "Table 1", digits = 3), include.rownames=FALSE)

#reading in qtl found for genotype random effects and GxE random effects
setwd("/Users/James/Documents/GitHub/sar_qtl/figures/qtl_table/")
qtl_table_geno <- lapply(list.files(pattern = "geno"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE))
qtl_table_geno_traits <- sapply(strsplit(list.files(pattern = "geno"), "_geno"), function(x) x[1])

#adding column with trait names for each table
for (i in 1:length(qtl_table_geno_traits)) {
    qtl_table_geno[[i]]$trait <- qtl_table_geno_traits[i]
}

#combining tables together
qtl_table_geno <- rbindlist(qtl_table_geno)

#changing column names and removing underscores
qtl_table_geno$avg_snp_betas_sq <- NULL
colnames(qtl_table_geno) <- c("QTL", "SNP PVE", "-log10p", "Left Bound", "Right Bound", "Trait")

qtl_table_geno <- subset(qtl_table_geno, select = c("Trait", "QTL", "SNP PVE", "-log10p", "Left Bound", "Right Bound"))

#generating xtable
xtable(qtl_table_geno, caption = "Table 2")

#reading in qtl found for GxE random effects
qtl_table_gxe <- lapply(list.files(pattern = "gxe"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE))
qtl_table_gxe_traits <- sapply(strsplit(list.files(pattern = "gxe"), "_gxe"), function(x) x[1])

#adding column with trait names for each table
for (i in 1:length(qtl_table_gxe_traits)) {
    qtl_table_gxe[[i]]$trait <- qtl_table_gxe_traits[i]
}

#combining tables together
qtl_table_gxe <- rbindlist(qtl_table_gxe)

#changing column names and removing underscores
qtl_table_gxe$avg_snp_betas_sq <- NULL
colnames(qtl_table_gxe) <- c("QTL", "SNP PVE", "-log10p", "Left Bound", "Right Bound", "Trait")

#reordering columns so trait is first
qtl_table_gxe <- subset(qtl_table_gxe, select = c("Trait", "QTL", "SNP PVE", "-log10p", "Left Bound", "Right Bound"))

#generating xtable
print(xtable(qtl_table_gxe, caption = "Table 2"), include.rownames = FALSE)

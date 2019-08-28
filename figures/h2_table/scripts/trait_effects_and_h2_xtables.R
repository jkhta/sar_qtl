#this script will generate a LaTeX output of different tables
library(data.table)
library(xtable)

rm(list = ls())

#table output for trait effects and heritabilities
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/h2_table/input/")

h2_table <- fread("trait_effects_and_h2_no_ci.csv", 
                  sep = ",", 
                  header = TRUE, 
                  stringsAsFactors = FALSE)

acc_h2_table <- fread("acc_trait_effects_and_h2_no_ci.csv",
                      sep = ",", 
                      header = TRUE, 
                      stringsAsFactors = FALSE)

#grabbing only the gxe heritability
acc_h2_table_subset <- subset(acc_h2_table, select = c(trait, acc_gxe_h2))

#nam stats with acc stats
h2_table_merge <- merge(h2_table, acc_h2_table_subset, by = "trait")

#changing column names and removing underscores
rownames(h2_table_merge) <- NULL
colnames(h2_table_merge) <- c("Trait", "Shelf", "Intercept", "Treatment", "Geno Var", "GxE Var", "Residual Var", "H2", "GxE PVE", "Acc GxE PVE")

#multiplying H2 and GxE PVE by 100 to represent percentages
h2_table_merge$H2 <- h2_table_merge$H2 * 100
h2_table_merge$`GxE PVE` <- h2_table_merge$`GxE PVE` * 100 
h2_table_merge$`Acc GxE PVE` <- h2_table_merge$`Acc GxE PVE` * 100 

#moreving underscores
h2_table_merge$Trait <- gsub("_", "", h2_table_merge$Trait)

#changing trait names
h2_table_merge$Trait <- revalue(h2_table_merge$Trait, c("bd" = "BD_SAR", "h3h1" = "IG_SAR", "idry" = "IB_SAR", "rdry" = "RB_SAR"))

#removing columns that are already included in Table 1
h2_table_merge$Intercept <- NULL
h2_table_merge$Treatment <- NULL

#generating xtable
print(xtable(h2_table_merge, label = ("S2_Table"), digits = 2), include.rownames=FALSE)

#table with sun intercept, treatment fixed effect, and coefficient of variation
h2_table_new <- data.frame(Trait = h2_table$trait,
                           Intercept = h2_table$int_fixef,
                           Treatment = h2_table$treatment_fixef,
                           CV = sqrt(h2_table$gxe_var)/h2_table$treatment_fixef,
                           stringsAsFactors = FALSE)
h2_table_new$Trait <- gsub("_", "", h2_table_new$Trait)

#changing trait names
h2_table_new$Trait <- revalue(h2_table_new$Trait, c("bd" = "BD_SAR", "h3h1" = "IG_SAR", "idry" = "IB_SAR", "rdry" = "RB_SAR"))

print(xtable(h2_table_new, label = ("Table 1"), caption = "Averages of the intercept and treatment fixed effects, and the coefficient of variation for plasticity (CV) for each trait. bd, bolting time; h3h1, inflorescence growth over 2 weeks; rdry, dry rosette biomass; idry, dry inflorescence biomass.", digits = 2), include.rownames=FALSE)


#reading in qtl found for genotype random effects and GxE random effects
setwd("/Users/James/Documents/GitHub/sar_qtl/figures/qtl_table/data/")
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

#multiplying SNP PVE by 100 to represent a percentage
qtl_table_geno$`SNP PVE` <- qtl_table_geno$`SNP PVE` * 100

#generating xtable
print(xtable(qtl_table_geno, caption = "Table 2"), include.rownames=FALSE)

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

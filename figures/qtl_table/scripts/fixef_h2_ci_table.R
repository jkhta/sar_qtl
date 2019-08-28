#this script will generate a table with the credible intervals for the different statistics
library(data.table)
library(xtable)
library(plyr)

rm(list = ls())

#reading in the credible interval table
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/h2_table/input/")

fixef_h2_ci_table <- fread("trait_effects_and_h2_ci.csv",
                           sep = ",",
                           header = TRUE,
                           stringsAsFactors = FALSE, 
                           data.table = FALSE)

for (i in 1:ncol(fixef_h2_ci_table)) {
  fixef_h2_ci_table_col <- unlist(fixef_h2_ci_table[, i])
  fixef_h2_ci_table_col <- gsub("to", "-", fixef_h2_ci_table_col)
  fixef_h2_ci_table[, i] <- fixef_h2_ci_table_col
}

#choosing only some columns
fixef_h2_ci_table_subset <- subset(fixef_h2_ci_table, select = c(1, 2, 3, 4, 8, 9))
colnames(fixef_h2_ci_table_subset) <- c("Trait", "Shelf", "Intercept","Treatment", "Geno H2", "GxE PVE")

#changing trait names
fixef_h2_ci_table_subset$Trait <- mapvalues(fixef_h2_ci_table_subset$Trait, from = c("bd_ci", "h3_h1_ci", "r_dry_ci", "i_dry_ci"), to = c("BD_SAR", "IG_SAR", "RB_SAR", "IB_SAR"))

print(xtable(fixef_h2_ci_table_subset), include.rownames=FALSE)

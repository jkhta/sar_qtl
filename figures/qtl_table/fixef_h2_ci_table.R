library(data.table)
library(xtable)

rm(list = ls())

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/h2_table/")

fixef_h2_ci_table <- fread("trait_effects_and_h2_ci.csv",
                           sep = ",",
                           header = TRUE,
                           stringsAsFactors = FALSE)

fixef_h2_ci_table_subset <- subset(fixef_h2_ci_table, select = c(1, 2, 3, 7, 8))
colnames(fixef_h2_ci_table_subset) <- c("Trait", "Shelf Fixef", "Treatment Fixef", "Geno H2", "GxE PVE")

print(xtable(fixef_h2_ci_table_subset), include.rownames=FALSE)

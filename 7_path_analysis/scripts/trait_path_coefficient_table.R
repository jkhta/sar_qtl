#this script will generate tables for the trait correlations in sun and shade conditions,
#estimated from the path models
library(data.table)
library(xtable)

rm(list = ls())

setwd("/Users/James/Documents/GitHub/sar_qtl/7_path_analysis/output/")

sun_trait_eff <- rbindlist(lapply(list.files(pattern = "trait_eff_alt_sun.csv"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE)))

shade_trait_eff <- rbindlist(lapply(list.files(pattern = "trait_eff_alt_shade.csv"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE)))

sun_trait_eff_reorder <- subset(sun_trait_eff, select = c(pop, grep("~", colnames(sun_trait_eff))))
shade_trait_eff_reorder <- subset(shade_trait_eff, select = c(pop, grep("~", colnames(shade_trait_eff))))

print(xtable(sun_trait_eff_reorder), include.rownames = FALSE, NA.string = "NA")
print(xtable(shade_trait_eff_reorder), include.rownames = FALSE, NA.string = "NA")

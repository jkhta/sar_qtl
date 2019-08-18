#this script will generate tables for the trait correlations in sun and shade conditions,
#estimated from the path models
library(data.table)
library(xtable)
library(plyr)

rm(list = ls())

setwd("/Users/James/Documents/GitHub/sar_qtl/7_path_analysis/data/")

sun_trait_eff <- rbindlist(lapply(list.files(pattern = "trait_eff_blups_sun.csv"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE)))

shade_trait_eff <- rbindlist(lapply(list.files(pattern = "trait_eff_blups_shade.csv"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE)))

sun_trait_eff_reorder <- subset(sun_trait_eff, select = c(pop, grep("~", colnames(sun_trait_eff))))
shade_trait_eff_reorder <- subset(shade_trait_eff, select = c(pop, grep("~", colnames(shade_trait_eff))))

colnames(sun_trait_eff_reorder)[1] <- "Population"
colnames(shade_trait_eff_reorder)[1] <- "Population"

sun_trait_eff_reorder$pop <- revalue(sun_trait_eff_reorder$pop, c("blh" = "Blh-1", "bur" = "Bur-0", "cvi" = "Cvi-0", "ita" = "Ita-0", "jea" = "Jea", "oy" = "Oy-0", "sha" = "Sha"))
shade_trait_eff_reorder$pop <- revalue(shade_trait_eff_reorder$pop, c("blh" = "Blh-1", "bur" = "Bur-0", "cvi" = "Cvi-0", "ita" = "Ita-0", "jea" = "Jea", "oy" = "Oy-0", "sha" = "Sha"))

print(xtable(sun_trait_eff_reorder), include.rownames = FALSE, NA.string = "NA")
print(xtable(shade_trait_eff_reorder), include.rownames = FALSE, NA.string = "NA")

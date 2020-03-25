#this script will generate tables for the trait correlations in sun and shade conditions,
#estimated from the path models
library(data.table)
library(xtable)
library(plyr)

rm(list = ls())

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/7_path_analysis/data/")

sun_trait_eff <- rbindlist(lapply(list.files(pattern = "trait_eff_blups_sun.csv"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE)))

shade_trait_eff <- rbindlist(lapply(list.files(pattern = "trait_eff_blups_shade.csv"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE)))

sun_trait_eff_reorder <- subset(sun_trait_eff, select = c(pop, grep("~", colnames(sun_trait_eff))))
shade_trait_eff_reorder <- subset(shade_trait_eff, select = c(pop, grep("~", colnames(shade_trait_eff))))

colnames(sun_trait_eff_reorder)[1] <- "Population"
colnames(shade_trait_eff_reorder)[1] <- "Population"

sun_trait_eff_reorder$Population <- revalue(sun_trait_eff_reorder$Population, c("blh" = "Blh-1", "bur" = "Bur-0", "cvi" = "Cvi-0", "ita" = "Ita-0", "jea" = "Jea", "oy" = "Oy-0", "sha" = "Sha"))
shade_trait_eff_reorder$Population <- revalue(shade_trait_eff_reorder$Population, c("blh" = "Blh-1", "bur" = "Bur-0", "cvi" = "Cvi-0", "ita" = "Ita-0", "jea" = "Jea", "oy" = "Oy-0", "sha" = "Sha"))

colnames(sun_trait_eff_reorder) <- gsub("rdry", "RB", colnames(sun_trait_eff_reorder))
colnames(sun_trait_eff_reorder) <- gsub("bd", "BD", colnames(sun_trait_eff_reorder))
colnames(sun_trait_eff_reorder) <- gsub("h3h1", "IG", colnames(sun_trait_eff_reorder))
colnames(sun_trait_eff_reorder) <- gsub("idry", "IB", colnames(sun_trait_eff_reorder))

colnames(shade_trait_eff_reorder) <- gsub("rdry", "RB", colnames(shade_trait_eff_reorder))
colnames(shade_trait_eff_reorder) <- gsub("bd", "BD", colnames(shade_trait_eff_reorder))
colnames(shade_trait_eff_reorder) <- gsub("h3h1", "IG", colnames(shade_trait_eff_reorder))
colnames(shade_trait_eff_reorder) <- gsub("idry", "IB", colnames(shade_trait_eff_reorder))

setwd("../figure/")
fwrite(sun_trait_eff_reorder, "table_s7_sun_trait_effects.csv", sep = ",", row.names = FALSE)
fwrite(shade_trait_eff_reorder, "table_s8_shade_trait_effects.csv", sep = ",", row.names = FALSE)

print(xtable(sun_trait_eff_reorder), include.rownames = FALSE, NA.string = "NA")
print(xtable(shade_trait_eff_reorder), include.rownames = FALSE, NA.string = "NA")

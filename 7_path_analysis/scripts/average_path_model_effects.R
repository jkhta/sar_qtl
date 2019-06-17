#this script will average the qtl effects and traits across populations
library(data.table)

rm(list = ls())

setwd("/Users/James/Documents/GitHub/sar_qtl/7_path_analysis/output/")

pop_qtl_eff_list_sun <- as.data.frame(rbindlist(lapply(list.files(pattern = "qtl_eff_alt_sun.csv"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE))))
pop_qtl_eff_list_sun[is.na(pop_qtl_eff_list_sun)] <- 0
pop_qtl_eff_list_sun_avg <- apply(pop_qtl_eff_list_sun[,grepl("~", colnames(pop_qtl_eff_list_sun))], 2, function(x) mean(x))

pop_qtl_eff_list_shade <- as.data.frame(rbindlist(lapply(list.files(pattern = "qtl_eff_alt_shade.csv"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE))))
pop_qtl_eff_list_shade[is.na(pop_qtl_eff_list_shade)] <- 0
pop_qtl_eff_list_shade_avg <- apply(pop_qtl_eff_list_shade[,grepl("~", colnames(pop_qtl_eff_list_shade))], 2, function(x) mean(x))

pop_trait_eff_list_sun <- as.data.frame(rbindlist(lapply(list.files(pattern = "trait_eff_alt_sun.csv"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE))))
pop_trait_eff_list_sun[is.na(pop_trait_eff_list_sun)] <- 0
pop_trait_eff_list_sun_avg <- apply(pop_trait_eff_list_sun[,grepl("~", colnames(pop_trait_eff_list_sun))], 2, function(x) mean(x))

pop_trait_eff_list_shade <- as.data.frame(rbindlist(lapply(list.files(pattern = "trait_eff_alt_shade.csv"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE))))
pop_trait_eff_list_shade[is.na(pop_trait_eff_list_shade)] <- 0
pop_trait_eff_list_shade_avg <- apply(pop_trait_eff_list_shade[,grepl("~", colnames(pop_trait_eff_list_shade))], 2, function(x) mean(x))


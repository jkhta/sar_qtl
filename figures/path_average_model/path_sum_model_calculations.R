library(data.table)

rm(list = ls())

setwd("/Users/James/Documents/GitHub/sar_qtl/figures/path_average_model/")

significant_paths <- fread("sun_shade_significant_paths.csv", 
                           sep = ",",
                           header = TRUE,
                           stringsAsFactors = FALSE)

sun_sig_paths <- subset(significant_paths, Condition == "Sun", select = grepl("~", colnames(significant_paths)))
sun_sig_paths_sum <- apply(sun_sig_paths, 2, function(x) sum(x))

shade_sig_paths <- subset(significant_paths, Condition == "Shade", select = grepl("~", colnames(significant_paths)))
shade_sig_paths_sum <- apply(shade_sig_paths, 2, function(x) sum(x))

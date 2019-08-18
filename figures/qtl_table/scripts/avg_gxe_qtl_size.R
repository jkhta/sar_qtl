#this script will generate estimates of the average QTL in Mb

library(data.table)

rm(list = ls())

setwd("/Users/James/Documents/GitHub/sar_qtl/figures/qtl_table/data/")

sar_qtl <- fread("sar_gxe_qtl_info_with_anno.csv",
                 sep = ",",
                 header = TRUE,
                 stringsAsFactors = FALSE)

qtl_window <- data.frame(left = as.numeric(sapply(strsplit(sar_qtl$`Left Bound`, split = "_"), function(x) x[length(x)])),
                         right = as.numeric(sapply(strsplit(sar_qtl$`Right Bound`, split = "_"), function(x) x[length(x)])))

qtl_window$size <- with(qtl_window, (right - left)/1000000)

qtl_window_avg_size <- mean(qtl_window$size)
                         

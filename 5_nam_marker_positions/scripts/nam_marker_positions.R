library(data.table)

rm(list = ls())

#this sript will grab the marker info for the latest marker set
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/4_geno_prob_array/input/")

#reading in the marker/phenotype data
marc_final_markers <- fread("merged_NAM_lines_all_chr_11_NRP_joint_linkage_map_CSVR.csv",
                            sep = ",",
                            header = TRUE,
                            stringsAsFactors = FALSE)

#outputting the marker information (first three columns from rows 23 onwards)
marc_final_marker_info <- marc_final_markers[23:nrow(marc_final_markers), 1:3]
colnames(marc_final_marker_info) <- c("snp", "chr", "cM")

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/5_nam_marker_positions/output/")

#writing the genotype file
fwrite(marc_final_marker_info, 
       "nam_marker_info_final.csv", 
       sep = ",", 
       row.names = FALSE, 
       col.names = TRUE)

#path analysis
#need to grab the snps m_4_407208 and m_5_3799350
library(data.table)

rm(list = ls())
#reading in the genotype data
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/6_push_button/")

geno_data <- readRDS("nam_rqtl_geno_prob_array_comp_11_all_rils_NEW.RDS")

#grabbing the qtl markers
geno_data_markers <- geno_data[c("4_407208", "5_3799350")]

#now need to 
geno_data_markers_col <- lapply(geno_data_markers, function(x) subset(x, select = 1))

geno_data_markers_both <- lapply(geno_data_markers_col, function(x) data.frame(A = x[,1], B = 1 - x[,1]))

geno_data_markers_both_df <- do.call(cbind, geno_data_markers_both)

colnames(geno_data_markers_both_df) <- paste("m", colnames(geno_data_markers_both_df), sep = "_")

geno_data_qtl_combos <- with(geno_data_markers_both_df, data.frame(geno = rownames(geno_data_markers_col[[1]]),
                                                                   A4 = m_4_407208.A,
                                                                   B4 = m_4_407208.B,
                                                                   A5 = m_5_3799350.A,
                                                                   B5 = m_5_3799350.B,
                                                                   AA = m_4_407208.A * m_5_3799350.A, 
                                                                   AB = m_4_407208.A * m_5_3799350.B,
                                                                   BA = m_4_407208.B * m_5_3799350.A,
                                                                   BB = m_4_407208.B * m_5_3799350.B))


#reading in phenotype data (plant data)
setwd("/Users/jkhta/Desktop/nam_cam_fixing/10.5 - nam_cam_other_rep_data_cleaning/output/")
pheno_data <- fread("nam_cam_data_combined_tformed_std_FINAL.csv",
                    sep = ",", header = TRUE, stringsAsFactors = FALSE)

#subsetting phenotype data that is not the accessions
unique(pheno_data$cross)
pheno_data_subset <- subset(pheno_data, cross != "accession")
geno_family <- paste(sapply(str_split(pheno_data_subset$geno, pattern = "RV"), function(x) x[1]), "RV", sep = "")
genotypes <- paste(geno_family, pheno_data_subset$geno, sep = "_")
pheno_data_subset$geno <- genotypes

geno_pheno_merge <- merge(pheno_data_subset, geno_data_qtl_combos, by = "geno")
setwd("/Users/jkhta/Desktop/nam_cam_fixing/34 - pop_path_analysis/output/")
fwrite(geno_pheno_merge, 
       "nam_pheno_qtl_combo_merge.csv",
       sep = ",", 
       row.names = FALSE, 
       na = "NA")

geno_pheno_merge_by_pop <- lapply(unique(geno_pheno_merge$cross), function(x) subset(geno_pheno_merge, cross == x))

for (i in 1:length(geno_pheno_merge_by_pop)) {
  file_name <- paste(unique(geno_pheno_merge_by_pop[[i]]$cross), "pheno_qtl_combo_merge.csv", sep = "_")
  fwrite(geno_pheno_merge_by_pop[[i]], 
         file_name, 
         sep = ",", 
         row.names = FALSE,
         na = "NA")
}
                             
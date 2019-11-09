#path analysis
#need to grab the snps m_4_41028, m_5_3799350, 4_8938713, and 5_25961748
library(data.table)
library(plyr)

rm(list = ls())
#reading in the genotype data
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/6_push_button/")

geno_data <- readRDS("nam_rqtl_geno_prob_array_comp_11_all_rils_NEW.RDS")

#grabbing the qtl markers
#i chose the markers that were significant, and are more in the middle of the QTL regions 
geno_data_markers <- geno_data[c("4_41028", "5_3799350", "4_8938713", "5_25961748")]

#now need to 
geno_data_markers_col <- lapply(geno_data_markers, function(x) subset(x, select = 1))

geno_data_markers_both <- lapply(geno_data_markers_col, function(x) data.frame(A = x[,1], B = 1 - x[,1]))

geno_data_markers_both_df <- do.call(cbind, geno_data_markers_both)

colnames(geno_data_markers_both_df) <- paste("m", colnames(geno_data_markers_both_df), sep = "_")

geno_data_qtl_combos <- with(geno_data_markers_both_df, data.frame(geno = rownames(geno_data_markers_col[[1]]),
                                                                   A4 = m_4_41028.A,
                                                                   B4 = m_4_41028.B,
                                                                   A5 = m_5_3799350.A,
                                                                   B5 = m_5_3799350.B,
                                                                   A42 = m_4_8938713.A,
                                                                   B42 = m_4_8938713.B,
                                                                   A52 = m_5_25961748.A,
                                                                   B52 = m_5_25961748.B,
                                                                   AA = m_4_41028.A * m_5_3799350.A, 
                                                                   AB = m_4_41028.A * m_5_3799350.B,
                                                                   BA = m_4_41028.B * m_5_3799350.A,
                                                                   BB = m_4_41028.B * m_5_3799350.B))


#reading in phenotype data (plant data)
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/3_lmm_blups/output/")
pheno_data <- fread("pop_sun_shade_blups.csv",
                    sep = ",", 
                    header = TRUE, 
                    stringsAsFactors = FALSE)

#subsetting phenotype data that is not the accessions
pheno_data_subset <- subset(pheno_data, grepl("RV", pheno_data$geno))
geno_family <- paste(sapply(strsplit(pheno_data_subset$geno, split = "RV"), function(x) x[1]), "RV", sep = "")
genotypes <- paste(geno_family, pheno_data_subset$geno, sep = "_")
pheno_data_subset$geno <- genotypes

geno_pheno_merge <- merge(pheno_data_subset, geno_data_qtl_combos, by = "geno")
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/7_path_analysis/data/")
fwrite(geno_pheno_merge, 
       "nam_blups_qtl_combo_merge.csv",
       sep = ",", 
       row.names = FALSE, 
       na = "NA")

geno_pheno_merge_by_pop <- lapply(paste(unique(geno_family), unique(geno_family), sep = "_"), function(x) subset(geno_pheno_merge, grepl(x, geno)))
geno_family_unique <- unique(geno_family)
geno_family_pop_names <- revalue(geno_family_unique, c("21RV" = "blh", "20RV" = "bur", "8RV" = "cvi", "29RV" = "ita", "28RV" = "jea", "27RV" = "oy", "13RV" = "sha"))
for (i in 1:length(geno_pheno_merge_by_pop)) {
  file_name <- paste(unique(geno_family_pop_names)[i], "blups_qtl_combo_merge.csv", sep = "_")
  fwrite(geno_pheno_merge_by_pop[[i]], 
         file_name, 
         sep = ",", 
         row.names = FALSE,
         na = "NA")
}

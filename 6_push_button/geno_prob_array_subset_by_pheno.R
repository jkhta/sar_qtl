#reading in the genotype array for all genotypes, then going to subset based
#on my phenotypes
all_ril_geno_probs <- readRDS("nam_rqtl_geno_prob_array_comp_11_all_rils_NEW.RDS")

#reading in my phenotype data
nam_blups <- fread("nam_blups_combined_univariate.csv",
                   sep = ",",
                   header = TRUE,
                   stringsAsFactors = FALSE)

#removing missing data
nam_blups <- nam_blups[complete.cases(nam_blups), ]

#matching up the genotype column of the phenotype data with the genotype probability array
nam_blups$geno <- with(nam_blups, paste(paste(sapply(strsplit(geno, split = "RV"), function(x) x[1]), "RV", sep = ""), geno, sep = "_"))

match(nam_blups$geno, rownames(all_ril_geno_probs[[1]]))
match(rownames(all_ril_geno_probs[[1]]), nam_blups$geno)

#now need to subset the marker array with only genotypes that have phenotypes
ril_geno_probs_subset <- mclapply(all_ril_geno_probs, function(x) subset(x, rownames(x) %in% nam_blups$geno), mc.cores = detectCores())

#need to make sure these match up
nam_blups_subset <- nam_blups[match(rownames(ril_geno_probs_subset[[1]]), nam_blups$geno), ]

#checking to make sure the phenotype line and genotype line order is the same
identical(nam_blups_subset$geno, rownames(ril_geno_probs_subset[[1]]))

#saving the file
saveRDS(ril_geno_probs_subset, "nam_rqtl_geno_prob_array_final_pheno_subset.RDS")

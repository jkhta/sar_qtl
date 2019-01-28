#reading in the pruned genotypes, and will subset the complete genotypes later by the markers remaining in the pruned set
Genotypes_pruned <- readRDS("nam_rqtl_geno_prob_array_final_0.99_no_kinship.rds")

#READIN IN THE MARKERS FOR genotypes only found within the blups data set
#reading in the unpruned data set and grabbing the markers that were not pruned
Genotypes_unpruned <- readRDS("nam_rqtl_geno_prob_array_final_pheno_subset.RDS")

#reading in all of the phenotype data; need to read this in first to subset genotypes and kinship matrix
nam_pheno <- fread("nam_blups_combined_univariate.csv", 
                   sep = ",", 
                   header = TRUE, 
                   stringsAsFactors = FALSE)

#trying to match up the phenotype data with the kinship matrix 
nam_pheno_pop <- paste(sapply(strsplit(nam_pheno$geno, split = "RV"), function(x) x[1]), "RV", sep = "")
nam_pheno$geno <- paste(nam_pheno_pop, nam_pheno$geno, sep = "_")
nam_pheno <- subset(nam_pheno, geno %in% rownames(Genotypes_unpruned[[1]]))

#now need to subset the genotypes 
Genotypes_unpruned <- mclapply(Genotypes_unpruned, function(x) x[match(nam_pheno$geno, rownames(x)), ])

#grabbing only the markers that are found in the pruned subset
Genotypes <- Genotypes_unpruned[names(Genotypes_unpruned) %in% names(Genotypes_pruned)]

#Genotypes_subset_no <- Genotypes_subset_no[order(Genotypes_subset_no)]
names(Genotypes) <- paste("m", names(Genotypes), sep = "_")
marker_names <- names(Genotypes)

#making a new column for the population type
nam_pheno$pop <- sapply(strsplit(nam_pheno$geno, split = "_"), function(x) x[1])
nam_pheno$pop_pop <- paste(nam_pheno$pop, nam_pheno$pop, sep = "_")

#function to calculate the correlations between current marker and all other markers
#using the genotype probabilities
marker_cor <- function(marker_matrix, marker_col) {
  #pulling out the testing marker
  testing_marker <- as.vector(marker_matrix[, marker_col])
  
  #calculating marker correlations
  marker_cors <- unlist(mclapply(1:ncol(marker_matrix), function(y) cor(testing_marker, as.vector(marker_matrix[, y])), mc.cores = detectCores()))
  marker_cors <- ifelse(marker_cors > cor_threshold, 1, 0)
  
  return(marker_cors)
}

#this function will pull out data for an individual population
pheno_and_geno_pop_data_puller <- function(pheno_data, marker_prob_array, pop_name) {
  #grabbing the population specific pheno data 
  pop_pheno_data <- subset(pheno_data, grepl(pop_name, geno))
  pop_pheno_data$geno_2 <- pop_pheno_data$geno
  
  #subsetting the marker array for just the population
  pop_marker_array <- mclapply(marker_prob_array, function(x) x[grepl(pop_name, rownames(x)), ])
  
  #now making a marker prob matrix for the single biparental pop
  pop_geno_list <- mclapply(pop_marker_array, function(x) as.data.frame(x[, 1]))
  pop_geno_mat <- as.matrix(do.call(cbind, pop_geno_list))
  rownames(pop_geno_mat) <- rownames(pop_marker_array[[1]])
  colnames(pop_geno_mat) <- names(marker_prob_array)
  
  #now need to make sure the pheno and the geno data are matching each other in terms of genotypes
  pop_pheno_subset <- subset(pop_pheno_data, geno %in% rownames(pop_geno_mat))
  pop_geno_mat_subset <- subset(pop_geno_mat, rownames(pop_geno_mat) %in% pop_pheno_subset$geno)
  
  #need to generate a proximal matrix for each population
  pop_marker_window_proximal <- centimorgan_proximal_matrix_generator(pop_geno_mat_subset, nam_marc_marker_info, cm_threshold)
  pop_marker_cor <- lapply(1:ncol(pop_geno_mat_subset), function(x) marker_cor(pop_geno_mat_subset, x))
  pop_marker_cors <- as.matrix(rbindlist(lapply(pop_marker_cor, function(x) as.data.frame(t(as.matrix(x))))))
  
  #returning the pheno and geno data
  return(list(pop_pheno = pop_pheno_data,
              pop_geno = pop_geno_mat,
              marker_cors = pop_marker_cors,
              cM_proximal = pop_marker_window_proximal))
}

#need to merge file with genetic distances to plot genetic distances
nam_marc_marker_info <- fread("nam_marker_info_final.csv",
                              sep = ",", 
                              header = TRUE, 
                              stringsAsFactors = FALSE)
colnames(nam_marc_marker_info)[1] <- "snp"
nam_marc_marker_info$snp <- paste("m", nam_marc_marker_info$snp, sep = "_")

#this function will return a row with 0's and 1's for removing markers based on a 
#cM window
proximal_matrix_row <- function(test_snp, marker_info, window_size) {
  #grabbing the chromosome of the test snp so that i can subset the marker info
  test_snp_chr <- sapply(strsplit(test_snp, split = "_"), function(x) x[2])
  
  #subsetting the markers from the same chr as the test snp
  marker_info_same_snp_chr <- subset(marker_info, chr == test_snp_chr)
  marker_info_other_chr_proximal <- rep(0, nrow(subset(marker_info, chr != test_snp_chr)))
  
  #getting the distance of the test snp in order to calculate snp distances
  test_snp_dist <- subset(marker_info_same_snp_chr, snp == test_snp)$cM
  marker_info_same_snp_chr$cM_from_test <- abs(test_snp_dist - marker_info_same_snp_chr$cM)
  
  #now need to find which snps are in the same interval
  markers_within_bounds <- subset(marker_info_same_snp_chr, cM_from_test < window_size)$snp
  marker_pos <- which(marker_info$snp %in% markers_within_bounds)
  snps_within_int <- rep(0, nrow(marker_info))
  snps_within_int[marker_pos] <- 1
  
  #making the row a data frame so that i can combine the different rows of the 
  #proximal matrix
  snps_within_int_row <- as.data.frame(t(as.data.frame(snps_within_int)))
  
  #returning the row of the proximal matrix
  return(snps_within_int_row)
}

#this function will generate a proximal matrix and remove markers within a certain interval  
centimorgan_proximal_matrix_generator <- function(pop_geno_mat, marker_info, window_size) {
  #this will generate a marker info subset based on the population markers
  marker_info_subset <- subset(marker_info, snp %in% colnames(pop_geno_mat))
  
  #now need to loop through each marker and generate a proximal row for that marker
  #using an mclapply in order to generate the proximal matrix
  gen_dist_proximal <- as.matrix(rbindlist(lapply(1:ncol(pop_geno_mat), function(x) proximal_matrix_row(colnames(pop_geno_mat)[x], marker_info_subset, window_size = cm_threshold))))
  
  return(gen_dist_proximal)
}


#this prepares all the data for each individual population
all_pop_data <- lapply(unique(nam_pheno$pop_pop), function(x) pheno_and_geno_pop_data_puller(pheno_data = nam_pheno,
                                                                                             marker_prob_array = Genotypes,
                                                                                             pop_name = x))
names(all_pop_data) <- unique(nam_pheno$pop_pop)

saveRDS(all_pop_data, file = "nam_all_traits_ind_pop_pheno_geno_proximal_final.RDS")
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

#need to merge file with genetic distances to plot genetic distances
colnames(nam_marc_marker_info)[1] <- "snp"
nam_marc_marker_info$snp <- paste("m", nam_marc_marker_info$snp, sep = "_")

#NEED TO SAVE THIS AS AN R OBJECT SO I DON'T HAVE TO KEEP RUNNING IT OVER AND OVER AGAIN

#grabbing the individual population data
pheno_name <- colnames(all_pop_data$`21RV_21RV`$pop_pheno)
pheno_name <- pheno_name[grepl(phenotype_list, pheno_name)]

#grabbing the markers to generate a marker map
pop_markers <- all_pop_data$`21RV_21RV`$pop_geno

#generating a marker map: it can be physical or genetic depending on the choices
marker_map <- data.frame(snp = colnames(pop_markers), stringsAsFactors = FALSE)


#generating markers depending on physical or genetic distance
if (pos_type == "phys_dist") {
  marker_map$Chr <- as.character(sapply(strsplit(marker_map$snp, split = "_"), function(x) x[2]))
  marker_map$pos <- as.numeric(sapply(strsplit(marker_map$snp, split = "_"), function(x) x[3]))
} else if (pos_type == "gen_dist") {
  marker_map <- merge(marker_map, nam_marc_marker_info, by = "snp")
  colnames(marker_map)[2:3] <- c("Chr", "pos")
}

#grabbing an individual phenotype
ind_pheno <- pheno_name[run]

#reading in the stepwise model
pheno_sw_model_files <- list.files(pattern = ind_pheno, 
                                   recursive = TRUE, 
                                   include.dirs = TRUE)
pheno_sw_model_files <- pheno_sw_model_files[grepl("stepwise_model", pheno_sw_model_files)]

test_model <- readRDS(pheno_sw_model_files)

#need to read in all of the markers
complete_markers <- Genotypes_unpruned
names(complete_markers) <- paste("m", names(complete_markers), sep = "_")
marker_info <- data.frame(marker = names(complete_markers), stringsAsFactors = FALSE)
marker_info$chr <- sapply(strsplit(marker_info$marker, split = "_"), function(x) as.numeric(x[2]))

#grabbing the vector of QTL
test_model_qtl <- test_model$found_qtl

#function that will return a data frame with markers on the same chromosome
same_chromosome_markers <- function(marker, marker_df) {
  #getting the marker chromosome
  marker_chr <- as.numeric(sapply(strsplit(marker, split = "_"), function(x) x[2]))
  
  #subsetting the marker data frame
  marker_df_subset <- subset(marker_df, chr == marker_chr)
  
  #returning the data frame
  return(marker_df_subset)
}

#function that will return the marker to the left or right of a marker
marker_shift <- function(marker, marker_chr_subset, shift_direction = "left", steps = 1) {
  #vector of markers
  marker_chr_vector <- marker_chr_subset$marker
  
  #finding the position of the testing marker in the vector
  marker_pos <- which(marker_chr_vector == marker)
  
  if (shift_direction == "left") {
    shift <- -1 * steps
  } else if (shift_direction == "right") {
    shift <- 1 * steps
  }
  
  #just checking to make sure the shifts don't return an error: marker shift positions 
  #are not out of the range of the chromosome
  marker_shift_pos <- marker_pos + shift
  
  if (marker_shift_pos < 1) {
    marker_shift_pos <- 1
  } else if (marker_shift_pos > length(marker_chr_vector)) {
    marker_shift_pos <- length(marker_chr_vector)
  }
  
  #grabbing the shifted marker
  marker_shift <- marker_chr_vector[marker_shift_pos]
  
  #if the shifted marker is the same as the marker, return nothing and stop the
  #pipeline later
  if (marker_shift == marker) {
    return(NULL)
  } else {
    return(marker_shift)
  }
}

#function that will take the complete marker set and then subset out for 
#individual populations, and the reduced and complete marker set
#outputs the phenotype data, appended with the genotype data
full_marker_data_subsetter <- function(full_marker_vector, marker_array, pheno_data) {
  #subsetting the marker array with only the qtl markers
  marker_array_subset <- marker_array[full_marker_vector]
  
  #now need to grab the genotypes in the population being tested
  marker_array_pop_subset <- mclapply(marker_array_subset, function(x) x[match(pheno_data$geno, rownames(x)), ])
  
  #grabbing only the col geno probabilities, and then cbinding into a data frame
  marker_pop_subset_df <- do.call(cbind, mclapply(marker_array_pop_subset, function(x) x[1]))
  
  #need to scale the snps or else the test marker (bound testing marker) will not
  #be scaled while the other markers will be
  marker_pop_subset_df <- scale_SNPs(marker_pop_subset_df)
  
  #naming the marker df with the qtl names
  colnames(marker_pop_subset_df) <- full_marker_vector
  
  #combinging the phenotype data with the genotype data
  pheno_with_geno <- cbind(pheno_data, marker_pop_subset_df)
  
  return(pheno_with_geno)
}

#function that will remove markers from the pre-specified kinship matrix
pop_kinship_generator_proximal <- function(pruned_geno_data, proximal_matrix) {
  pruned_geno_proximal_subset <- pruned_geno_data[, proximal_matrix]
  pop_kinship_after_proximal <- tcrossprod(scale_SNPs(pruned_geno_proximal_subset, center = T, scale = T)) / ncol(pruned_geno_proximal_subset)
  return(pop_kinship_after_proximal)
}

#need to change the proximal matrices in all_pop_data for Dan's new GridLMM_GWAS function
for (j in 1:length(all_pop_data)) {
  #list of markers to remove in the kinship matrix
  all_pop_data[[j]]$cM_proximal_list <- apply(all_pop_data[[j]]$cM_proximal, 1, function(x) which(x == 1))
  
  #scaling the snps so that i don't have to scale them in the function
  all_pop_data[[j]]$pop_geno <- scale_SNPs(all_pop_data[[j]]$pop_geno)
}

#making a data frame to hold all the qtl confidence intervals for a particular trait
trait_qtl_cis <- c()

#a for loop that will generate a confidence interval for each marker
for (i in 1:length(test_model_qtl)) {
  #grabbing the marker
  test_marker <- test_model_qtl[i]
  
  print(paste("=============", "Currently generating confidence intervals for marker", test_marker, "=============", sep = " "))
  
  #need to get the marker position for later model fitting (downdat X's)
  test_marker_pos <- which(colnames(pop_markers) == test_marker)
  
  #generating a df of markers on the same chromosome as the testing_marker
  test_marker_same_chr <- same_chromosome_markers(test_marker, marker_info)
  
  #grabbing the proximal matrices for all pops just for the test_marker
  all_pop_test_marker_proximal_matrices <- lapply(all_pop_data, function(x) {
    #grabbing the proximal matrix
    pop_proximal <- as.data.frame(x$cM_proximal)
    
    #naming the rows of the proximal matrix
    rownames(pop_proximal) <- colnames(pop_markers)
    
    #grabbing only the test marker from the pro
    test_marker_proximal <- subset(pop_proximal, rownames(pop_proximal) == test_marker)
    
    #need to change the test_marker_proximal matrices to fit with Dan's new GridLMM_GWAS function
    test_marker_proximal <- which(as.vector(test_marker_proximal[1, ]) == 1)
    
    return(test_marker_proximal)
  })
  
  #formula for the first scan
  first_scan_formula <- as.formula(paste(ind_pheno, " ~ (1|geno)", sep = ""))
  
  #running the scan over all markers over all populations
  first_marker_scan <- lapply(1:length(all_pop_data), function(x) GridLMM_GWAS(first_scan_formula, ~1, ~0,
                                                                               data = all_pop_data[[x]]$pop_pheno,
                                                                               X = all_pop_data[[x]]$pop_geno,
                                                                               X_ID = 'geno',
                                                                               proximal_markers = all_pop_data[[x]]$cM_proximal_list,
                                                                               h2_step = 0.01,
                                                                               max_steps = 100,
                                                                               mc.cores = 1,
                                                                               centerX = FALSE,
                                                                               scaleX = FALSE,
                                                                               method = "ML"))
  
  #grabbing the single marker for each individual population
  all_pop_test_marker <- lapply(all_pop_data, function(x) subset(x$pop_geno, select = test_marker))
  
  #kinship matrix 
  #all_pop_kinship_matrix <- lapply(1:length(all_pop_data), function(x) tcrossprod(scale_SNPs(all_pop_data[[x]]$pop_geno, center = T, scale = F)) / ncol(all_pop_data[[x]]$pop_geno))
  
  #first finding the left side
  left_side <- FALSE
  
  step_size <- 1
  
  print("============= Slide to the left =============")
  while (!left_side) {
    #finding markers just to the right or left of a marker
    shifted_marker <- marker_shift(test_marker, test_marker_same_chr, "left", steps = step_size)
    
    #generating a vector of markers for the full and reduced models; the formulas
    #aren't generated until the markers are assigned the the global environment
    full_qtl_vector <- c(test_model_qtl, shifted_marker)
    
    #generating a vector of markers for the reduced model
    reduced_qtl_vector <- test_model_qtl
    reduced_qtl_vector[[which(test_model_qtl == test_marker)]] <- shifted_marker
    
    #generating information for all populations
    all_pop_pheno_with_geno <- mclapply(all_pop_data, function(x) full_marker_data_subsetter(full_marker_vector = reduced_qtl_vector,
                                                                                             marker_array = complete_markers,
                                                                                             pheno_data = x$pop_pheno))
    
    #generating a formula for the reduced qtl list
    reduced_qtl_formula <- paste(reduced_qtl_vector, collapse = " + ")
    reduced_model_formula <- as.formula(paste(paste(ind_pheno, reduced_qtl_formula, sep = " ~ "), "(1|geno)", sep = " + "))
    
    #printing out the marker name; this is both to see the progress of the confidence intervals
    #and also to see if there is an error with the script or Dan's package
    print(paste("=============", paste("The current marker is ", shifted_marker, sep = " "), "=============", sep = ""))
    
    #reduced model fit testing all markers; i don't know how to test only 1 marker 
    all_pop_reduced_tests <- lapply(1:length(all_pop_pheno_with_geno), function(x) GridLMM_GWAS(reduced_model_formula, ~1, ~0,
                                                                                                data = all_pop_pheno_with_geno[[x]],
                                                                                                X = all_pop_test_marker[[x]],
                                                                                                X_ID = 'geno',
                                                                                                h2_step = 0.01,
                                                                                                max_steps = 100,
                                                                                                mc.cores = 1,
                                                                                                centerX = FALSE,
                                                                                                scaleX = FALSE,
                                                                                                V_setup = first_marker_scan[[x]]$setup$V_setup,
                                                                                                proximal_markers = all_pop_data[[x]]$cM_proximal_list[test_marker_pos],
                                                                                                proximal_Xs = list(all_pop_data[[x]]$pop_geno),
                                                                                                method = "ML"))
    
    #grabbing the results for all markers
    pop_reduced_tests_dt <- rbindlist(lapply(all_pop_reduced_tests, function(x) subset(x$results)))
    test_marker_sig <- with(pop_reduced_tests_dt, pchisq(2 * (sum(ML_logLik) - sum(ML_Reduced_logLik)), df = length(all_pop_data), lower.tail = FALSE))
    
    #if the test is significant for the testing marker, or if the shifted marker is the same as
    #the testing marker, then stop the while loop
    if (test_marker_sig < 0.05 | shifted_marker == test_marker_same_chr$marker[1]) {
      left_side <- TRUE
      left_bound <- shifted_marker
    } else {
      step_size <- step_size + 1
    }
  }
  
  #grabbing the marker
  test_marker <- test_model_qtl[i]
  
  #need to get the marker position for later model fitting (downdat X's)
  test_marker_pos <- which(colnames(pop_markers) == test_marker)
  
  #generating a df of markers on the same chromosome as the testing_marker
  test_marker_same_chr <- same_chromosome_markers(test_marker, marker_info)
  
  #grabbing the proximal matrices for all pops just for the test_marker
  all_pop_test_marker_proximal_matrices <- lapply(all_pop_data, function(x) {
    #grabbing the proximal matrix
    pop_proximal <- as.data.frame(x$cM_proximal)
    
    #naming the rows of the proximal matrix
    rownames(pop_proximal) <- colnames(pop_markers)
    
    #grabbing only the test marker from the pro
    test_marker_proximal <- subset(pop_proximal, rownames(pop_proximal) == test_marker)
    
    #need to change the test_marker_proximal matrices to fit with Dan's new GridLMM_GWAS function
    test_marker_proximal <- which(as.vector(test_marker_proximal[1, ]) == 1)
    
    return(test_marker_proximal)
  })
  
  #first finding the left side
  right_side <- FALSE
  
  step_size <- 1
  
  print("============= Slide to the right =============")
  
  while (!right_side) {
    #finding markers just to the right or left of a marker
    shifted_marker <- marker_shift(test_marker, test_marker_same_chr, "right", steps = step_size)
    
    #generating a vector of markers for the full and reduced models; the formulas
    #aren't generated until the markers are assigned the the global environment
    full_qtl_vector <- c(test_model_qtl, shifted_marker)
    
    #generating a vector of markers for the reduced model
    reduced_qtl_vector <- test_model_qtl
    reduced_qtl_vector[[which(test_model_qtl == test_marker)]] <- shifted_marker
    
    #generating information for all populations
    all_pop_pheno_with_geno <- mclapply(all_pop_data, function(x) full_marker_data_subsetter(full_marker_vector = reduced_qtl_vector,
                                                                                             marker_array = complete_markers,
                                                                                             pheno_data = x$pop_pheno))
    
    #generating a formula for the reduced qtl list
    reduced_qtl_formula <- paste(reduced_qtl_vector, collapse = " + ")
    reduced_model_formula <- as.formula(paste(paste(ind_pheno, reduced_qtl_formula, sep = " ~ "), "(1|geno)", sep = " + "))
    
    print(paste("=============", paste("The current marker is ", shifted_marker, sep = " "), "=============", sep = ""))
    
    #reduced model fit testing all markers; i don't know how to test only 1 marker 
    all_pop_reduced_tests <- lapply(1:length(all_pop_pheno_with_geno), function(x) GridLMM_GWAS(reduced_model_formula, ~1, ~0,
                                                                                                data = all_pop_pheno_with_geno[[x]],
                                                                                                X = all_pop_test_marker[[x]],
                                                                                                X_ID = 'geno',
                                                                                                h2_step = 0.01,
                                                                                                max_steps = 100,
                                                                                                mc.cores = 1,
                                                                                                centerX = FALSE,
                                                                                                scaleX = FALSE,                                                                                                V_setup = first_marker_scan[[x]]$setup$V_setup,
                                                                                                proximal_markers = all_pop_data[[x]]$cM_proximal_list[test_marker_pos],
                                                                                                proximal_Xs = list(all_pop_data[[x]]$pop_geno),
                                                                                                method = "ML"))
    
    #grabbing the results for all markers
    pop_reduced_tests_dt <- rbindlist(lapply(all_pop_reduced_tests, function(x) subset(x$results)))
    test_marker_sig <- with(pop_reduced_tests_dt, pchisq(2 * (sum(ML_logLik) - sum(ML_Reduced_logLik)), df = length(all_pop_data), lower.tail = FALSE))
    
    #if the test is significant for the testing marker, or if the shifted marker is the same as
    #the testing marker, then stop the while loop
    if (test_marker_sig < 0.05 | shifted_marker == test_marker_same_chr$marker[nrow(test_marker_same_chr)]) {
      right_side <- TRUE
      right_bound <- shifted_marker
    } else {
      step_size <- step_size + 1
    }
  }
  #creating a data frame with the qtl and its bounds, and then adding it to the 
  #dt of qtl cis
  qtl_ci <- data.frame(qtl = test_marker, left_bound = left_bound, right_bound = right_bound)
  trait_qtl_cis <- rbindlist(list(trait_qtl_cis, qtl_ci))
}

#now saving the confidence intervals as a csv file
trait_qtl_ci_file_name <- paste(ind_pheno, "qtl_ci.csv", sep = "_")
fwrite(trait_qtl_cis, trait_qtl_ci_file_name, sep = ",", row.names = FALSE, col.names = TRUE)
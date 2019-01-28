#reading in the genotype probability array
Genotypes <- readRDS("nam_rqtl_geno_prob_array_final_pheno_subset.RDS")

#reading in the genotype line names from the genotypes array rows
geno_line_names <- rownames(Genotypes[[1]])

#function that generates a kinship matrix using the genotype probabilities
kinship_matrix_maker <- function(marker_probability_array, pop_pop_name) {
  #grabbing the rows for individual populations 
  pop_pop_marker_array <- mclapply(marker_probability_array, function(x) subset(x, grepl(pop_pop_name, rownames(x))))
  
  #now grabbing just the col-0 states from each of the arrays
  pop_pop_marker_array_col <- as.matrix(do.call(cbind, mclapply(pop_pop_marker_array, function(x) as.data.frame(x[, 1]), mc.cores = detectCores())))
  
  #naming the rows of the newly created col df for the population
  rownames(pop_pop_marker_array_col) <- rownames(pop_pop_marker_array[[1]])
  
  #now need to center and sweep
  pop_pop_marker_array_kinship_cen <- sweep(pop_pop_marker_array_col, 2, colMeans(pop_pop_marker_array_col), '-')
  pop_pop_marker_array_kinship <- (pop_pop_marker_array_kinship_cen %*% t(pop_pop_marker_array_kinship_cen)) / ncol(pop_pop_marker_array_kinship_cen)
  
  return(pop_pop_marker_array_kinship)
}
  
kinship <- "no"

Genotypes_copy <- mclapply(Genotypes, function(x) as.matrix(x) %*% contr.sum(8))
  
#creating an empty vector that will tell me the bad vectors with perfectly correlated markers
bad_markers <- c()

#running a for loop that will calculate the corrleation of the marker with other markers
#and then drop the markers that are perfectly correlated with the one currently being tested
#the idea is that we will create a vector of bad markers 
cor_threshold <- 0.99

#setting a cor test list to keep track of markers and their correlations
cor_vector_list <- list()

#testing out the correlated marker for loop
for (i in names(Genotypes_copy)) {
  print(paste("Currently on marker", i, sep = " "))
  
  #grabbing the marker for testing
  test_marker <- as.vector(Genotypes_copy[[i]])
  
  #grabbing the position to omit it from the analysis
  test_marker_position <- which(names(Genotypes_copy) == i)
  
  #getting the position in the original matrix so that i can reference the position
  #when making a new list
  test_marker_position_original <- which(names(Genotypes) == i)
  
  #grabbing the subset 
  Genotypes_subset <- Genotypes_copy[-test_marker_position]
  
  #calculating the correlation of the marker with the other markers
  cor_vector <- unlist(mclapply(1:length(Genotypes_subset), function(x) cor(test_marker, as.vector(Genotypes_subset[[x]])), mc.cores = detectCores()))
  
  #this will create a new element using the position of the current testing marker
  #and store the correlation vector for later testing
  cor_vector_list[[test_marker_position_original]] <- cor_vector
  
  #if any of the markers are perfectly correlated with the markers, besides the 
  if (any(abs(cor_vector) >= cor_threshold, na.rm = TRUE)) {
    bad_markers <- c(bad_markers, i)
    Genotypes_copy <- Genotypes_copy[-test_marker_position]
  }
}

#naming the markers
names(cor_vector_list) <- names(Genotypes)

#list with original probability array, pruned probability array, list of correlated markers,
#and correlation lists
bad_markers_list <- list(original_genotypes_rotated = Genotypes,
                         pruned_genotypes_rotated = Genotypes_copy,
                         bad_markers = bad_markers,
                         correlation_list = cor_vector_list)

#generating the file names and then saving them as RDS files
file_name <- paste("bad_markers_list_final_", cor_threshold, "_", paste(kinship, "kinship", sep = "_"), ".rds", sep = "")
marker_name <- paste("nam_rqtl_geno_prob_array_final_", cor_threshold, "_", paste(kinship, "kinship", sep = "_"), ".rds", sep = "")

saveRDS(bad_markers_list, file_name)
saveRDS(Genotypes_copy, marker_name)


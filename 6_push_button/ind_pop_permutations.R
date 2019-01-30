#reading in the pruned genotypes, and will subset the complete genotypes later by the markers remaining in the pruned set
Genotypes_pruned <- readRDS("nam_rqtl_geno_prob_array_final_0.99_no_kinship.rds")

#READIN IN THE MARKERS FOR genotypes only found within the blups data set
#reading in the unpruned data set and grabbing the markers that were not pruned
Genotypes_unpruned <- readRDS("nam_rqtl_geno_prob_array_final_pheno_subset.RDS")

#reading in all of the phenotype data; need to read this in first to subset genotypes and kinship matrix
pheno_type <- "univariate"

#reading in different blup files based on if i want to use a univariate or a multivariate approach
if (pheno_type == "univariate") {
  #reading in the blups file, that were generated from a linear mixed model that was only
  #using the univariate mixed models in brms
  nam_pheno <- fread("nam_blups_combined_univariate.csv",
                     sep = ",",
                     header = TRUE,
                     stringsAsFactors = FALSE)
} else if (pheno_type == "multivariate") {
  #reading in the blups file, that are based on a multivariate model in brms; this model
  #has correlated residuals between the different traits
  nam_pheno <- fread("nam_blups_combined_final.csv", 
                     sep = ",", 
                     header = TRUE, 
                     stringsAsFactors = FALSE)
}

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
nam_marc_marker_info <- fread("nam_marker_info_final.csv",
                              sep = ",", 
                              header = TRUE, 
                              stringsAsFactors = FALSE)
colnames(nam_marc_marker_info)[1] <- "snp"
nam_marc_marker_info$snp <- paste("m", nam_marc_marker_info$snp, sep = "_")

#NEED TO SAVE THIS AS AN R OBJECT SO I DON'T HAVE TO KEEP RUNNING IT OVER AND OVER AGAIN
all_pop_data <- readRDS("nam_all_traits_ind_pop_pheno_geno_proximal_final.RDS")

#grabbing the individual population data
pop_name <- names(all_pop_data)[[run]]
pop_data <- all_pop_data[[pop_name]]

#grabbing the markers and then doing stuff
pop_markers <- pop_data$pop_geno

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

#making proximal matrix specific to population
if (proximal_type == "cor") {
  #grabbing the proximal matrix based off marker correlations
  pop_marker_cors <- pop_data$marker_cors
} else if (proximal_type == "window") {
  #grabbing the proximal matrix based off a sliding window
  pop_marker_cors <- pop_data$cM_proximal
}

#this is only needed if not done previously
#pop_marker_cors <- as.matrix(rbindlist(lapply(pop_marker_cors, function(x) as.data.frame(t(as.matrix(x))))))
colnames(pop_marker_cors) <- marker_map$snp

pop_pheno <- pop_data$pop_pheno
pop_pheno_only <- subset(pop_pheno, select = grepl("_geno|_gxe", colnames(pop_pheno)))
pop_pheno_info <- subset(pop_pheno, select = geno)

try(dir.create(pop_name))
setwd(pop_name)

#need to add this step since the new gemma gwas
pop_marker_cors <- apply(pop_marker_cors, 1, function(x) which(x == 1))

#for each phenotype do 1000 permutations, and save the output
for (j in 1:ncol(pop_pheno_only)) {
  #grabbing the phenotype name
  pheno_name <- names(pop_pheno_only)[j]
  
  #generating results for an individual population
  pop_perms <- vector("list", length = perm_no)
  pop_ind_pheno <- pop_pheno_only[,..j]
  
  pop_pheno_subset <- cbind(pop_pheno_info, pop_ind_pheno)
  
  model_formula <- as.formula(paste(pheno_name, " ~ (1|geno)", sep = ""))
  
  #remove correlated markers proximal matrix
  m1_normal_scan <- GridLMM_GWAS(model_formula, ~1, ~0,
                                 data = pop_pheno_subset,
                                 X = pop_markers,
                                 X_ID = 'geno',
                                 X_map = marker_map,
                                 proximal_markers = pop_marker_cors,
                                 h2_step = .01,
                                 max_steps = 100,
                                 mc.cores = 1,
                                 centerX = FALSE,
                                 scaleX = FALSE,
                                 method = "ML")
  
  m1_normal_scan$setup$V_setup = set_p_test(m1_normal_scan$setup$V_setup, 
                                            list(geno = 0),
                                            list(geno = ncol(pop_markers)))
  
  pop_markers_scrambled <- pop_markers
  
  registerDoParallel(max_cores)
  start_time <- Sys.time()
  
  pop_perms <- foreach(i = 1:perm_no) %dopar% {
    print(paste("Perm no", i, sep = " "))

    model_formula <- as.formula(paste(pheno_name, " ~ (1|geno)", sep = ""))
    
    #changing the permutation type - phenotype or genotype
    if (perm_type == "pheno") {
      #this part is for the phenotypic permutations
      pop_pheno_subset[, 2] <- pop_pheno_subset[sample(1:nrow(pop_pheno_subset), size = nrow(pop_pheno_subset)), 2]
      
      #remove correlated markers proximal matrix
      m1_cor_mark <- GridLMM_GWAS(model_formula, ~1, ~0,
                                  data = pop_pheno_subset,
                                  X = pop_markers,
                                  X_ID = 'geno',
                                  X_map = marker_map,
                                  #proximal_markers = NULL, NULL only if we're doing genotype permutations
                                  proximal_markers = pop_marker_cors,
                                  h2_step = .01,
                                  max_steps = 100,
                                  mc.cores = 1,
                                  V_setup = m1_normal_scan$setup$V_setup,
                                  centerX = FALSE,
                                  scaleX = FALSE,
                                  method = "ML")
    } else if (perm_type == "geno") {
      #scrambling the phenotypes
      rownames(pop_markers_scrambled) <- sample(rownames(pop_markers_scrambled))
      
      #remove correlated markers proximal matrix
      m1_cor_mark <- GridLMM_GWAS(model_formula, ~1, ~0,
                                  data = pop_pheno_subset,
                                  X = pop_markers_scrambled,
                                  X_ID = 'geno',
                                  X_map = marker_map,
                                  proximal_markers = NULL, #NULL only if we're doing genotype permutations
                                  h2_step = .01,
                                  max_steps = 100,
                                  mc.cores = 1,
                                  V_setup = m1_normal_scan$setup$V_setup,
                                  centerX = FALSE,
                                  scaleX = FALSE,
                                  method = "ML")
    }
    pop_results <- m1_cor_mark$results
    return(pop_results)
  }
  end_time <- Sys.time()
  print(paste("Perms", "took", as.numeric(end_time - start_time), sep = " "))
  pop_perms_name <- paste(pop_name, pheno_name, pheno_type, paste(perm_no, "x", sep = ""), proximal_type, perm_type, "perms.RDS", sep = "_")
  saveRDS(pop_perms, pop_perms_name)
}

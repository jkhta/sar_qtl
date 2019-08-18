library(dplyr)
library(lme4)
library(lmerTest)
library(parallel)
library(data.table)
library(lme4qtl)
library(qtl)
library(pryr)
library(GridLMM)
library(ggplot2)

rm(list = ls())

max_cores <- 5
run <- as.numeric(commandArgs(t = T)[1])

###START FROM HERE AFTER SAVING THE GENOTYPE PROBABILITY ARRAY
#reading in the kinship matrix - FROM TASSEL; this is using the centered IBS method
setwd("/home/jkta/projects/col_sha/complete_pipeline/input_output")

#reading in the pruned genotypes, and will subset the complete genotypes later by the markers remaining in the pruned set
Genotypes_pruned <- readRDS("nam_rqtl_geno_prob_array_final_0.99_no_kinship.rds")

#READIN IN THE MARKERS FOR ALL GENOTYPES
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
pheno_name <- colnames(all_pop_data$`21RV_21RV`$pop_pheno)
pheno_name <- pheno_name[grepl("_geno|_gxe", pheno_name)]

#grabbing the markers to generate a marker map
pop_markers <- all_pop_data$`21RV_21RV`$pop_geno

#generating a marker map: it can be physical or genetic depending on the choices
marker_map <- data.frame(snp = colnames(pop_markers), stringsAsFactors = FALSE)

pos_type = "gen_dist"

#generating markers depending on physical or genetic distance
if (pos_type == "phys_dist") {
  marker_map$Chr <- as.character(sapply(strsplit(marker_map$snp, split = "_"), function(x) x[2]))
  marker_map$pos <- as.numeric(sapply(strsplit(marker_map$snp, split = "_"), function(x) x[3]))
} else if (pos_type == "gen_dist") {
  marker_map <- merge(marker_map, nam_marc_marker_info, by = "snp")
  colnames(marker_map)[2:3] <- c("Chr", "pos")
}

#this function runs GxEMMA with a stepwise algorithm
GridLMM_scan <- function(pheno, pheno_name, geno, map, proximal, sig_qtl) {
  pheno_copy <- pheno
  
  #since the new GridLMM update, we need to change the proximal markers into a list
  #instead of a matrix
  proximal <- apply(proximal, 1, function(x) which(x == 1))
  
  #if we already have significant QTL, then add them as fixed effects in the model
  if (length(sig_qtl) < 1) {
    model_formula <- as.formula(paste(pheno_name, " ~ ", "(1|geno)", sep = ""))
  } else if (length(sig_qtl) >= 1) {
    pheno_copy <- cbind(pheno_copy, subset(geno, select = colnames(geno) %in% sig_qtl))
    sig_qtl_formula <- paste(sig_qtl, collapse = " + ")
    model_formula <- as.formula(paste(pheno_name, " ~ ", sig_qtl_formula, "+ (1|geno)", sep = ""))
  }
  
  #running the GxEMMA code
  m1_cor_mark <- GridLMM_GWAS(model_formula, ~1, ~0,
                              data = pheno_copy,
                              X = geno,
                              X_ID = 'geno',
                              X_map = map,
                              proximal_markers = proximal,
                              h2_step = 0.01,
                              max_steps = 100,
                              mc.cores = max_cores,
                              centerX = TRUE, 
                              method = "ML")
  
  #getting the results and the most significant snp
  qtl_df <- m1_cor_mark$results
  
  return(qtl_df)
}

#this function runs GxEMMA with a stepwise algorithm
GridLMM_stepwise <- function(list_of_pop_things, pheno_name, threshold, max_markers) {
  #making an empty score function
  qtl_df <- c()
  
  sig_qtl <- c()
  
  qtl_max_scores <- c()
  #signifying if thi is the first step
  first_step <- TRUE
  
  print(paste("Currently have", length(sig_qtl), "QTL", sep = " "))
  
  while (first_step | (length(sig_qtl) < max_markers & max(qtl_df$score) > threshold)) {
    #running individual population scans for all populations: without a kinship matrix
    ind_pop_scans <- lapply(list_of_pop_things, function(x) GridLMM_scan(pheno = x$pop_pheno, 
                                                                         pheno_name = pheno_name, 
                                                                         geno = x$pop_geno, 
                                                                         map = marker_map, 
                                                                         proximal = x$marker_cors,
                                                                         sig_qtl = sig_qtl))
    
    #data frame with all scores 
    qtl_df <- data.frame(snp = ind_pop_scans[[1]]$snp,
                         chr = ind_pop_scans[[1]]$Chr,
                         pos = ind_pop_scans[[1]]$pos,
                         score = -log10(pchisq(rowSums(do.call(cbind, lapply(ind_pop_scans, function(x) as.data.frame(2 * (x$ML_logLik - x$ML_Reduced_logLik))))), df = 7, lower.tail = FALSE)),
                         stringsAsFactors = FALSE)
    
    #we already ran the model for the first time, so need to set to false 
    first_step <- FALSE
    
    #getting the results and the most significant snp
    qtl_df_max <- qtl_df[which.max(qtl_df$score), ]
    
    #if the snp significance is over the threshold, add it to the vector of significant snps 
    if (qtl_df_max$score > threshold) {
      sig_qtl <- c(sig_qtl, qtl_df_max$snp)
    }
    qtl_max_scores <- rbindlist(list(qtl_max_scores, qtl_df_max))
  }
  #subsetting the list of found qtl by the ones that are above the threshold
  qtl_max_scores <- subset(qtl_max_scores, score > threshold)
  
  #returning last model scans, last scan scores, all found qtl scores, and a vector of found QTL
  return(list(last_scan_models = ind_pop_scans,
              last_scan = qtl_df,
              qtl_scores = qtl_max_scores,
              found_qtl = sig_qtl))
}

#grabbing an individual phenotype
ind_pheno <- pheno_name[run]

setwd("/home/jkta/projects/col_sha/complete_pipeline/input_output/permutations")

#reading all of the permutations into a list
pheno_perms_list <- list.files(pattern = ind_pheno, 
                               recursive = TRUE, 
                               include.dirs = TRUE)

perm_type <- "geno"

if (perm_type == "geno") {
  #grabbing only the permutations that permuted the markers
  pheno_perms_list <- pheno_perms_list[grepl("geno_perms", pheno_perms_list)]
} else if (perm_type == "pheno") {
  #grabbing only the permutations that permuted the phenotypes
  pheno_perms_list <- pheno_perms_list[grepl("window_pheno", pheno_perms_list)]
}

#reading the permutations for each population
pheno_perms <- lapply(pheno_perms_list, function(x) readRDS(x))

#calculating the difference between the full and reduced model log-likelihoods
pheno_perm_lik <- lapply(pheno_perms, function(x) lapply(x, function(y) as.data.frame(2 * (y$ML_logLik - y$ML_Reduced_logLik))))

#grabbing the corresponding permutation for from each population, all 1st, 2nd, 3rd permutations, etc.
pheno_perms_together <- lapply(1:length(pheno_perm_lik[[1]]), function(x) do.call(cbind, lapply(pheno_perm_lik, function(y) as.data.frame(y[[x]]))))

#calculating the sum of all log-likelihood scores
pheno_perm_scores <- -log10(pchisq(unlist(lapply(pheno_perms_together, function(x) max(rowSums(x)))), df = length(pheno_perms), lower.tail = FALSE))
pheno_perm_threshold <- quantile(pheno_perm_scores, probs = 0.95)

pheno_stepwise_model <- GridLMM_stepwise(all_pop_data, 
                                         pheno_name = ind_pheno,
                                         threshold = pheno_perm_threshold, 
                                         max_markers = 10)

setwd("/home/jkta/projects/col_sha/complete_pipeline/input_output/stepwise_models")
pheno_stepwise_file_name <- paste(ind_pheno, pheno_type, "GridLMM_stepwise_model.RDS", sep = "_")

saveRDS(pheno_stepwise_model, pheno_stepwise_file_name)
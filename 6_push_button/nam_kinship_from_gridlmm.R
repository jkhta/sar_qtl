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
library(foreach)
library(doParallel)

rm(list = ls())

#directory with genotype probability array and blups
setwd("/home/jkta/projects/col_sha/complete_pipeline/")

#-----STEP 1 START-----#

#-----INPUT-----#
#reading in the pruned genotypes, and will subset the complete genotypes later by the markers remaining in the pruned set
Genotypes_pruned <- readRDS("nam_rqtl_geno_prob_array_james_0.99_pruned.rds")

#READING IN THE MARKERS FOR ALL GENOTYPES
#reading in the unpruned data set and grabbing the markers that were not pruned
Genotypes_unpruned <- readRDS("nam_rqtl_geno_prob_array_james_pheno_subset.RDS")

#reading in all of the phenotype data; need to read this in first to subset genotypes and kinship matrix
pheno_type <- "james"

#reading in different blup files based on if i want to use a univariate or a multivariate approach
nam_pheno <- fread("nam_blups_combined_univariate.csv",
                   sep = ",",
                   header = TRUE,
                   stringsAsFactors = FALSE)

#need to merge file with genetic distances to plot genetic distances
nam_marc_marker_info <- fread("nam_marker_info_final.csv",
                              sep = ",", 
                              header = TRUE, 
                              stringsAsFactors = FALSE)

all_pop_data <- readRDS("nam_all_traits_ind_pop_pheno_geno_proximal_james.RDS")

#-----PARAMETERS-----#
#number of cores to run the permutations
max_cores <- 5

#run is a number that corresponds to a trait; goes from 1-8 for
#the 8 different traits; each number is run on a different cluster
run <- as.numeric(commandArgs(t = T)[1])

phenotype_list <- "_geno|_gxe"

#generating a marker map: it can be physical or genetic depending on the choices
pos_type = "gen_dist"

perm_type <- "geno"

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

#this function runs GxEMMA with a stepwise algorithm
GridLMM_scan <- function(pheno, geno, map, proximal, model_formula) {
  pheno_copy <- pheno
  
  #since the new GridLMM update, we need to change the proximal markers into a list
  #instead of a matrix
  proximal <- apply(proximal, 1, function(x) which(x == 1))
  
  #running the GxEMMA code
  m1_cor_mark <- GridLMM_GWAS(model_formula, ~1, ~0,
                              data = pheno_copy,
                              X = geno,
                              X_ID = 'geno',
                              X_map = map,
                              proximal_markers = proximal,
                              h2_step = 0.01,
                              max_steps = 100,
                              mc.cores = 1,
                              centerX = TRUE, 
                              method = "ML")
  
  #getting the results and the most significant snp
  qtl_df <- m1_cor_mark$results
  qtl_df <- m1_cor_mark
  
  return(qtl_df)
}

ind_pop_scans <- lapply(all_pop_data, function(x) GridLMM_scan(pheno = x$pop_pheno, 
                                                               geno = x$pop_geno, 
                                                               map = marker_map, 
                                                               proximal = x$marker_cors,
                                                               model_formula = "r_dry_gxe ~ bd_gxe + (1|geno)"))

kinship_list <- lapply(ind_pop_scans, function(x) x$setup$V_setup$RE_setup$geno$K)

saveRDS(kinship_list, "nam_GridLMM_kinship.RDS")


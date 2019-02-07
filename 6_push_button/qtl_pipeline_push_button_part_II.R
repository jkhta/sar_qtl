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

#PART II

#directory with genotype probability array and blups
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/6_push_button/")

#-----STEP 1 START-----#

#-----INPUT-----#
#reading in the pruned genotypes, and will subset the complete genotypes later by the markers remaining in the pruned set
Genotypes_pruned <- readRDS("nam_rqtl_geno_prob_array_final_0.99_no_kinship.rds")

#READIN IN THE MARKERS FOR genotypes only found within the blups data set
#reading in the unpruned data set and grabbing the markers that were not pruned
Genotypes_unpruned <- readRDS("nam_rqtl_geno_prob_array_james_pheno_subset.RDS")

pheno_type <- "james"

#reading in all of the phenotype data; need to read this in first to subset genotypes and kinship matrix
nam_pheno <- fread("nam_blups_combined_univariate.csv",
                   sep = ",",
                   header = TRUE,
                   stringsAsFactors = FALSE)

#need to merge file with genetic distances to plot genetic distances
nam_marc_marker_info <- fread("nam_marker_info_final.csv",
                              sep = ",", 
                              header = TRUE, 
                              stringsAsFactors = FALSE)

#reading in list of phenotype, genotype, and proximal data
all_pop_data <- readRDS("nam_all_traits_ind_pop_pheno_geno_proximal_james.RDS")


#-----PARAMETERS-----#
#number of cores to run the permutations
max_cores <- detectCores()

#run is a number that corresponds to a population; goes from 1-7 for
#the 7 different populations; each number is run on a different cluster
run <- as.numeric(commandArgs(t = T)[1])

#number of permutations
perm_no <- 1000

#permutation type - genotype or phenotype permutations
perm_type <- "geno"

#controls output distance; physical = phys_dist or genetic = gen_dist
pos_type = "gen_dist"

#proximal type: cor for removing correlated markers in the calculation of the
#kinship matrix or window - based on a centimorgan distance
proximal_type <- "window"

#regular expression for phenotypes which will subset the data frame
phenotype_list <- "_geno|_gxe"

#script for running permutations for each population/trait
source("ind_pop_permutations.R")

#-----STEP 1 END-----#
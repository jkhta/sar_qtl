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

#PART III

#directory with genotype probability array and blups
setwd("/home/jkta/projects/col_sha/complete_pipeline_path_analysis")

#-----STEP 1 START-----#

#-----INPUT-----#
#reading in the pruned genotypes, and will subset the complete genotypes later by the markers remaining in the pruned set
#Genotypes_pruned <- readRDS("nam_rqtl_geno_prob_array_james_0.99_pruned.rds")
Genotypes_pruned <- readRDS("nam_rqtl_geno_prob_array_james_0.99_pruned_revised.rds")

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

#-----RUNNING THE SCRIPT-----#
source("ind_trait_stepwise_QTL.R")

#-----STEP 1 END-----#

#-----STEP 2 START-----#
rm(list = ls())

#-----INPUT-----#
#reading in the pruned genotypes, and will subset the complete genotypes later by the markers remaining in the pruned set
#Genotypes_pruned <- readRDS("nam_rqtl_geno_prob_array_james_0.99_pruned.rds")
Genotypes_pruned <- readRDS("nam_rqtl_geno_prob_array_james_0.99_pruned_revised.rds")

#READING IN THE MARKERS FOR ALL GENOTYPES
#reading in the unpruned data set and grabbing the markers that were not pruned
Genotypes_unpruned <- readRDS("nam_rqtl_geno_prob_array_james_pheno_subset.RDS")

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

pos_type = "gen_dist"

pheno_type <- "james"

#-----RUNNING THE SCRIPT-----#
source("qtl_confidence_intervals.R")

#-----STEP 2 END-----#


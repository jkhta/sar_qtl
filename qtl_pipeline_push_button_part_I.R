library(data.table)
library(parallel)
library(Matrix)
library(dplyr)
library(lme4)
library(lmerTest)
library(lme4qtl)
library(qtl)
library(pryr)
library(GridLMM)
library(ggplot2)

rm(list = ls())

#PART I

#directory with genotype probability array and blups
setwd("/home/jkta/projects/col_sha/complete_pipeline/")

#-----STEP 1 START-----#

#-----INPUT-----#
#phenotype file
nam_blups <- fread("nam_blups_combined_univariate.csv",
                   sep = ",",
                   header = TRUE,
                   stringsAsFactors = FALSE)

#genotype probability array
all_ril_geno_probs <- readRDS("nam_rqtl_geno_prob_array_comp_11_all_rils_NEW.RDS")

#-----PARAMETERS-----#
#output file name 
file_output_name <- "nam_rqtl_geno_prob_array_james_pheno_subset.RDS"

#maximum number of cores used to subset data
max_cores <- 5

#script for subsetting the genotype probability array by only genotypes that
#have phenotypes
source("geno_prob_array_subset_by_pheno.R")

#-----STEP 1 END-----#

#-----STEP 2 START-----#
rm(list = ls())

#-----INPUT-----#
#reading in the genotype probability array
Genotypes <- readRDS("nam_rqtl_geno_prob_array_james_pheno_subset.RDS")

#-----PARAMETERS-----#
#running a for loop that will calculate the corrleation of the marker with other markers
#and then drop the markers that are above a certain threshold with the one currently being tested
#threshold for dropping markers 
cor_threshold <- 0.99

#maximum number of cores to calculate marker correlations
max_cores <- 5

#generating the file names and then saving them as RDS files
file_output_name <- paste("nam_rqtl_geno_prob_array_james_", cor_threshold, "_", "pruned.rds", sep = "")

#script for pruning the markers so that the stepwise qtl method doesn't 
#include correlated markers
source("geno_prob_array_prune.R")

#-----STEP 2 END-----#

#-----STEP 3 START-----#
rm(list = ls())

#generating the list of phenotype, genotype, and proximal information 
#for the stepwise 

#-----INPUT-----#
#reading in the pruned genotypes, and will subset the complete genotypes later by the markers remaining in the pruned set
Genotypes_pruned <- readRDS("nam_rqtl_geno_prob_array_james_0.99_pruned.rds")

#READIN IN THE MARKERS FOR genotypes only found within the blups data set
#reading in the unpruned data set and grabbing the markers that were not pruned
Genotypes_unpruned <- readRDS("nam_rqtl_geno_prob_array_james_pheno_subset.RDS")

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

#-----PARAMETERS-----#
#maximum number of cores to run at once
max_cores <- 5

#threshold for proximal matrix based on correlations
cor_threshold <- 0.9

#threshold for proximal matrix based on genetic distance of testing marker
cm_threshold <- 10

#output name for file
file_output_name <- "nam_all_traits_ind_pop_pheno_geno_proximal_james.RDS"

#script for subsetting each populations phenotype and genotype data, and
#generating a matrix of markers to remove from the kinship matrix, to
#correct for proximal contamination
source("pheno_geno_proximal_lists.R")

#-----STEP 3 END-----#
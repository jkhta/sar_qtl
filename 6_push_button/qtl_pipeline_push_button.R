#this script will generate the genotype file subsetted by what we have phenotyped
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

#directory with genotype probability array and blups
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/6_push_button/")

#script for subsetting the genotype probability array by only genotypes that
#have phenotypes
source("geno_prob_array_subset_by_pheno.R")

rm(list = ls())

#script for pruning the markers so that the stepwise qtl method doesn't 
#include correlated markers
source("geno_prob_array_prune.R")

rm(list = ls())

#generating the list of phenotype, genotype, and proximal information 
#for the stepwise 

#-----PARAMETERS-----#
#maximum number of cores to run at once
max_cores <- 5

#threshold for proximal matrix based on correlations
cor_threshold <- 0.9

#threshold for proximal matrix based on genetic distance of testing marker
cm_threshold <- 10
#-----PARAMETERS-----#

source("pheno_geno_proximal_lists.R")

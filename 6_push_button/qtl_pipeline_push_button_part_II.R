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
#-----PARAMETERS-----#

#script for running permutations for each population/trait
source("ind_pop_permutations.R")

rm(list = ls())

#-----PARAMETERS-----#
#number of cores to run the permutations
max_cores <- 5

#run is a number that corresponds to a trait; goes from 1-8 for
#the 8 different traits; each number is run on a different cluster
run <- as.numeric(commandArgs(t = T)[1])

#-----PARAMETERS-----#
source("ind_trait_stepwise_QTL.R")

rm(list = ls())

#-----PARAMETERS-----#
#number of cores to run the permutations
max_cores <- 5

#run is a number that corresponds to a trait; goes from 1-8 for
#the 8 different traits; each number is run on a different cluster
run <- as.numeric(commandArgs(t = T)[1])

#-----PARAMETERS-----#

source("qtl_confidence_intervals.R")

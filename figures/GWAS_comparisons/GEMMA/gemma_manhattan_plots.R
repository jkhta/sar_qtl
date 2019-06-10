#this file will generate the GEMMA table for the supporting information
library(data.table)
library(ggplot2)
library(ggsci)
library(dplyr)
library(qqman)

rm(list = ls())

#HAVE TO MANUALLY SAVE PLOT
setwd("/Users/James/Documents/GitHub/sar_qtl/figures/GWAS_comparisons/GEMMA/")

#first the results without kinship matrix
gemma_bd_lm_ouput <- fread("nam_bd_geno_lm.assoc.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gemma_bd_lm_ouput$chr <- with(gemma_bd_lm_ouput, as.integer(sapply(strsplit(rs, split = "_"), function(x) x[1]))) 
gemma_bd_lm_ouput$ps <- with(gemma_bd_lm_ouput, as.integer(sapply(strsplit(rs, split = "_"), function(x) x[2]))) 

#manhattan plot of results for GEMMA linear model
manhattan(gemma_bd_lm_ouput, chr = "chr", bp = "ps", snp = "rs", p = "p_lrt", suggestiveline = FALSE, genomewideline = FALSE) 

#now plotting results using kinship matrix
gemma_bd_lmm_ouput <- fread("nam_bd_geno_lmm.assoc.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gemma_bd_lmm_ouput$chr <- with(gemma_bd_lmm_ouput, as.integer(sapply(strsplit(rs, split = "_"), function(x) x[1]))) 
gemma_bd_lmm_ouput$ps <- with(gemma_bd_lmm_ouput, as.integer(sapply(strsplit(rs, split = "_"), function(x) x[2]))) 

#manhattan plot of results for GEMMA mixed model
manhattan(gemma_bd_lmm_ouput, chr = "chr", bp = "ps", snp = "rs", p = "p_lrt", suggestiveline = FALSE, genomewideline = FALSE) 


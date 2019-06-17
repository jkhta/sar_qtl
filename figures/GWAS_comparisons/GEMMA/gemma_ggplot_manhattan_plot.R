#this file will generate the GEMMA table for the supporting information
library(data.table)
library(ggplot2)
library(ggsci)
library(dplyr)
library(qqman)

rm(list = ls())

setwd("/Users/James/Documents/GitHub/sar_qtl/figures/GWAS_comparisons/GEMMA/")

#first the results without kinship matrix
gemma_bd_geno_ouput <- fread("nam_bd_geno_lm.assoc.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gemma_bd_geno_ouput$chr <- with(gemma_bd_geno_ouput, as.integer(sapply(strsplit(rs, split = "_"), function(x) x[1]))) 
gemma_bd_geno_ouput$ps <- with(gemma_bd_geno_ouput, as.integer(sapply(strsplit(rs, split = "_"), function(x) x[2]))) 
head(gemma_bd_geno_ouput)

gemma_bd_geno_breaks <- gemma_bd_geno_ouput %>% group_by(chr) %>% summarize(center = (max(ps) + min(ps)) / 2)
gemma_bd_geno_breaks$chr <- as.integer(gemma_bd_geno_breaks$chr)
gemma_bd_geno_breaks$center <- as.numeric(gemma_bd_geno_breaks$center)

ggplot(gemma_bd_geno_ouput, aes(x = ps, y = -log10(p_lrt))) + 
    geom_point(aes(color = as.factor(chr))) +
    scale_color_npg() +
    facet_grid(~ chr, scales = "free") +
    scale_x_continuous(label = gemma_bd_geno_breaks$chr, breaks = gemma_bd_geno_breaks$center) +
    theme(panel.spacing = unit(0, "lines")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.title = element_blank(), legend.key = element_blank()) +
    theme(axis.text.x = element_blank()) +
    xlab("Chromosome") + 
    ylab("-log10(p)")

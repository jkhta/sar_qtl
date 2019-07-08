#this file will generate the GEMMA table for the supporting information
library(data.table)
library(ggplot2)
library(ggsci)
library(dplyr)
library(qqman)

rm(list = ls())

setwd("/Users/James/Documents/GitHub/sar_qtl/figures/GWAS_comparisons/GEMMA/")

#first the results without kinship matrix
gemma_bd_geno_lm_ouput <- fread("nam_bd_geno_lm.assoc.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gemma_bd_geno_lm_ouput$chr <- with(gemma_bd_geno_lm_ouput, as.integer(sapply(strsplit(rs, split = "_"), function(x) x[1]))) 
gemma_bd_geno_lm_ouput$ps <- with(gemma_bd_geno_lm_ouput, as.integer(sapply(strsplit(rs, split = "_"), function(x) x[2]))) 
head(gemma_bd_geno_lm_ouput)

gemma_bd_geno_lm_ouput_change <- gemma_bd_geno_lm_ouput %>% 
    
    # Compute chromosome size
    group_by(chr) %>% 
    summarise(chr_len=max(ps)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(gemma_bd_geno_lm_ouput, ., by=c("chr"="chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(chr, ps) %>%
    mutate( PScum=ps+tot)

gemma_bd_geno_breaks <- gemma_bd_geno_lm_ouput_change %>% group_by(chr) %>% summarize(center = (max(PScum) + min(PScum)) / 2)
gemma_bd_geno_breaks$chr <- as.integer(gemma_bd_geno_breaks$chr)
gemma_bd_geno_breaks$center <- as.numeric(gemma_bd_geno_breaks$center)

ggplot(gemma_bd_geno_lm_ouput_change, aes(x = PScum, y = -log10(p_lrt))) + 
    geom_point(aes(color = as.factor(chr))) +
    scale_color_npg() +
    #facet_grid(~ chr, scales = "free") +
    scale_x_continuous(label = gemma_bd_geno_breaks$chr, breaks = gemma_bd_geno_breaks$center) +
    theme(panel.spacing = unit(0, "lines")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.position = "none") +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20)) +
    xlab("Chromosome") + 
    ylab("-log10(p)")
ggsave("gemma_bd_geno_lm_manhattan_plot.png", device = "png", width = 10, height = 8)

#the results with the kinship matrix
gemma_bd_geno_lmm_ouput <- fread("nam_bd_geno_lmm.assoc.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gemma_bd_geno_lmm_ouput$chr <- with(gemma_bd_geno_lmm_ouput, as.integer(sapply(strsplit(rs, split = "_"), function(x) x[1]))) 
gemma_bd_geno_lmm_ouput$ps <- with(gemma_bd_geno_lmm_ouput, as.integer(sapply(strsplit(rs, split = "_"), function(x) x[2]))) 
head(gemma_bd_geno_lmm_ouput)

gemma_bd_geno_lmm_ouput_change <- gemma_bd_geno_lmm_ouput %>% 
    
    # Compute chromosome size
    group_by(chr) %>% 
    summarise(chr_len=max(ps)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(gemma_bd_geno_lmm_ouput, ., by=c("chr"="chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(chr, ps) %>%
    mutate(PScum = ps+tot)

gemma_bd_geno_breaks <- gemma_bd_geno_lmm_ouput_change %>% group_by(chr) %>% summarize(center = (max(PScum) + min(PScum)) / 2)
gemma_bd_geno_breaks$chr <- as.integer(gemma_bd_geno_breaks$chr)
gemma_bd_geno_breaks$center <- as.numeric(gemma_bd_geno_breaks$center)

ggplot(gemma_bd_geno_lmm_ouput_change, aes(x = PScum, y = -log10(p_lrt))) + 
    geom_point(aes(color = as.factor(chr))) +
    scale_color_npg() +
    #facet_grid(~ chr, scales = "free") +
    scale_x_continuous(label = gemma_bd_geno_breaks$chr, breaks = gemma_bd_geno_breaks$center) +
    theme(panel.spacing = unit(0, "lines")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.position = "none") +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20)) +
    xlab("Chromosome") + 
    ylab("-log10(p)")
ggsave("gemma_bd_geno_lmm_manhattan_plot.png", device = "png", width = 10, height = 8)

#bd_gxe results with kinship matrix
gemma_bd_gxe_ouput <- fread("nam_bd_gxe_lmm.assoc.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gemma_bd_gxe_ouput$chr <- with(gemma_bd_gxe_ouput, as.integer(sapply(strsplit(rs, split = "_"), function(x) x[1]))) 
gemma_bd_gxe_ouput$ps <- with(gemma_bd_gxe_ouput, as.integer(sapply(strsplit(rs, split = "_"), function(x) x[2]))) 
head(gemma_bd_gxe_ouput)

gemma_bd_gxe_ouput_change <- gemma_bd_gxe_ouput %>% 
    
    # Compute chromosome size
    group_by(chr) %>% 
    summarise(chr_len=max(ps)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(gemma_bd_gxe_ouput, ., by=c("chr"="chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(chr, ps) %>%
    mutate( PScum=ps+tot)

gemma_bd_gxe_breaks <- gemma_bd_gxe_ouput_change %>% group_by(chr) %>% summarize(center = (max(PScum) + min(PScum)) / 2)
gemma_bd_gxe_breaks$chr <- as.integer(gemma_bd_gxe_breaks$chr)
gemma_bd_gxe_breaks$center <- as.numeric(gemma_bd_gxe_breaks$center)

ggplot(gemma_bd_gxe_ouput_change, aes(x = PScum, y = -log10(p_lrt))) + 
    geom_point(aes(color = as.factor(chr))) +
    scale_color_npg() +
    #facet_grid(~ chr, scales = "free") +
    scale_x_continuous(label = gemma_bd_gxe_breaks$chr, breaks = gemma_bd_geno_breaks$center) +
    theme(panel.spacing = unit(0, "lines")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.position = "none") +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20)) +
    xlab("Chromosome") + 
    ylab("-log10(p)")
ggsave("gemma_bd_gxe_lmm_manhattan_plot.png", device = "png", width = 10, height = 8)


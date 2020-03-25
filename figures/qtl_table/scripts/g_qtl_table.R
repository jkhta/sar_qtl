#this script will make a table for the G QTL found
library(data.table)
library(xtable)
library(plyr)

rm(list = ls())

#reading in tables for QTL found for genotype random effects (G QTL)
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/qtl_table/data/")

g_qtl_list <- lapply(list.files(pattern = "geno_qtl_ci.csv"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE))
g_qtl_phenotypes <- sapply(strsplit(list.files(pattern = "geno_qtl_ci.csv"), split = "_geno"), function(x) x[1])

#adding trait names to data and then list -> data table
for (i in 1:length(g_qtl_list)) {
  g_qtl_list[[i]]$trait <- g_qtl_phenotypes[i]
}

g_qtl_comb <- rbindlist(g_qtl_list)

#changing trait names and adding chromosome columns
g_qtl_comb$trait <- mapvalues(g_qtl_comb$trait, from = c("bd", "h3_h1", "i_dry", "r_dry"), to = c("BD", "IG", "IB", "RB"))
g_qtl_comb$chromosome <- sapply(strsplit(g_qtl_comb$qtl, split = "_"), function(x) x[2])

g_qtl_comb_reordered <- subset(g_qtl_comb, select = c(trait, qtl, avg_snp_pve, score, chromosome, left_bound, right_bound))
colnames(g_qtl_comb_reordered) <- c("Trait", "QTL Marker", "SNP PVE", "-log10p", "Chromosome","Left Bound", "Right Bound")

g_qtl_comb_reordered_complete <- c()

#for each trait subset the QTL for that trait
for (i in unique(g_qtl_comb_reordered$Trait)) {
  g_qtl_comb_reordered_subset <- subset(g_qtl_comb_reordered, Trait == i)
  
  #for each chromosome with QTL, subset the data by chromosome and create new QTL names
  for (j in unique(g_qtl_comb_reordered_subset$Chromosome)) {
    g_qtl_comb_reordered_chr_subset <- subset(g_qtl_comb_reordered_subset, Chromosome == j)
    g_qtl_comb_reordered_chr_subset$Base <- as.numeric(sapply(strsplit(g_qtl_comb_reordered_chr_subset$`QTL Marker`, split = "_"), function(x) x[3]))
    g_qtl_comb_reordered_chr_subset <- g_qtl_comb_reordered_chr_subset[order(g_qtl_comb_reordered_chr_subset$Base)]
    g_qtl_comb_reordered_chr_subset$QTL_no <- 1:nrow(g_qtl_comb_reordered_chr_subset)
    g_qtl_comb_reordered_chr_subset$QTL <- paste(paste(g_qtl_comb_reordered_chr_subset$Trait, g_qtl_comb_reordered_chr_subset$Chromosome, sep = ""), g_qtl_comb_reordered_chr_subset$QTL_no, sep = "_")
    g_qtl_comb_reordered_complete <- rbindlist(list(g_qtl_comb_reordered_complete, g_qtl_comb_reordered_chr_subset))
  }
}

#removing unnecessary columns
g_qtl_comb_reordered_complete$QTL_no <- NULL
g_qtl_comb_reordered_complete$`-log10p` <- NULL

#reordering variables again
g_qtl_comb_reordered_complete <- subset(g_qtl_comb_reordered_complete, select = c(Trait, QTL, `SNP PVE`, `QTL Marker`, Chromosome, `Left Bound`, `Right Bound`))

print(xtable(g_qtl_comb_reordered_complete), include.rownames = FALSE)

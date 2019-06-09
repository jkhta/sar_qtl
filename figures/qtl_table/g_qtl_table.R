#this script will make a table for the G QTL found
library(data.table)
library(xtable)

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/qtl_table/")

g_qtl_list <- lapply(list.files(pattern = "geno_qtl_ci.csv"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE))
g_qtl_phenotypes <- sapply(strsplit(list.files(pattern = "geno_qtl_ci.csv"), split = "_geno"), function(x) x[1])

for (i in 1:length(g_qtl_list)) {
  g_qtl_list[[i]]$trait <- g_qtl_phenotypes[i]
}

g_qtl_comb <- rbindlist(g_qtl_list)

g_qtl_comb_reordered <- subset(g_qtl_comb, select = c(trait, qtl, avg_snp_pve, score, left_bound, right_bound))
colnames(g_qtl_comb_reordered) <- c("Trait", "QTL", "SNP PVE", "-log10p", "Left Bound", "Right Bound")

print(xtable(g_qtl_comb_reordered), include.rownames = FALSE)

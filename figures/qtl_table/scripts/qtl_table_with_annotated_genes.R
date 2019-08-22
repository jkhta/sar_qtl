#this script will try to make the annotated confidence interval tables look a bit better
library(data.table)
library(plyr)
library(xtable)

rm(list = ls())

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/gene_annotation/data/")

sar_gene_counts <- fread("sar_gxe_gene_counts.csv", 
                         sep = ",", 
                         header = TRUE,
                         stringsAsFactors = FALSE)

colnames(sar_gene_counts) <- c("trait", "qtl", "num_genes")

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/qtl_table/data/")
sar_gxe_qtl_list <- lapply(list.files(pattern = "gxe_qtl_ci.csv"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE))
sar_gxe_qtl_pheno <- sapply(strsplit(list.files(pattern = "gxe_qtl_ci.csv"), split = "_qtl_ci"), function(x) x[1])

#labeling qtl tables with trait names
for (i in 1:length(sar_gxe_qtl_list)) {
  sar_gxe_qtl_list[[i]]$trait <- sar_gxe_qtl_pheno[i]
}

sar_gxe_qtl_dt <- rbindlist(sar_gxe_qtl_list)
sar_gene_counts$trait <- paste(sar_gene_counts$trait, "gxe", sep = "_")
sar_gxe_qtl_dt_with_anno <- merge(sar_gxe_qtl_dt, sar_gene_counts, by = c("trait", "qtl"))
sar_gxe_qtl_dt_with_anno_subset <- subset(sar_gxe_qtl_dt_with_anno, select = c(trait, qtl, avg_snp_pve, score, left_bound, right_bound, num_genes))
sar_gxe_qtl_dt_with_anno_subset[is.na(sar_gxe_qtl_dt_with_anno_subset$num_genes), ]$num_genes  <- 0
colnames(sar_gxe_qtl_dt_with_anno_subset) <- c("Trait", "QTL", "SNP PVE", "-log10(p)", "Left Bound", "Right Bound", "# Annotated Genes")

#multiple SNP PVE by 100 to represent a percentage
sar_gxe_qtl_dt_with_anno_subset$`SNP PVE` <- sar_gxe_qtl_dt_with_anno_subset$`SNP PVE` * 100

#removing the underscore and gxe
sar_gxe_qtl_dt_with_anno_subset$Trait <- gsub("_|gxe", "", sar_gxe_qtl_dt_with_anno_subset$Trait)

#changing trait name
sar_gxe_qtl_dt_with_anno_subset$Trait <- revalue(sar_gxe_qtl_dt_with_anno_subset$Trait, c("bd" = "BD", "h3h1" = "IG", "idry" = "IB", "rdry" = "RB"))

#creating new QTL names
sar_gxe_qtl_dt_with_anno_subset$`QTL Marker` <- sar_gxe_qtl_dt_with_anno_subset$QTL
sar_gxe_qtl_dt_with_anno_complete <- c()

for (i in unique(sar_gxe_qtl_dt_with_anno_subset$Trait)) {
  sar_gxe_qtl_dt_with_anno_trait_subset <- subset(sar_gxe_qtl_dt_with_anno_subset, Trait == i)
  sar_gxe_qtl_dt_with_anno_trait_subset$QTL_no <- 1:nrow(sar_gxe_qtl_dt_with_anno_trait_subset)
  sar_gxe_qtl_dt_with_anno_trait_subset$QTL <- paste(sar_gxe_qtl_dt_with_anno_trait_subset$Trait, sar_gxe_qtl_dt_with_anno_trait_subset$QTL_no, sep = "_")
  sar_gxe_qtl_dt_with_anno_complete <- rbindlist(list(sar_gxe_qtl_dt_with_anno_complete, sar_gxe_qtl_dt_with_anno_trait_subset))
}

#removing and rearranging columns
sar_gxe_qtl_dt_with_anno_complete$QTL_no <- NULL
sar_gxe_qtl_dt_with_anno_complete <- subset(sar_gxe_qtl_dt_with_anno_complete, select = c(Trait, QTL, `SNP PVE`, `-log10(p)`, `QTL Marker`, `Left Bound`, `Right Bound`, `# Annotated Genes`))

#writing data
fwrite(sar_gxe_qtl_dt_with_anno_complete, "sar_gxe_qtl_info_with_anno.csv", sep = ",", row.names = FALSE)
fwrite(subset(sar_gxe_qtl_dt_with_anno_complete, select = c(Trait, `QTL Marker`, QTL)), "sar_gxe_matched_qtl_name.csv", sep = ",", row.names = FALSE)

print(xtable(sar_gxe_qtl_dt_with_anno_complete, digits = 2), include.rownames = FALSE)

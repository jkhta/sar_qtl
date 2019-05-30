#this script will try to make the annotated confidence interval tables look a bit better
library(data.table)
library(plyr)

rm(list = ls())

setwd("/Users/jkhta/Desktop/nam_cam_fixing/32 - gene_annotation/output/")

#reading all of the annotated ci tables
annotated_ci_phenotypes <- sapply(strsplit(list.files(), split = "_"), function(x) paste(x[1], x[2], sep = "_"))
annotated_ci_tables <- lapply(list.files(), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE))

#putting a trait column for each ci table for each trait
for (i in 1:length(annotated_ci_tables)) {
  annotated_ci_tables[[i]]$trait <- annotated_ci_phenotypes[i]
}

#combining the annotated ci tables
annotated_ci_tables_comb <- rbindlist(annotated_ci_tables)
annotated_ci_tables_comb$qtl_id <- with(annotated_ci_tables_comb, paste(trait, qtl, sep = "_"))

#making a separate data frame of all QTL, even ones without any annotated genes
qtl_id_list <- data.frame(qtl_id = unique(annotated_ci_tables_comb$qtl_id))

#removing the QTL that don't have annotated genes
annotated_ci_tables_comb_subset <- subset(annotated_ci_tables_comb, symbol != NA | symbol != "")

#making sure which qtl have 0 annotations, and which to change the frequency
#from 1 to 0 later when using the count function
annotated_ci_tables_comb[is.na(annotated_ci_tables_comb$Alias), ]

#maybe what would be better would be a count of genes annotated
annotated_ci_table_count <- count(annotated_ci_tables_comb_subset, vars = "qtl_id")
annotated_ci_table_count_comp <- merge(annotated_ci_table_count, qtl_id_list, by = "qtl_id", all = TRUE)

#changing the name for later merging
colnames(annotated_ci_table_count) <- c("qtl_id", "annotated_gene_count")

#picking the columns for merging
names(annotated_ci_tables_comb)
annotated_ci_tables_comb_subset <- subset(annotated_ci_tables_comb, select = c(qtl_id, qtl, trait, chr, left_bound, right_bound, trait))

#only getting the unique values otherwise the merge will be a mess
annotated_ci_tables_comb_unique <- unique(annotated_ci_tables_comb_subset)

#merging
ci_table_gene_count_comb <- merge(annotated_ci_table_count, annotated_ci_tables_comb_unique, by = c("qtl_id"), all = TRUE)
ci_table_gene_count_comb$qtl_id <- NULL

ci_table_gene_count_comb_reorder <- subset(ci_table_gene_count_comb, select = c(qtl, trait, annotated_gene_count, chr, left_bound, right_bound))
ci_table_gene_count_comb_reorder$chr <- sapply(strsplit(ci_table_gene_count_comb_reorder$qtl, split = "_"), function(x) paste("Chr", x[2], sep = ""))

setwd("/Users/jkhta/Desktop/nam_cam_fixing/32 - gene_annotation/output/")
fwrite(ci_table_gene_count_comb_reorder, 
       file = "qtl_gene_annotation_count.csv", 
       sep = ",", 
       row.names = FALSE, 
       col.names = TRUE, 
       na = NA)

geno_ci_table_gene_count_comb_reorder <- subset(ci_table_gene_count_comb_reorder, grepl("geno", ci_table_gene_count_comb_reorder$trait))
gxe_ci_table_gene_count_comb_reorder <- subset(ci_table_gene_count_comb_reorder, grepl("gxe", ci_table_gene_count_comb_reorder$trait))

#only writing for geno or gxe
fwrite(geno_ci_table_gene_count_comb_reorder,
       file = "geno_qtl_gene_annotation_count.csv", 
       sep = ",", 
       row.names = FALSE, 
       col.names = TRUE, 
       na = NA)

fwrite(gxe_ci_table_gene_count_comb_reorder,
       file = "gxe_qtl_gene_annotation_count.csv", 
       sep = ",", 
       row.names = FALSE, 
       col.names = TRUE, 
       na = NA)

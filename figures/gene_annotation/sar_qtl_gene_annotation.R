#this script will try to make the annotated confidence interval tables look a bit better
library(data.table)
library(plyr)
library(xtable)

rm(list = ls())

setwd("/Users/James/Documents/GitHub/sar_qtl/figures/gene_annotation/")

#reading all of the annotated ci tables
annotated_ci_phenotypes <- sapply(strsplit(list.files(pattern = "gxe_qtl_ci"), split = "_"), function(x) paste(x[1], x[2], sep = "_"))
annotated_ci_tables <- lapply(list.files(pattern = "gxe_qtl_ci"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE))

#putting a trait column for each ci table for each trait
for (i in 1:length(annotated_ci_tables)) {
    annotated_ci_tables[[i]]$trait <- annotated_ci_phenotypes[i]
}

annotated_ci_comp <- rbindlist(annotated_ci_tables)

annotated_ci_comp_subset <- subset(annotated_ci_comp, select = c("trait", "qtl", "ID", "chr", "start", "end", "symbol", "Note"))

fwrite(annotated_ci_comp_subset, file = "sar_gxe_gene_annotation.csv", sep = ",", row.names = FALSE)

#now to generate sar gene annotation count for a smaller table
annotated_ci_comp_subset$trait_qtl <- with(annotated_ci_comp_subset, paste(trait, qtl, sep = "_"))

annotated_ci_comp_count <- count(annotated_ci_comp_subset, vars = c("trait_qtl"))
annotated_ci_comp_count$trait <- sapply(strsplit(annotated_ci_comp_count$trait_qtl, split = "_m"), function(x) x[1])
annotated_ci_comp_count$qtl <- sapply(strsplit(annotated_ci_comp_count$trait_qtl, split = "_m"), function(x) paste("m", x[2], sep = ""))

annotated_ci_comp_count <- subset(annotated_ci_comp_count, select = c("trait", "qtl", "freq"))
colnames(annotated_ci_comp_count) <- c("Trait", "QTL", "Freq")
annotated_ci_comp_count$Trait <- gsub("_gxe", "", annotated_ci_comp_count$Trait)

fwrite(annotated_ci_comp_count, file = "sar_gxe_gene_counts.csv", sep = ",", row.names = FALSE)
nrow(annotated_ci_comp_count)

annotated_ci_comp_count_1 <- annotated_ci_comp_count[1:10, ]
annotated_ci_comp_count_2 <- annotated_ci_comp_count[11:20, ]
annotated_ci_comp_count_comb <- cbind(annotated_ci_comp_count[1:10, ], annotated_ci_comp_count[11:20, ])

print(xtable(annotated_ci_comp_count_comb, caption = "Table 3"), include.rownames = FALSE)



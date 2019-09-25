#this script will attempt to merge the tair files and the shade avoidance response files to
#overlay plots with genes

#first need to take kazu's gene list and get the unique genes from the 1hr and the 4hr treatment
library(data.table)

rm(list = ls())

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/gene_annotation/data/raw_gene_list_and_tair/")

#reading in the petiole DE genes separately, and then going to combine them into one
#list
kazu_h1 <- fread("kazu_shade_avoidance_gene_list_1hr.csv",
                 sep = ",",
                 header = TRUE,
                 stringsAsFactors = FALSE)
colnames(kazu_h1)[1] <- "ID"

kazu_h4 <- fread("kazu_shade_avoidance_gene_list_4hr.csv",
                 sep = ",",
                 header = TRUE,
                 stringsAsFactors = FALSE)
colnames(kazu_h4)[1] <- "ID"

#reading in the leaf/apical shade responsive genes from Kazu's 2015 supplemental information
kazu_leaf_SAR_genes <- fread("kazu_leaf_apical_shade_responsive_genes.csv",
                             sep = ",",
                             header = TRUE,
                             stringsAsFactors = FALSE)
colnames(kazu_leaf_SAR_genes)[1] <- "ID"

#reading in the core genes list from the romano sellaro study
romano_core_SAR_genes <- fread("romano_core_genes.csv",
                               sep = ",",
                               header = TRUE,
                               stringsAsFactors = FALSE)
colnames(romano_core_SAR_genes)[1] <- "ID"

#grabbing the unique values from each of the lists in order to compile a list of unique genes
kazu_comb_genes <- rbind(kazu_h1, kazu_h4)


#checking what is the max p-value in kazu's data
max(kazu_comb_genes$table.PValue, na.rm = TRUE)

#so all of the genes are significant; now getting the unique list of genes
#so there are 239 unique genes in the set
kazu_unique_genes <- unique(c(kazu_comb_genes$ID, kazu_leaf_SAR_genes$ID, romano_core_SAR_genes$ID))

#reading in the other gene lists that i manually curated and from kazu's mutant and phenotype list
kazu_manual_genes <- fread("manually_curated_kazu_2015_plos_genes.csv",
                           sep = ',',
                           header = TRUE, 
                           stringsAsFactors = FALSE)

head(kazu_manual_genes)

kazu_mutant_genes <- fread("kazu_mutant_gene_list.csv",
                           sep = ',',
                           header = TRUE, 
                           stringsAsFactors = FALSE)

head(kazu_mutant_genes)

#putting all of the gene lists into one file
kazu_gene_list <- list(kazu_h1 = kazu_h1,
                       kazu_h4 = kazu_h4,
                       kazu_manual = kazu_manual_genes,
                       kazu_mutant = kazu_mutant_genes,
                       kazu_unique_genes = data.frame(ID = kazu_unique_genes))


#now need to subset all of the lists so that i am just getting the id and symbol
kazu_gene_list_subset <- lapply(kazu_gene_list, function(x) subset(x, select = ID))

#combining all of the kazu's genes together in order to read it in with tair
kazu_genes_comb <- rbindlist(kazu_gene_list_subset)

#grabbing only the unique genes
kazu_gene_comb_unique <- as.character(unique(kazu_genes_comb$ID))

#now need to read in the araport 11 data in order to get more information and combine
#the chromsome position with the listed gene
araport_11_first_lines <- readLines("Araport11_GFF3_genes_transposons.201606.gff",
                                    n = 50)

#reading in the araport data; using araport 11 as the base instead of tair 10 because it is more recent
tair_info <- as.data.frame(data.table(read.delim("Araport11_GFF3_genes_transposons.201606.gff", 
                                                 header = F, 
                                                 comment.char = "#"), 
                                      keep.rownames = TRUE))


#changing the tair infor into a character for easier parsing through;
#also changing the names of the columns so that they are easier to subset
tair_info$V9 <- as.character(tair_info$V9)
names(tair_info) <- c("feature", "chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

tair_genes_mrna <- data.table(subset(tair_info, grepl("gene|mRNA", type, ignore.case = TRUE)))

gene_mrna_annotation_extraction <- function(tair_data) {
  single_gene_information <- tair_data["attributes"]
  single_gene_row <- tair_data["feature"]
  single_gene_information_strsplit <- strsplit(single_gene_information, split = ";")
  
  #taking the wanted gene information (gene name, etc.)
  single_gene_information_subset <- lapply(single_gene_information_strsplit, function(x) x[grepl("ID=|Note=|symbol=|Alias=|full_name=|curator_summary=", x, ignore.case = TRUE)])
  
  #now splitting the gene information again by = in order to generate column names and gene information 
  single_gene_information_col_info_strsplit <- lapply(single_gene_information_subset, function(x) strsplit(x, split = "="))
  
  #now generating a data frame with the needed information (row number, column name (variable), and the information from the single gene information strsplit)
  single_gene_df <- data.frame(feature = single_gene_row, 
                               col_name = unlist(sapply(single_gene_information_col_info_strsplit, function(x) lapply(x, function(y) y[1]))),
                               info = unlist(sapply(single_gene_information_col_info_strsplit, function(x) lapply(x, function(y) y[2]))), 
                               row.names = NULL)
}

#using the apply function to see if that speeds up grabbing the information
#this speeds things up tremendously; it takes about 1 minutes vs. 15-20 minutes for the for loop
start_time <- Sys.time()
gene_mrna_annotation <- apply(tair_genes_mrna, 1, function(x) gene_mrna_annotation_extraction(x))
gene_mrna_annotation_dt <- rbindlist(gene_mrna_annotation)
end_time <- Sys.time()
print(end_time - start_time)

#changing from long format to wide format
gene_mrna_annotation_dt_wide <- dcast(gene_mrna_annotation_dt, feature ~ col_name, value.var = "info")
names(gene_mrna_annotation_dt_wide)[1] <- "feature"

#merging with the subsetted tair data frame
tair_annotated_dt <- merge(tair_genes_mrna, gene_mrna_annotation_dt_wide, by = "feature")

tair_annotated_dt$chr_num <- gsub("chr", "", tair_annotated_dt$chr, ignore.case = TRUE)

tair_annotation <- subset(tair_annotated_dt, chr_num %in% as.character(1:5))
tair_annotation$chr_num <- as.numeric(tair_annotation$chr_num)

kazu_unique_genes_df <- data.frame(ID = kazu_gene_comb_unique)

kazu_tair_merge <- merge(kazu_unique_genes_df, tair_annotation, by = "ID")

head(kazu_tair_merge)

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/gene_annotation/data/")

fwrite(kazu_tair_merge, 
       file = "kazu_romano_tair_merge.csv", 
       sep = ",", 
       row.names = FALSE)
fwrite(tair_annotation, 
       file = "araport11_tair_annotation_parsed.csv",
       sep = ",",
       row.names = FALSE)

#now after merging kazu's gene list with the tair data in order to get the chromosome information
#and the base pair information, now need to find out a way that will let me plot out the chromsomes
#and the genes in kazu's list, and then plot out my qtls and their confidence intervals
kazu_tair_genes <- fread("kazu_tair_merge.csv",
                         sep = ",",
                         header = TRUE,
                         stringsAsFactors = FALSE)

head(kazu_tair_genes)

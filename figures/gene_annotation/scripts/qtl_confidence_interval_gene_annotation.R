#this script will annotate the confidence intervals with genes
library(data.table)

rm(list = ls())

#reading in the merged gene list from kazu's 2015 paper; merged with script
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/gene_annotation/data/")
kazu_shade_genes <- fread("kazu_romano_tair_merge.csv",
                          sep = ",",
                          header = TRUE,
                          stringsAsFactors = FALSE)

for (j in list.files(pattern = "qtl_ci.csv")) {
  #grabbing the phenotype name
  pheno_name <- sapply(strsplit(j, split = "_qtl_ci"), function(x) x[1])
  
  #need to put this here again because at the end of the for loop i change the directory
  setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/qtl_table/data/")
  #reading in an example data frame for confidence intervals
  qtl_confidence_intervals <- fread(j, 
                                    sep = ",", 
                                    header = TRUE, 
                                    stringsAsFactors = FALSE)
  
  #information for the confidence intervals
  qtl_confidence_intervals$chr <- sapply(strsplit(qtl_confidence_intervals$qtl, split = "_"), function(x) as.numeric(x[2]))
  qtl_confidence_intervals$left_pos <- sapply(strsplit(qtl_confidence_intervals$left_bound, split = "_"), function(x) as.numeric(x[3]))
  qtl_confidence_intervals$right_pos <- sapply(strsplit(qtl_confidence_intervals$right_bound, split = "_"), function(x) as.numeric(x[3]))
  
  
  trait_qtl_annotation <- list()
  
  #now need to loop over each qtl, and then annotate the confidence intervals with known QTL
  for (i in 1:nrow(qtl_confidence_intervals)) {
    #grabbing the qtl chromosome so that i can subset kazu's gene list by chromosome
    qtl_chr <- qtl_confidence_intervals[i, ]
    
    #subsetting kazu's gene list by the chromosome
    kazu_chr_subset <- subset(kazu_shade_genes, chr_num == qtl_chr$chr)
    
    #now need to see what genes are within the bounds
    kazu_chr_subset_wi_bounds <- subset(kazu_chr_subset, start > qtl_chr$left_pos & start < qtl_chr$right_pos)
    
    if (nrow(kazu_chr_subset_wi_bounds) == 0) {
      #filling the list with qtl
      kazu_chr_subset_wi_bounds_w_qtl_info <- cbind(t(as.data.frame(rep(NA, times = ncol(kazu_chr_subset_wi_bounds)))), qtl_chr)
      colnames(kazu_chr_subset_wi_bounds_w_qtl_info)[1:ncol(kazu_chr_subset_wi_bounds)] <- colnames(kazu_chr_subset_wi_bounds)
      trait_qtl_annotation[[i]] <- kazu_chr_subset_wi_bounds_w_qtl_info
    } else {
      #need the name of the qtl so that i know what annotations belong to what
      kazu_chr_subset_wi_bounds_w_qtl_info <- cbind(kazu_chr_subset_wi_bounds, qtl_chr)
      
      #filling the list with qtl
      trait_qtl_annotation[[i]] <- kazu_chr_subset_wi_bounds_w_qtl_info
    }
    
    #combining all of the qtl genes together
    trait_qtl_annotation_df <- rbindlist(trait_qtl_annotation, fill = TRUE)
    
    #generating the file name for the confidence intervals
    file_name <- paste(pheno_name, "annotated_qtl_ci.csv", sep = "_")
    
    #writing a file for the qtl confidence intervals 
    setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/gene_annotation/data/")
    fwrite(trait_qtl_annotation_df, 
           file = file_name, 
           sep = ",", 
           row.names = FALSE, 
           col.names = TRUE, na = "NA")
  }
}


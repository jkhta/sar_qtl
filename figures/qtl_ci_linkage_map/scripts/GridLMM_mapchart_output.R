#this file will generate a mapchart output using the GridLMM output
library(qtl)
library(data.table)
library(plyr)

rm(list = ls())

#reading in the nam marker information
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/6_push_button/data")
nam_markers <- fread("nam_marker_info_final.csv",
                 sep = ",",
                 header = TRUE,
                 stringsAsFactors = FALSE)
colnames(nam_markers) <- c("marker", "chr", "pos")
nam_markers$marker <- paste("m", nam_markers$marker, sep = "_")

trait_type <- "covariate"

#reading in an example trait 95% CI using lme4qtl
setwd(paste("/Users/jkhta/Documents/GitHub/sar_qtl/figures/qtl_ci_linkage_map/", trait_type, "/", sep = ""))

phenotype_ci_list <- lapply(list.files(pattern = "gxe"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE))
phenotype_names <- sapply(strsplit(list.files(pattern = "gxe"), split = "_gxe"), function(x) x[1])
phenotype_names <- gsub("_", "", phenotype_names)
names(phenotype_ci_list) <- phenotype_names

for (i in 1:length(phenotype_ci_list)) {
  single_phenotype <- phenotype_ci_list[[i]]
  single_phenotype$phenotype_name <- phenotype_names[i]
  single_phenotype$color <- i
  single_phenotype$qtl_no <- 1:nrow(single_phenotype)
  phenotype_ci_list[[i]] <- single_phenotype
}

phenotype_ci_dt <- rbindlist(phenotype_ci_list)
phenotype_ci_dt$chr <- sapply(strsplit(phenotype_ci_dt$qtl, split = "_"), function(x) x[2])

#replacing trait names 
phenotype_ci_dt$phenotype_name <- revalue(phenotype_ci_dt$phenotype_name, c("bd" = "BD_SAR", "h3h1" = "IG_SAR", "idry" = "IB_SAR", "rdry" = "RB_SAR"))

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/qtl_ci_linkage_map/img/")

if (trait_type == "no_covariate") {
  sink("GridLMM_gxe_no_cov_cis.mct")
} else if (trait_type == "covariate") {
  sink("GridLMM_gxe_cov_cis.mct")
}

for (i in 1:5) {
  chr_markers <- subset(nam_markers, chr == i)
  
  cat(sprintf("group %s \n \n", i))
  
  for (j in 1:nrow(chr_markers)) {
    #grabbing an individual marker at time
    chr_single_marker <- chr_markers[j, ]
    chr_single_marker_id <- chr_single_marker$marker
    chr_single_marker_pos <- chr_single_marker$pos
    cat(sprintf('"" %s \n', chr_single_marker_pos))
  }
  
  chr_qtls <- subset(phenotype_ci_dt, chr == i)
  
  if (nrow(chr_qtls) > 0) {
    cat("\n")
    cat("qtls ;\n")
    
    #adding individual qtl information
    #grab each qtl at a time and send it to the text file
    for (j in 1:nrow(chr_qtls)) {
      #grabbing all of the qtl information needed
      chr_single_qtl <- chr_qtls[j, ]
      chr_single_phenotype_name <- chr_single_qtl$phenotype_name
      chr_single_qtl_start_pos <- unlist(subset(chr_markers, marker == chr_single_qtl$left_bound, select = pos))
      chr_single_qtl_end_pos <- unlist(subset(chr_markers, marker == chr_single_qtl$right_bound, select = pos))
      chr_single_qtl_color <- chr_single_qtl$color
      
      #now need to print all of the information for the qtl
      #the order is qtl name, start start start end color 
      cat(sprintf("%s %s %s %s %s B C%s\n \n", 
                  chr_single_phenotype_name, 
                  chr_single_qtl_start_pos,
                  chr_single_qtl_start_pos,
                  chr_single_qtl_start_pos,
                  chr_single_qtl_end_pos,
                  chr_single_qtl_color))
    }
  } else if (nrow(chr_qtls) < 1) {
    next
  }
}
sink()



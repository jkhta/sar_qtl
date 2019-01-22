#getting the blups from the univariate models and then mapping qtl with those BLUPs
library(brms)
library(data.table)
library(plyr)

rm(list = ls())

#setting the working directory
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/2_lmm_models/output/brms_models/")

traits <- c("bd", "h3_h1", "i_dry", "r_dry")
pop_names <- c("blh_col", "bur_col", "cvi_col", "ita_col", "jea_col", "oy_col", "sha_col")
pheno_names <- paste(rep(traits, each = 2), rep(c("geno", "gxe"), times = length(traits)), sep = "_")

#function to get the random effect estimates (BLUPs) for each trait and population
pop_ranefs <- function(trait_brm_list) {
  #getting the list of random effects
  pop_ranef_list <- lapply(pop_trait_brm_list, function(x) ranef(x))
  
  #changing the names of the random effects
  pop_ranef_df_list <- lapply(pop_ranef_list, function(x) 
    do.call(cbind, lapply(alply(x$geno, 3), function(y) data.frame(y[, 1]))))
  
  #replacing all of the column names for each df in the list
  for (i in 1:length(pop_ranef_df_list)) {
    #grabbing an individual ranef df
    pop_ranef_df_single <- data.table(pop_ranef_df_list[[i]], keep.rownames = TRUE)
    
    pheno_names_subset <- pheno_names[grepl(traits[i], pheno_names)]
    
    #renaming each of the columns in the data frame
    colnames(pop_ranef_df_single) <- c("geno", pheno_names_subset)
    
    #replacing the df's in the list 
    pop_ranef_df_list[[i]] <- pop_ranef_df_single
  }
  
  pop_trait_list <- Reduce(function(x, y) merge(x, y, by = "geno", all = TRUE), pop_ranef_df_list)
  return(pop_trait_list)
}

#iterates through each population and then generates a data frame with population blups
#population blups are appended sequentially
all_pop_ranef_dfs <- c()

#for each population
for (i in pop_names) {
  #grabbing only the files for a population
  pop_files <- list.files(pattern = i)
  
  #grabbing only the traits we are analyzing in the path model
  pop_trait_files <- pop_files[grepl(paste(traits, collapse = "|"), pop_files)]
  
  #reading in the brms files
  pop_trait_brm_list <- lapply(pop_trait_files, function(x) readRDS(x))
  
  #extracting out the blups for the population and formatting it into a data table
  pop_ranefs_df <- pop_ranefs(pop_trait_brm_list)
  
  #combining the population random effects together into a single data frame
  all_pop_ranef_dfs <- rbindlist(list(all_pop_ranef_dfs, pop_ranefs_df))
}

#now need to write out the blups into a csv file for the stepwise regression to run
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/3_lmm_blups/output/")

#writing out the file as a csv file
fwrite(all_pop_ranef_dfs, 
       na = "NA",
       file = "nam_blups_combined_univariate.csv", 
       sep = ",", 
       row.names = FALSE, 
       col.names = TRUE)

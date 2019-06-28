#this fille will combine the blups together
library(brms)
library(data.table)

rm(list = ls())

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/sun_and_shade_blups/input/")
pop_names <- c("blh_col", "bur_col", "cvi_col", "ita_col", "jea_col", "oy_col", "sha_col")

all_pop_blups <- c()

for (i in pop_names) {
  #reading in brms models for sun condition
  brms_sun_models <- lapply(list.files(include.dirs = TRUE, recursive = TRUE)[grepl(i, list.files(include.dirs = TRUE, recursive = TRUE)) & grepl("sun.rds", list.files(include.dirs = TRUE, recursive = TRUE)) & grepl("contr_trt_stud_10k", list.files(include.dirs = TRUE, recursive = TRUE))], function(x) readRDS(x))
  
  #splitting up file name to grab the phenotype
  file_names_sun <- list.files(include.dirs = TRUE, recursive = TRUE)[grepl(i, list.files(include.dirs = TRUE, recursive = TRUE)) & grepl("sun.rds", list.files(include.dirs = TRUE, recursive = TRUE)) & grepl("contr_trt_stud_10k", list.files(include.dirs = TRUE, recursive = TRUE))]
  file_names_sun_split_1 <- sapply(strsplit(file_names_sun, split = "/"), function(x) x[2])
  file_names_sun_split_2 <- sapply(strsplit(file_names_sun_split_1, split = "_brm_contr_"), function(x) x[1])
  pheno_names <- sapply(strsplit(file_names_sun_split_2, split = "_col_"), function(x) x[2])
  
  #grabbing the trait blups for the sun condition
  ranef_sun <- lapply(brms_sun_models, function(x) subset(ranef(x)$geno[, , 1], select = Estimate))
  
  #relabeling columns
  for (j in 1:length(ranef_sun)) {
    ranef_sun[[j]] <- as.data.table(ranef_sun[[j]], keep.rownames = TRUE)
    colnames(ranef_sun[[j]]) <- c("geno", pheno_names[j])
  }
  
  #merging sun blups together
  ranef_sun_merge <- Reduce(function(x, y) merge(x, y, by = "geno", all = TRUE), ranef_sun)
  ranef_sun_merge$treatment <- "Sun"
  
  #now working on shade blups
  brms_shade_models <- lapply(list.files(include.dirs = TRUE, recursive = TRUE)[grepl(i, list.files(include.dirs = TRUE, recursive = TRUE)) & grepl("shade.rds", list.files(include.dirs = TRUE, recursive = TRUE)) & grepl("contr_trt_stud_10k", list.files(include.dirs = TRUE, recursive = TRUE))], function(x) readRDS(x))
  file_names_shade <- list.files(include.dirs = TRUE, recursive = TRUE)[grepl(i, list.files(include.dirs = TRUE, recursive = TRUE)) & grepl("shade.rds", list.files(include.dirs = TRUE, recursive = TRUE)) & grepl("contr_trt_stud_10k", list.files(include.dirs = TRUE, recursive = TRUE))]
  file_names_shade_split_1 <- sapply(strsplit(file_names_shade, split = "/"), function(x) x[2])
  file_names_shade_split_2 <- sapply(strsplit(file_names_shade_split_1, split = "_brm_contr_"), function(x) x[1])
  pheno_names <- sapply(strsplit(file_names_shade_split_2, split = "_col_"), function(x) x[2])
  
  ranef_shade <- lapply(brms_shade_models, function(x) subset(ranef(x)$geno[, , 1], select = Estimate))
  
  for (j in 1:length(ranef_shade)) {
    ranef_shade[[j]] <- as.data.table(ranef_shade[[j]], keep.rownames = TRUE)
    colnames(ranef_shade[[j]]) <- c("geno", pheno_names[j])
  }
  
  ranef_shade_merge <- Reduce(function(x, y) merge(x, y, by = "geno", all = TRUE), ranef_shade)
  ranef_shade_merge$treatment <- "Shade"
  
  pop_sun_shade_ranef_comb <- rbindlist(list(ranef_sun_merge, ranef_shade_merge))
  pop_sun_shade_ranef_comb$pop <- i
  
  all_pop_blups <- rbindlist(list(all_pop_blups, pop_sun_shade_ranef_comb))
}

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/sun_and_shade_blups/output/")
fwrite(all_pop_blups, "pop_sun_shade_blups.csv", sep = ",", row.names = FALSE, col.names = TRUE)

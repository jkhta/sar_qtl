#this script will generate the lme4qtl models to get the standard errors of the effects
library(lme4)
library(lme4qtl)
library(data.table)
library(plyr)

rm(list = ls())

#reading in the kinship information
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/path_analysis_GridLMM/data/")
nam_kinship <- readRDS("nam_GridLMM_kinship.RDS")

#reading in the GridLMM stepwise models 
trait_qtl_fixef <- c()
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/nam_allelic_series/input/")
gridlmm_models <- list.files(pattern = "james_GridLMM_stepwise_model.RDS")

#function to merge the genotype and phenotype information into a single data frame
geno_pheno_merger <- function(pop_name, nam_data, sig_qtl) {
  pop_geno_pheno <- nam_data[[pop_name]]
  pop_pheno <- pop_geno_pheno$pop_pheno
  pop_geno <- as.data.table(pop_geno_pheno$pop_geno, keep.rownames = TRUE)
  colnames(pop_geno)[1] <- "geno"
  geno_col_subset <- which(colnames(pop_geno) %in% c("geno", gridlmm_qtl))
  pop_geno_qtl <- subset(pop_geno, select = geno_col_subset)
  pop_geno_pheno_merge <- merge(pop_pheno, pop_geno_qtl, by = "geno")
  return(pop_geno_pheno_merge)
}

#generating a table for the fixed effects and standard errors 
fixef_std_tabler <- function(pop_lme4qtl_model, nam_population, trait_of_interest, sig_qtl) {
  model_fixef <- fixef(pop_lme4qtl_model)
  model_fixef_std <- coef(summary(pop_lme4qtl_model))[, "Std. Error"]
  qtl_coef_ind <- tail(1:length(model_fixef), length(sig_qtl))
  qtl_fixef <- model_fixef[qtl_coef_ind]
  qtl_fixef_std <- model_fixef_std[qtl_coef_ind]
  qtl_coef_df <- data.frame(qtl = names(qtl_fixef),
                            fixef = qtl_fixef,
                            fixef_std = qtl_fixef_std,
                            trait = trait_of_interest,
                            pop = nam_population, 
                            stringsAsFactors = FALSE)
  return(qtl_coef_df)
}

#generating a kinship matrix with the QTL removed from the kinship estimation
proximal_kinship_generator <- function(geno_pheno_data, sig_qtl) {
  pop_geno <- geno_pheno_data$pop_geno
  pop_proximal <- geno_pheno_data$cM_proximal
  rownames(pop_proximal) <- colnames(pop_proximal) <- colnames(pop_geno)
  pop_proximal_subset <- pop_proximal[rownames(pop_proximal) %in% sig_qtl, ]
  
  #if trait only has 1 qtl then keep the vector the same; if the trait
  #has more than one than apply the sum function
  if (is.null(dim(pop_proximal_subset)) == TRUE) {
    pop_proximal_comb <- pop_proximal_subset
  } else {
    pop_proximal_comb <- unlist(apply(pop_proximal_subset, 2, function(x) sum(x)))
  }
  pop_proximal_tf <- ifelse(pop_proximal_comb > 0, FALSE, TRUE)
  pop_geno_subset <- pop_geno[, pop_proximal_tf]
  X  <- pop_geno_subset
  X2 <- sweep(X, 2, colMeans(X),"-")
  proximal_kinship <- tcrossprod(X2)/ncol(X2)
  return(proximal_kinship)
}

#reading in the models
for (i in 1:length(gridlmm_models)) {
  #grabbing the trait
  setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/nam_allelic_series/input/")
  traits <- sapply(strsplit(list.files(pattern = "james_GridLMM_stepwise_model.RDS"), split = "_james_"), function(x) x[1])
  trait <- traits[i]
  
  #grabbing the right model
  gridlmm_trait_model_name <- list.files(pattern = "james_GridLMM_stepwise_model.RDS")[i]
  gridlmm_trait_model <- readRDS(gridlmm_trait_model_name)
  
  #grabbing the qtl found
  gridlmm_qtl <- gridlmm_trait_model$found_qtl
  
  #grabbing the genotype and phenotype
  setwd("/Users/jkhta/Documents/GitHub/sar_qtl/6_push_button/")
  nam_geno_pheno <- readRDS("nam_all_traits_ind_pop_pheno_geno_proximal_james.RDS")
  
  #merging the genotype and phenotype information
  nam_geno_pheno_merge_list <- lapply(names(nam_geno_pheno), function(x) geno_pheno_merger(x, nam_geno_pheno, gridlmm_qtl))
  
  #generating formula for the responses; the shade responses are different because they have covariates  
  model_formula <- as.formula(paste(trait, " ~ ", paste(gridlmm_qtl, collapse = " + "), " + ", "(1|geno)", sep = ""))
  
  #generating the population counts 
  nam_pop_counts <- data.frame(pop = unlist(lapply(nam_geno_pheno_merge_list, function(x) unique(x$pop_pop))), 
                               count = unlist(lapply(nam_geno_pheno_merge_list, function(x) nrow(x))))
  
  #estimating a kinship matrix without the qtl 
  nam_kinship_proximal <- lapply(nam_geno_pheno, function(x) proximal_kinship_generator(x, gridlmm_qtl))
  
  #generating the linear mixed models to grab the fixed effects and their standard error estimates
  lme4qtl_model_list <- lapply(1:length(nam_geno_pheno_merge_list), function(x) relmatLmer(model_formula, nam_geno_pheno_merge_list[[x]], relmat = list(geno = nam_kinship_proximal[[x]])))
  names(lme4qtl_model_list) <- names(nam_geno_pheno)
  
  #generating a table of fixed effects for each trait and population
  fixef_table_list <- rbindlist(lapply(1:length(lme4qtl_model_list), function(x) fixef_std_tabler(pop_lme4qtl_model = lme4qtl_model_list[[x]],
                                                                                                  nam_population = names(lme4qtl_model_list)[[x]], 
                                                                                                  trait_of_interest = trait, 
                                                                                                  sig_qtl = gridlmm_qtl)))
  #mergin the fixed effects table with the counts to run the anova with sufficient statistics
  fixef_table_list_with_count <- merge(fixef_table_list, nam_pop_counts, by = "pop")
  
  trait_qtl_fixef <- rbindlist(list(trait_qtl_fixef, fixef_table_list_with_count))
} 

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/nam_allelic_series/data/")
trait_qtl_fixef$pop <- revalue(trait_qtl_fixef$pop, c("21RV_21RV" = "Blh-1",
                                                      "20RV_20RV" = "Bur-0",
                                                      "8RV_8RV" = "Cvi-0",
                                                      "29RV_29RV" = "Ita-0",
                                                      "28RV_28RV" = "Jea",
                                                      "27RV_27RV" = "Oy-0",
                                                      "13RV_13RV" = "Sha"))

fwrite(trait_qtl_fixef, "trait_lme4qtl_fixef_no_cov.csv", sep = ",", row.names = FALSE, col.names = TRUE)




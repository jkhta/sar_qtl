#this script will generate the accession heritability values for each trait
library(data.table)
library(brms)

rm(list = ls())

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/2_lmm_models/output/brms_models/accessions/")

table_formatter <- function(brms_object) {
  #grabbing the posterior samples
  pop_ps <- posterior_samples(brms_object, pars = required_parnames)
  
  #turning sd into variance
  pop_ps$sd_geno__Intercept <- pop_ps$sd_geno__Intercept^2
  pop_ps$sd_geno__treatment2 <- pop_ps$sd_geno__treatment2^2
  pop_ps$sigma <- pop_ps$sigma^2
  
  #renaming the variables
  colnames(pop_ps) <- c("shelf_fixef", "int_fixef","treatment_fixef", "geno_var", "gxe_var", "res")
  
  #calculating heritabilities
  pop_ps$geno_h2 <- with(pop_ps, geno_var / (geno_var + gxe_var + res))
  pop_ps$gxe_h2 <- with(pop_ps, gxe_var / (geno_var + gxe_var + res))
  
  return(pop_ps)
}

post_mean_and_ci <- function(post_vector) {
  #mean of posterior
  post_mean <- mean(post_vector)
  post_mean <- round(post_mean, digits = 4)
  
  #credible interval of posteriors
  post_ci <- quantile(post_vector, probs = c(0.025, 0.975))
  post_ci <- round(post_ci, digits = 4)
  post_ci_string <- paste(post_ci[1], post_ci[2], sep = " to ")
  
  return(data.frame(c(post_mean, post_ci_string)))
}

phenotype <- c("bd", "h3_h1", "r_dry", "i_dry")
required_parnames <- c("b_shelf", "b_Intercept","b_treatment2", "sd_geno__Intercept", "sd_geno__treatment2", "sigma")

all_trait_table <- c()

for (i in 1:length(phenotype)) {
  chosen_phenotype <- phenotype[i]

  brm_file_name <- list.files(pattern = chosen_phenotype, recursive = TRUE)
  
  acc_brms <- readRDS(brm_file_name)
  
  #want to have a table of credible intervals for block (shelf) effects, 
  #treatment effects, line (genotype) effects, and line by treatment effects (GxE)
  
  #need to pull out the samples for those variables
  acc_trait_table <- table_formatter(acc_brms)
  
  acc_trait_post_mean_table <- apply(acc_trait_table, 2, function(x) mean(x))
  names(acc_trait_post_mean_table) <- paste("acc", names(acc_trait_post_mean_table), sep = "_")
  acc_trait_post_mean_table <- as.data.frame(t(acc_trait_post_mean_table))
  acc_trait_post_mean_table$trait <- chosen_phenotype
  
  all_trait_table <- rbind(all_trait_table, acc_trait_post_mean_table)
}

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/h2_table/input/")
fwrite(all_trait_table, "acc_trait_effects_and_h2_no_ci.csv", sep = ",", row.names = FALSE, col.names = TRUE)

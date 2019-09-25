#this script will generate a table for the significances of the fixed, random, and heritabilities
#i will generate credible intervals for each effect
library(brms)
library(data.table)

rm(list = ls())

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/2_lmm_models/output/brms_models/nam/")

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
  pop_ps$geno_h2 <- with(pop_ps, (geno_var / (geno_var + gxe_var + res)) * 100)
  pop_ps$gxe_h2 <- with(pop_ps, (gxe_var / (geno_var + gxe_var + res)) * 100)
  pop_ps$e_h2 <- with(pop_ps, (res / (geno_var + gxe_var + res)) * 100)
  
  #calculating CV
  pop_ps$CV <- with(pop_ps, sqrt(gxe_var)/abs(treatment_fixef))
  
  return(pop_ps)
}

post_mean_and_ci <- function(post_vector) {
  #mean of posterior
  post_mean <- mean(post_vector)
  post_mean <- round(post_mean, digits = 2)
  
  #credible interval of posteriors
  post_ci <- quantile(post_vector, probs = c(0.025, 0.975))
  post_ci <- round(post_ci, digits = 2)
  post_ci_string <- paste(post_ci[1], post_ci[2], sep = " - ")
  
  return(data.frame(c(post_mean, post_ci_string)))
}

phenotype <- c("bd", "h3_h1", "r_dry", "i_dry")
required_parnames <- c("b_shelf", "b_Intercept", "b_treatment2", "sd_geno__Intercept", "sd_geno__treatment2", "sigma")

all_trait_table <- c()

for (i in 1:length(phenotype)) {
  chosen_phenotype <- phenotype[i]
  phenotype_brm <- paste(chosen_phenotype, "adapt_0.99_brm", sep = "_")
  
  brm_file_names <- list.files(pattern = phenotype_brm, recursive = TRUE)
  
  pop_brms <- lapply(brm_file_names, function(x) readRDS(x))
  
  #want to have a table of credible intervals for block (shelf) effects, 
  #treatment effects, line (genotype) effects, and line by treatment effects (GxE)
  
  #need to pull out the samples for those variables
  pop_tables <- lapply(pop_brms, function(x) table_formatter(x))
  
  pop_all_avg_values_list <- list()
  
  #now need to average the values
  for (j in colnames(pop_tables[[1]])) {
    pop_values <- do.call(cbind, lapply(pop_tables, function(x) subset(x, select = j)))
    pop_avg_values <- apply(pop_values, 1, function(x) mean(x))
    pop_avg_value_df <- data.frame(pop_avg_values)
    colnames(pop_avg_value_df) <- j
    pop_all_avg_values_list[[j]] <- pop_avg_value_df
  }
  
  #putting all effects into same df
  pop_all_avg_values_df <- do.call(cbind, pop_all_avg_values_list)
  
  #generating means and credible intervals
  trait_mean_and_ci <- do.call(cbind, apply(pop_all_avg_values_df, 2, function(x) post_mean_and_ci(x)))
  colnames(trait_mean_and_ci) <- colnames(pop_tables[[1]])
  rownames(trait_mean_and_ci) <- c(chosen_phenotype, paste(chosen_phenotype, "ci", sep = "_"))
  
  all_trait_table <- rbind(all_trait_table, trait_mean_and_ci)
}

#subsetting different data frames for different outputs
all_trait_table <- data.frame(data.table(all_trait_table, keep.rownames = TRUE), stringsAsFactors = FALSE)
colnames(all_trait_table)[1] <- "trait"
all_trait_table_ci_only <- subset(all_trait_table, grepl("_ci", trait))
all_trait_table_no_ci <- subset(all_trait_table, !grepl("_ci", trait))

all_trait_table_post_mean_and_ci <- c()

#creating data frame with both the posterior mean and credible intervals
for (i in 1:length(phenotype)) {
  phenotype_name <- phenotype[i]
  phenotype_post_mean <- unlist(subset(all_trait_table_no_ci, trait == phenotype_name)[, 2:ncol(all_trait_table_no_ci)])
  phenotype_post_ci <- unlist(subset(all_trait_table_ci_only, grepl(phenotype_name, trait))[, 2:ncol(all_trait_table_ci_only)])
  phenotype_post_mean_and_ci <- paste(phenotype_post_mean, " (", phenotype_post_ci, ")", sep = "")
  all_trait_table_post_mean_and_ci <- rbindlist(list(all_trait_table_post_mean_and_ci, data.frame(t(phenotype_post_mean_and_ci), stringsAsFactors = FALSE)))
}

all_trait_table_post_mean_and_ci <- cbind(data.frame(trait = phenotype, stringsAsFactors = FALSE), all_trait_table_post_mean_and_ci)
colnames(all_trait_table_post_mean_and_ci) <- colnames(all_trait_table_no_ci)
  
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/h2_table/input/")
fwrite(all_trait_table, file = "trait_effects_and_h2.csv", sep = ",", row.names = FALSE, col.names = TRUE)
fwrite(all_trait_table_ci_only, file = "trait_effects_and_h2_ci.csv", sep = ",", row.names = FALSE, col.names = TRUE)
fwrite(all_trait_table_no_ci, file = "trait_effects_and_h2_no_ci.csv", sep = ",", row.names = FALSE, col.names = TRUE)
fwrite(all_trait_table_post_mean_and_ci, file = "trait_effects_and_h2_post_mean_and_ci.csv", sep = ",", row.names = FALSE, col.names = TRUE)


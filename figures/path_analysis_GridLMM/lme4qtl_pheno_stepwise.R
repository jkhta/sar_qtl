#this script will generate a stepwise algorithm for the path analysis QTL scan
library(lme4)
library(lmerTest)
library(data.table)
library(lme4qtl)

rm(list = ls())

setwd("/Users/James/Documents/GitHub/sar_qtl/figures/path_analysis_GridLMM/")

kinship_mat <- readRDS("nam_GridLMM_kinship.RDS")

nam_pheno <- fread("nam_blups_combined_univariate.csv",
                   sep = ",",
                   header = TRUE)
nam_pheno$family <- paste(sapply(strsplit(nam_pheno$geno, "RV"), function(x) x[1]), "RV", sep = "")
nam_pheno$geno <- paste(nam_pheno$family, nam_pheno$geno, sep = "_")

nam_pheno_list <- lapply(unique(nam_pheno$family), function(x) as.data.frame(subset(nam_pheno, family == x)))

for (i in 1:length(unique(nam_pheno$family))) {
  nam_pheno_list_copy <- nam_pheno_list[[i]][match(rownames(kinship_mat[[i]]), nam_pheno_list[[i]]$geno), ]
  nam_pheno_list[[i]] <- nam_pheno_list_copy
}

nam_pheno_order <- c("bd", "r_dry", "h3_h1", "i_dry")


lme4qtl_test <- function(formula_red, formula_full, pheno_data, kinship) {
  kinship <- kinship
  reduced_model <- relmatLmer(formula_red,
                              data = pheno_data,
                              relmat = list(geno = kinship),
                              REML = F,
                              verbose = FALSE)
  
  full_model <- relmatLmer(formula_full, 
                           data = pheno_data, 
                           relmat = list(geno = kinship), 
                           REML = F, 
                           verbose = FALSE)
  
  return(anova(reduced_model, full_model))
}


test <- lapply(1:length(kinship_mat), function(x) lme4qtl_test(formula_red = "i_dry_gxe ~ (1|geno)", 
                                                               formula_full = "i_dry_gxe ~ bd_gxe + r_dry_gxe + h3_h1_gxe + (1|geno)", 
                                                               pheno_data = nam_pheno_list[[x]], 
                                                               kinship = kinship_mat[[x]]))

significance_tester <- function(pop_anova) {
  anova_table <- as.data.frame(pop_anova)
  lik_diff <- 2 * (pop_anova$logLik[2] - pop_anova$logLik[1])
}
               
pchisq(sum(unlist(lapply(test, function(x) significance_tester(x)))), df = length(kinship_mat), lower.tail = FALSE)
               
model_formulas <- list(r_dry = "r_dry ~ bd",
                       h3_h1 = "h3_h1 ~ bd + r_dry",
                       i_dry = "i_dry ~ bd + r_dry + h3_h1")
               
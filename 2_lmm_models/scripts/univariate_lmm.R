library(brms)
library(data.table)

rm(list = ls())

#run for each population; each cluster is a different run or population
run = as.numeric(commandArgs(t=T)[1])

#reading in data
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/1_phenotype_standardization/output/")
nam_cam_data_transformed <- read.csv("nam_cam_data_combined_tformed_std.csv",
                                     header = TRUE,
                                     stringsAsFactors = FALSE)

#subsetting the data based on population
nam_pops <- c("blh_col", "bur_col", "cvi_col", "ita_col", "jea_col", "oy_col", "sha_col")
nam_pop <- nam_pops[run]
nam_cam_data_subset <- subset(nam_cam_data_transformed, cross == nam_pop)
nam_cam_data_subset$treatment <- factor(nam_cam_data_subset$treatment, levels = c("Sun", "Shade"))

#choosing the phenotypes to generate a bayesian mixed model for
nam_cam_data_subset_phenotypes <- c("bd", "h3_h1", "i_dry", "r_dry")
contrasts(nam_cam_data_subset$treatment) <- contr.treatment

#setting a different working directory for each population
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/2_lmm_models/output/")
try(dir.create(nam_pop))
setwd(nam_pop)

#input values for number of iterations and the adapt value
iterations <- 10000
adapt_value <- 0.99

#for each phenotype in the population, generate a brms model, or bayesian mixed model
for (i in nam_cam_data_subset_phenotypes) {
  single_brm_model_name <- paste(nam_pop, i, "adapt", adapt_value, "brm_contr_trt_stud_10k.rds", sep = "_")
  single_brm_formula <- as.formula(paste(i, "~ shelf + treatment + (1 + treatment|geno)", sep = " "))
  single_brm_model <- brm(single_brm_formula,
                          data = nam_cam_data_subset,
                          family = student(),
                          iter = iterations,
                          cores = 4,
                          seed = 13,
                          control = list(adapt_delta = adapt_value))
  saveRDS(single_brm_model, single_brm_model_name)
}

#this script fits a bayesian linear mixed model to each popuation, for each trait, for sun and shade conditions separately
#this is to generate the BLUPs used later for the path analysis
library(brms)

rm(list = ls())

#run for each population; each cluster is a different run or population
run = as.numeric(commandArgs(t=T)[1])

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/1_phenotype_standardization/output/")
nam_cam_data_transformed <- read.csv("nam_cam_data_combined_tformed_std.csv",
                                     header = TRUE,
                                     stringsAsFactors = FALSE)

#subsetting by population
nam_pops <- c("blh_col", "bur_col", "cvi_col", "ita_col", "jea_col", "oy_col", "sha_col")
nam_pop <- nam_pops[run]
nam_cam_data_subset <- subset(nam_cam_data_transformed, cross == nam_pop)
nam_cam_data_subset$treatment <- factor(nam_cam_data_subset$treatment, levels = c("Sun", "Shade"))

#subsetting based on treatment
nam_cam_data_subset_sun <- subset(nam_cam_data_subset, treatment == "Sun")
nam_cam_data_subset_shade <- subset(nam_cam_data_subset, treatment == "Shade")

#choosing the phenotypes to generate a bayesian mixed model for
nam_cam_data_subset_phenotypes <- c("bd", "h3_h1", "i_dry", "r_dry")
contrasts(nam_cam_data_subset$treatment) <- contr.treatment

#setting a different working directory for each population
setwd("/home/jkta/projects/col_sha/brms/output/brms_models")
dir_name <- paste(nam_pop, "sun_shade", sep = "_")
try(dir.create(dir_name))
setwd(nam_pop)

iterations <- 10000
adapt_value <- 0.99

#fitting bayesian mixed models using the brm package to the sun and shade data separately
for (i in nam_cam_data_subset_phenotypes) {
  sun_brm_model_name <- paste(nam_pop, i, "brm_contr_trt_stud_10k_sun.rds", sep = "_")
  shade_brm_model_name <- paste(nam_pop, i, "brm_contr_trt_stud_10k_shade.rds", sep = "_")
  
  single_brm_formula <- as.formula(paste(i, "~ shelf + (1|geno)", sep = " "))
  
  #sun model
  sun_brm_model <- brm(single_brm_formula,
                       data = nam_cam_data_subset_sun,
                       family = student(),
                       iter = iterations,
                       cores = 4,
                       seed = 13)

  saveRDS(sun_brm_model, sun_brm_model_name)
  
  #shade model
  shade_brm_model <- brm(single_brm_formula,
                         data = nam_cam_data_subset_shade,
                         family = student(),
                         iter = iterations,
                         cores = 4,
                         seed = 13)
  saveRDS(shade_brm_model, shade_brm_model_name)
}



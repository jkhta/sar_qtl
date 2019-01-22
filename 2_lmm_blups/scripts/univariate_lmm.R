library(brms)

rm(list = ls())
run = as.numeric(commandArgs(t=T)[1])

setwd("/home/jkta/projects/col_sha/brms/input")

data_type <- "not_bd_std"

if (data_type == "bd_std") {
  nam_cam_data_transformed <- read.csv("nam_cam_data_combined_tformed_std.csv",
                                       header = TRUE,
                                       stringsAsFactors = FALSE)
} else if (data_type == "not_bd_std") {
  nam_cam_data_transformed <- read.csv("nam_cam_data_combined_tformed_std_FINAL.csv",
                                       header = TRUE,
                                       stringsAsFactors = FALSE)
}

nam_pops <- c("blh", "bur", "cvi", "ita", "jea", "oy", "sha")
nam_pop <- nam_pops[run]
nam_cam_data_subset <- subset(nam_cam_data_transformed, subset = grepl(nam_pop, cross))
nam_cam_data_subset$treatment <- factor(nam_cam_data_subset$treatment, levels = c("Sun", "Shade"))

nam_cam_data_subset_sun <- subset(nam_cam_data_subset, treatment == "Sun")
nam_cam_data_subset_shade <- subset(nam_cam_data_subset, treatment == "Shade")

nam_cam_data_subset_phenotypes <- c("bd", "h2_h1", "h3_h2", "h3_h1", "i_dry", "r_dry")
contrasts(nam_cam_data_subset$treatment) <- contr.treatment

#setting a different working directory for each population
setwd("/home/jkta/projects/col_sha/brms/output/brms_models")
try(dir.create(nam_pop))
setwd(nam_pop)

iterations <- 10000
adapt_value <- 0.99

#THIS ONE RUN FOR SUN AND SHADE SEPARATELY
# for (i in nam_cam_data_subset_phenotypes) {
#   sun_brm_model_name <- paste(nam_pop, i, "brm_contr_trt_stud_10k_sun.rds", sep = "_")
#   shade_brm_model_name <- paste(nam_pop, i, "brm_contr_trt_stud_10k_shade.rds", sep = "_")
# 
#   single_brm_formula <- as.formula(paste(i, "~ shelf + (1|geno)", sep = " "))
# 
#   sun_brm_model <- brm(single_brm_formula,
#                        data = nam_cam_data_subset_sun,
#                        family = student(),
#                        iter = iterations,
#                        cores = 4,
#                        seed = 13)
# 
#   saveRDS(sun_brm_model, sun_brm_model_name)
# 
#   shade_brm_model <- brm(single_brm_formula,
#                          data = nam_cam_data_subset_shade,
#                          family = student(),
#                          iter = iterations,
#                          cores = 4,
#                          seed = 13)
#   saveRDS(shade_brm_model, shade_brm_model_name)
# }
# 
# 

#THIS ONE RUN FOR UNIVARIATE
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

# 
# multivariate_brm_sun_model_name <- paste(nam_pop, "multivar_brm_contr_trt_stud_10k_sun.rds", sep = "_")
# 
# multivariate_brm_sun_model <- brm(cbind(bd, h2_h1, h3_h2, h3_h1, i_dry, r_dry) ~ shelf + (1|1|geno), 
#                                   data = nam_cam_data_subset_sun, 
#                                   family = student(), 
#                                   iter = 10000, 
#                                   cores = 4, 
#                                   seed = 13)
# 
# saveRDS(multivariate_brm_sun_model, multivariate_brm_sun_model_name)
# 
# multivariate_brm_shade_model_name <- paste(nam_pop, "multivar_brm_contr_trt_stud_10k_shade.rds", sep = "_")
# 
# multivariate_brm_shade_model <- brm(cbind(bd, h2_h1, h3_h2, h3_h1, i_dry, r_dry) ~ shelf + (1|1|geno), 
#                                     data = nam_cam_data_subset_shade, 
#                                     family = student(), 
#                                     iter = 10000, 
#                                     cores = 4, 
#                                     seed = 13)
# 
# saveRDS(multivariate_brm_shade_model, multivariate_brm_shade_model_name)

#THIS ONE RUN FOR MULTIVARIATE
# iteration_name <- paste(iterations/1000, "k", sep = "")
# univariate_brm_model_name <- paste(nam_pop, data_type, iteration_name, "var_brm_contr_trt_stud.rds", sep = "_")
# multivariate_brm_model_name <- paste(nam_pop, data_type, iteration_name, "multivar_brm_contr_trt_stud.rds", sep = "_")
# 
# multivariate_brm_model <- brm(cbind(bd, h3_h1, i_dry, r_dry) ~ shelf + treatment + (1 + treatment|p|geno), 
#                               data = nam_cam_data_subset, 
#                               family = student(), 
#                               iter = iterations,
#                               chains = 4,
#                               cores = 4, 
#                               seed = 13)
# 
# saveRDS(multivariate_brm_model, multivariate_brm_model_name)

#now running models for individual traits

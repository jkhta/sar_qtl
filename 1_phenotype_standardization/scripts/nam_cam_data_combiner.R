#now the ultimate nam cam testing
library(data.table)
library(lme4)
library(car)
library(lmerTest)
library(brms)
library(sjmisc)
library(optimx)

rm(list = ls())
setwd("/Users/jkhta/Desktop/nam_cam_fixing/10.5 - nam_cam_other_rep_data_cleaning/output/")

nam_cam_data_files <- list.files(pattern = "*phenotypes_combined.csv")
nam_cam_data_list <- lapply(nam_cam_data_files, function(x) read.csv(x, header = TRUE, stringsAsFactors = FALSE))

for (i in 1:length(nam_cam_data_list)) {
  nam_cam_data_rep <- nam_cam_data_list[[i]]
  nam_cam_data_rep$assigned_no <- NULL
  nam_cam_data_list[[i]] <- nam_cam_data_rep
}

nam_cam_data <- as.data.frame(rbindlist(nam_cam_data_list))
nam_cam_data$cross <- "accession"
nam_cam_data[grepl("^8RV", nam_cam_data$geno), ]$cross <- "cvi_col"
nam_cam_data[grepl("^13RV", nam_cam_data$geno), ]$cross <- "sha_col"
nam_cam_data[grepl("^20RV", nam_cam_data$geno), ]$cross <- "bur_col"
nam_cam_data[grepl("^21RV", nam_cam_data$geno), ]$cross <- "blh_col"
nam_cam_data[grepl("^27RV", nam_cam_data$geno), ]$cross <- "oy_col"
nam_cam_data[grepl("^28RV", nam_cam_data$geno), ]$cross <- "jea_col"
nam_cam_data[grepl("^29RV", nam_cam_data$geno), ]$cross <- "ita_col"

nam_cam_info <- c("plant_id", "tray", "pot_pos", "geno", "geno2", "shelf", "treatment", "envelope", "rep", "cross")
nam_cam_info_data <- nam_cam_data[, (colnames(nam_cam_data) %in% nam_cam_info)]
nam_cam_data_phenotypes <- nam_cam_data[, !(colnames(nam_cam_data) %in% nam_cam_info)]
apply(nam_cam_data_phenotypes, 2, function(x) any(is.infinite(x)))
apply(nam_cam_data_phenotypes, 2, function(x) any(x < 0, na.rm = TRUE))
apply(nam_cam_data_phenotypes, 2, function(x) any(x == 0, na.rm = TRUE))

nam_cam_data_phenotypes_clean <- do.call(data.frame, lapply(nam_cam_data_phenotypes, function(x) replace(x, is.infinite(x),NA)))
nam_cam_data_phenotypes_clean <- do.call(data.frame, lapply(nam_cam_data_phenotypes_clean, function(x) replace(x, which(x < 0),NA)))
nam_cam_data_phenotypes_clean <- do.call(data.frame, lapply(nam_cam_data_phenotypes_clean, function(x) replace(x, which(x == 0),NA)))

apply(nam_cam_data_phenotypes_clean, 2, function(x) any(is.infinite(x)))
apply(nam_cam_data_phenotypes_clean, 2, function(x) any(x < 0, na.rm = TRUE))
apply(nam_cam_data_phenotypes_clean, 2, function(x) any(x == 0, na.rm = TRUE))

nam_cam_clean <- nam_cam_data
nam_cam_clean[, !(colnames(nam_cam_clean) %in% nam_cam_info)] <- nam_cam_data_phenotypes_clean

write.csv(nam_cam_clean, "nam_cam_data_combined.csv", row.names = FALSE)

nam_cam_clean$treatment <- factor(nam_cam_clean$treatment, levels = c("Sun", "Shade"))
contrasts(nam_cam_clean$treatment) <- contr.sum

nam_cam_clean_phenotypes <- nam_cam_clean[, !(colnames(nam_cam_clean) %in% nam_cam_info)]

#do NOT standardize each other; that's where the path analysis comes in
#nam_cam_clean_phenotypes$r_dry <- with(nam_cam_clean_phenotypes, r_dry/days)
#nam_cam_clean_phenotypes$i_dry <- with(nam_cam_clean_phenotypes, i_dry/i_growth_days)

nam_cam_good_phenotypes_names <- c("days",
                                   "h2_h1_height_diff_avg",
                                   "h3_h2_height_diff_avg",
                                   "h3_h1_height_diff_avg",
                                   "i_dry",
                                   "r_dry")

nam_cam_good_phenotypes <- nam_cam_clean_phenotypes[, nam_cam_good_phenotypes_names]
names(nam_cam_good_phenotypes) <- c("bd", "h2_h1", "h3_h2", "h3_h1", "i_dry", "r_dry")

for (i in 1:ncol(nam_cam_good_phenotypes)) {
  phenotype_data <- nam_cam_good_phenotypes[, i]
  phenotype_model <- lmer(phenotype_data ~ shelf + treatment + (1 + treatment|geno), 
                          data = nam_cam_clean, 
                          REML = TRUE)
  phenotype_transform <- powerTransform(phenotype_model, family = "bcPower")
  phenotype_transformed_data <- bcPower(phenotype_data, phenotype_transform$lambda)
  phenotype_transformed_standardized_data <- std(phenotype_transformed_data)
  nam_cam_good_phenotypes[, i] <- phenotype_transformed_standardized_data
}

nam_cam_clean_tformed_std <- cbind(nam_cam_info_data, nam_cam_good_phenotypes)

write.csv(nam_cam_clean_tformed_std, "nam_cam_data_combined_tformed_std_FINAL.csv", row.names = FALSE)


#------------THIS SECTION WILL USE BRMS TO TRANSFORM THE PHENOTYPES-------------

#the brm model takes too long; later on i might want to try to use this model to powertransform
#instead of the lmer model
#generating a data frame with the experimental information and the renamed phenotypes
#need to generate another data frame with the renamed phenotypes and other information
#to fit the brms models
nam_cam_info <- subset(nam_cam_clean, select = c(shelf, treatment, geno))
nam_cam_clean_comp <- cbind(nam_cam_info, nam_cam_good_phenotypes)

phenotype_data <- nam_cam_good_phenotypes[, i]
phenotype_name <- colnames(nam_cam_good_phenotypes)[i]
phenotype_formula <- as.formula(paste(phenotype_name, "~ shelf + treatment + (1 + treatment|geno)", sep = " "))
#the brms model takes about ~ 10 minutes to run

#brms output
phenotype_model_1 <- brm(phenotype_formula, 
                       data = nam_cam_clean_comp, 
                       family = student(), 
                       chains = 4,
                       cores = 4, 
                       seed = 13)
phenotype_transform_1 <- powerTransform(phenotype_model_1, family = "bcPower")
phenotype_transformed_data_1 <- bcPower(phenotype_data, phenotype_transform_1$lambda)
phenotype_transformed_standardized_data_1 <- std(phenotype_transformed_data_1)

#comparing to lmer output
phenotype_data <- nam_cam_good_phenotypes[, i]
phenotype_model_2 <- lmer(phenotype_data ~ shelf + treatment + (1 + treatment|geno), 
                        data = nam_cam_clean, 
                        REML = TRUE)
phenotype_transform_2 <- powerTransform(phenotype_model_2, family = "bcPower")
phenotype_transformed_data_2 <- bcPower(phenotype_data, phenotype_transform_2$lambda)
phenotype_transformed_standardized_data_2 <- std(phenotype_transformed_data)

plot(phenotype_transformed_standardized_data_1, phenotype_transformed_standardized_data_2)



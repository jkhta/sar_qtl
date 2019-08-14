#combines data from all reps, removes erroneous data, and then transforms and standardizes the data
library(data.table)
library(lme4)
library(car)
library(lmerTest)
library(brms)
library(sjmisc)
library(optimx)

rm(list = ls())
setwd("/Users/James/Documents/GitHub/sar_qtl/1_phenotype_standardization/input/")

#reading in the phenotype data for all experiments as a list
nam_cam_data_files <- list.files(pattern = "*phenotypes_combined.csv")
nam_cam_data_list <- lapply(nam_cam_data_files, function(x) read.csv(x, header = TRUE, stringsAsFactors = FALSE))

#getting rid of the assigned numbers column because it serves no purpose and interferes
#with the rbindlist later
for (i in 1:length(nam_cam_data_list)) {
  nam_cam_data_rep <- nam_cam_data_list[[i]]
  nam_cam_data_rep$assigned_no <- NULL
  nam_cam_data_list[[i]] <- nam_cam_data_rep
}

#adding a cross name column
nam_cam_data <- as.data.frame(rbindlist(nam_cam_data_list))
nam_cam_data$cross <- "accession"
nam_cam_data[grepl("^8RV", nam_cam_data$geno), ]$cross <- "cvi_col"
nam_cam_data[grepl("^13RV", nam_cam_data$geno), ]$cross <- "sha_col"
nam_cam_data[grepl("^20RV", nam_cam_data$geno), ]$cross <- "bur_col"
nam_cam_data[grepl("^21RV", nam_cam_data$geno), ]$cross <- "blh_col"
nam_cam_data[grepl("^27RV", nam_cam_data$geno), ]$cross <- "oy_col"
nam_cam_data[grepl("^28RV", nam_cam_data$geno), ]$cross <- "jea_col"
nam_cam_data[grepl("^29RV", nam_cam_data$geno), ]$cross <- "ita_col"

#subsetting the data by specific columns
nam_cam_info <- c("plant_id", "tray", "pot_pos", "geno", "geno2", "shelf", "treatment", "envelope", "rep", "cross")
nam_cam_info_data <- nam_cam_data[, (colnames(nam_cam_data) %in% nam_cam_info)]

#grabbing only the phenotype data
nam_cam_data_phenotypes <- nam_cam_data[, !(colnames(nam_cam_data) %in% nam_cam_info)]

#checking for bad data points: infinite, negative, or 0; FALSE is good
apply(nam_cam_data_phenotypes, 2, function(x) any(is.infinite(x)))
apply(nam_cam_data_phenotypes, 2, function(x) any(x < 0, na.rm = TRUE))
apply(nam_cam_data_phenotypes, 2, function(x) any(x == 0, na.rm = TRUE))

#replacing bad data points with an NA
nam_cam_data_phenotypes_clean <- do.call(data.frame, lapply(nam_cam_data_phenotypes, function(x) replace(x, is.infinite(x),NA)))
nam_cam_data_phenotypes_clean <- do.call(data.frame, lapply(nam_cam_data_phenotypes_clean, function(x) replace(x, which(x < 0),NA)))
nam_cam_data_phenotypes_clean <- do.call(data.frame, lapply(nam_cam_data_phenotypes_clean, function(x) replace(x, which(x == 0),NA)))

#checking again to see if bad data points are there; FALSE is good
apply(nam_cam_data_phenotypes_clean, 2, function(x) any(is.infinite(x)))
apply(nam_cam_data_phenotypes_clean, 2, function(x) any(x < 0, na.rm = TRUE))
apply(nam_cam_data_phenotypes_clean, 2, function(x) any(x == 0, na.rm = TRUE))

#combining the information with the cleaned phenotypes
nam_cam_clean <- cbind(nam_cam_info_data, nam_cam_data_phenotypes_clean)

#writing out the cleaned phenotype data
#setwd("/Users/jkhta/Documents/GitHub/sar_qtl/1_phenotype_standardization/output/")
#write.csv(nam_cam_clean, "nam_cam_data_combined.csv", row.names = FALSE)

#reordering the factors in the data set
nam_cam_clean$treatment <- factor(nam_cam_clean$treatment, levels = c("Sun", "Shade"))

#changing the contrasts so it's sun + shade/2
contrasts(nam_cam_clean$treatment) <- contr.sum

#subsetting the phenotype only
nam_cam_clean_phenotypes <- nam_cam_clean[, !(colnames(nam_cam_clean) %in% nam_cam_info)]

#subsetting by these phenotypes only
nam_cam_good_phenotypes_names <- c("days",
                                   "h3_h1_height_diff_avg",
                                   "i_dry",
                                   "r_dry")
nam_cam_good_phenotypes <- nam_cam_clean_phenotypes[, nam_cam_good_phenotypes_names]

#changing the phenotype names so they're easier to use
names(nam_cam_good_phenotypes) <- c("bd", "h3_h1", "i_dry", "r_dry")

#using a power transform on the data
for (i in 1:ncol(nam_cam_good_phenotypes)) {
  phenotype_data <- nam_cam_good_phenotypes[, i]
  phenotype_model <- lm(phenotype_data ~ shelf + treatment * geno, 
                        data = nam_cam_clean, 
                        REML = TRUE)
  phenotype_transform <- powerTransform(phenotype_model, family = "bcPower")
  phenotype_transformed_data <- bcPower(phenotype_data, phenotype_transform$lambda)
  phenotype_transformed_standardized_data <- std(phenotype_transformed_data)
  nam_cam_good_phenotypes[, i] <- phenotype_transformed_standardized_data
}

#combining the data information with the transformed and standardized phenotypes
nam_cam_clean_tformed_std <- cbind(nam_cam_info_data, nam_cam_good_phenotypes)

#writing out the transformed and standardized data set
write.csv(nam_cam_clean_tformed_std, "nam_cam_data_combined_tformed_std.csv", row.names = FALSE)



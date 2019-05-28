#this script will calculate the GxE in natural accessions
library(data.table)
library(ggplot2)
library(lme4)
library(sjPlot)
library(sjmisc)

rm(list = ls())

setwd("/Users/jkhta/Desktop/nam_cam_fixing/10.5 - nam_cam_other_rep_data_cleaning/output/")

data <- fread("nam_cam_data_combined_tformed_std_FINAL.csv",
              sep = ",",
              header = TRUE,
              stringsAsFactors = FALSE)
data$treatment <- factor(data$treatment, levels = c("Sun", "Shade"))
data$treatment <- as.numeric(data$treatment)

head(data)

data_acc <- subset(data, cross == "accession")

acc_lmer <- lmer(bd ~ treatment + (1 + treatment|geno), 
                 data = data_acc)

summary(acc_lmer)

0.0002897/(0.8212640 + 0.0002897 + 0.4265306)

plot_model(acc_lmer, 
           type = "pred", 
           terms = c("treatment", "geno"),
           show.legend = FALSE) + 
  theme(axis.text=element_text(size = 15),
        axis.title=element_text(size = 30)) +
  ylab("Bolting time")

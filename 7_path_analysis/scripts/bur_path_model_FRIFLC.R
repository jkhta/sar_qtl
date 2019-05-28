#path analysis on bur group
library(data.table)
library(lavaan)
library(semPlot)

rm(list = ls())

setwd("/Users/jkhta/Desktop/nam_cam_fixing/34 - pop_path_analysis/output/")

bur_data <- fread("bur_col_pheno_qtl_combo_merge.csv",
                  sep = ",",
                  header = TRUE,
                  stringsAsFactors = FALSE)
setnames(bur_data, c("r_dry", "h3_h1", "i_dry"), c("rdry", "h3h1", "idry"))
bur_data$treatment <- factor(bur_data$treatment, levels = c("Sun", "Shade"))
head(bur_data)

bur_data[, c("A4", "B4", "A5", "B5", "AA", "AB", "BB")] <- NULL

m_best <- "
rdry ~ c(A, a) * bd + c(B, b) * BA
h3h1 ~ c(C, c) * bd + c(D, d) * rdry
idry ~ c(E, e) * bd + c(F, f) * rdry + c(G, g) * h3h1 + c(H, h) * BA

#FRIFLC

#bd
FRIFLC_bd_dir_sun := 0
FRIFLC_bd_dir_shade := 0
FRIFLC_bd_dir_diff := 0

FRIFLC_bd_ind_sun := 0
FRIFLC_bd_ind_shade := 0
FRIFLC_bd_ind_diff := 0

#rdry
FRIFLC_rdry_dir_sun := b
FRIFLC_rdry_dir_shade := B
FRIFLC_rdry_dir_diff := B - b

FRIFLC_rdry_ind_sun := 0
FRIFLC_rdry_ind_shade := 0
FRIFLC_rdry_ind_diff := 0

#h3h1
FRIFLC_h3h1_dir_sun := 0
FRIFLC_h3h1_dir_shade := 0
FRIFLC_h3h1_dir_diff := 0

FRIFLC_h3h1_ind_sun := b * d
FRIFLC_h3h1_ind_shade := B * D
FRIFLC_h3h1_ind_diff := (B * D) - (b * d)

#idry
FRIFLC_idry_dir_sun := h
FRIFLC_idry_dir_shade := H
FRIFLC_idry_dir_diff := H - h

FRIFLC_idry_ind_sun := (b * f) + (b * d * g)
FRIFLC_idry_ind_shade := (B * F) + (B * D * G)
FRIFLC_idry_ind_diff := ((B * F) + (B * D * G)) - ((b * f) + (b * d * g))
"

fit_best <- sem(m_best, data = bur_data, group = "treatment")
summary(fit_best, fit.measures = TRUE)

semPaths(fit_best, 
         what = "path", 
         whatLabels = "name", 
         layout = "tree", 
         edge.label.cex = 1.5,
         residuals = FALSE,
         exoCov = FALSE,
         exoVar = FALSE,
         intercepts = FALSE, 
         #edge.color = FRI_sha_FLC_col_sun_coef, 
         fade = FALSE, 
         combineGroups = FALSE,
         nCharNodes = 0)

path_parameters <- parameterestimates(fit_best, standardized = TRUE)
path_parameters_eff <- subset(path_parameters, grepl("sun|shade|diff", lhs))
path_parameters_eff$allele_comb <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[1])
path_parameters_eff$trait <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[2])
path_parameters_eff$effect_type <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[3])
path_parameters_eff$env <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[4])
nrow(path_parameters_eff)
path_parameters_eff$pop <- "bur"

setwd("/Users/jkhta/Desktop/nam_cam_fixing/34 - pop_path_analysis/output/")
fwrite(path_parameters_eff, "bur_path_eff_FRIFLC.csv", sep = ",", row.names = FALSE)

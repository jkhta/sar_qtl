#path analysis on blh group
library(data.table)
library(lavaan)
library(semPlot)

rm(list = ls())

setwd("/Users/jkhta/Desktop/nam_cam_fixing/34 - pop_path_analysis/output/")

blh_data <- fread("blh_col_pheno_qtl_combo_merge.csv",
                  sep = ",",
                  header = TRUE,
                  stringsAsFactors = FALSE)
setnames(blh_data, c("r_dry", "h3_h1", "i_dry"), c("rdry", "h3h1", "idry"))
blh_data$treatment <- factor(blh_data$treatment, levels = c("Sun", "Shade"))
head(blh_data)

blh_data[, c("A4", "A5", "B4", "B5", "AA", "AB")] <- NULL
blh_data$FRIFLC <- ifelse(blh_data$BA > 0.5 | blh_data$BB > 0.5, 1, 0)

m_best <- "
bd ~ c(a, A) * FRIFLC
rdry ~ c(b, B) * bd + c(c, C) * FRIFLC
h3h1 ~ c(d, D) * bd
idry ~ c(e, E) * bd + c(f, F) * rdry + c(g, G) * h3h1 + c(h, H) * FRIFLC

#FRIFLC

#bd
FRIFLC_bd_dir_sun := a
FRIFLC_bd_dir_shade := A
FRIFLC_bd_dir_diff := A - a

FRIFLC_bd_ind_sun := 0
FRIFLC_bd_ind_shade := 0
FRIFLC_bd_ind_diff := 0

#rdry
FRIFLC_rdry_dir_sun := c
FRIFLC_rdry_dir_shade := C
FRIFLC_rdry_dir_diff := C - c

FRIFLC_rdry_ind_sun := a * b
FRIFLC_rdry_ind_shade := A * B
FRIFLC_rdry_ind_diff := (A * B) - (a * b)

#h3h1
FRIFLC_h3h1_dir_sun := 0
FRIFLC_h3h1_dir_shade := 0
FRIFLC_h3h1_dir_diff := 0

FRIFLC_h3h1_ind_sun := a * d
FRIFLC_h3h1_ind_shade := A * D
FRIFLC_h3h1_ind_diff := (A * D) - (a * d)

#idry
FRIFLC_idry_dir_sun := h
FRIFLC_idry_dir_shade := H
FRIFLC_idry_dir_diff := H - h

FRIFLC_idry_ind_sun := (a * e) + (a * b * f) + (c * f)
FRIFLC_idry_ind_shade := (A * E) + (A * B * F) + (C * F)
FRIFLC_idry_ind_diff := ((A * E) + (A * B * F) + (C * F)) - ((a * e) + (a * b * f) + (c * f))
"

fit_best <- sem(m_best, data = blh_data, group = "treatment")
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
path_parameters_eff$pop <- "blh"

setwd("/Users/jkhta/Desktop/nam_cam_fixing/34 - pop_path_analysis/output/")
fwrite(path_parameters_eff, "blh_path_eff_FRIFLC.csv", sep = ",", row.names = FALSE)

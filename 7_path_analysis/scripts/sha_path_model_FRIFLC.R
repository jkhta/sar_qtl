#path analysis on sha group
library(data.table)
library(lavaan)
library(semPlot)

rm(list = ls())

setwd("/Users/jkhta/Desktop/nam_cam_fixing/34 - pop_path_analysis/output/")

sha_data <-fread("sha_col_pheno_qtl_combo_merge.csv",
                 sep = ",",
                 header = TRUE,
                 stringsAsFactors = FALSE)
setnames(sha_data, c("r_dry", "h3_h1", "i_dry"), c("rdry", "h3h1", "idry"))
sha_data$treatment <- factor(sha_data$treatment, levels = c("Sun", "Shade"))
head(sha_data)

sha_data[, c("A4", "B4", "A5", "B5", "AA", "AB", "BB")] <- NULL

m_best <- "
bd ~ c(a, A) * BA
rdry ~ c(b, B) * bd + c(c, C) * BA
h3h1 ~ c(d, D) * bd + c(e, E) * rdry + c(f, F) * BA
idry ~ c(g, G) * bd + c(h, H) * rdry + c(i, I) * h3h1

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
FRIFLC_h3h1_dir_sun := f
FRIFLC_h3h1_dir_shade := F
FRIFLC_h3h1_dir_diff := F - f

FRIFLC_h3h1_ind_sun := (a * d) + (a * b * e) + (c * e)
FRIFLC_h3h1_ind_shade := (A * D) + (A * B * E) + (C * E)
FRIFLC_h3h1_ind_diff := ((A * D) + (A * B * E) + (C * E)) - ((a * d) + (a * b * e) + (c * e))

#idry
FRIFLC_idry_dir_sun := 0
FRIFLC_idry_dir_shade := 0
FRIFLC_idry_dir_diff := 0

FRIFLC_idry_ind_sun := (a * g) + (a * b * h) + (a * d * i) + (a * b * e * i) + (c * h) + (c * e * i) + (f * i)
FRIFLC_idry_ind_shade := (A * G) + (A * B * H) + (A * D * I) + (A * B * E * I) + (C * H) + (C * E * I) + (F * I)
FRIFLC_idry_ind_diff := ((A * G) + (A * B * H) + (A * D * I) + (A * B * E * I) + (C * H) + (C * E * I) + (F * I)) - ((a * g) + (a * b * h) + (a * d * i) + (a * b * e * i) + (c * h) + (c * e * i) + (f * i))
"

fit_best <- sem(m_best, data = sha_data, group = "treatment")
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
path_parameters_eff$pop <- "sha"

setwd("/Users/jkhta/Desktop/nam_cam_fixing/34 - pop_path_analysis/output/")
fwrite(path_parameters_eff, "sha_path_eff_FRIFLC.csv", sep = ",", row.names = FALSE)





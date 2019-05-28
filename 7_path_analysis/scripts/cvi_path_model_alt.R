#path analysis on bur group
library(data.table)
library(lavaan)
library(semPlot)

rm(list = ls())

setwd("/Users/jkhta/Desktop/nam_cam_fixing/34 - pop_path_analysis/output/")

cvi_data <- fread("cvi_col_pheno_qtl_combo_merge.csv",
                  sep = ",",
                  header = TRUE,
                  stringsAsFactors = FALSE)
setnames(cvi_data, c("r_dry", "h3_h1", "i_dry"), c("rdry", "h3h1", "idry"))
cvi_data$treatment <- factor(cvi_data$treatment, levels = c("Sun", "Shade"))
head(cvi_data)

cvi_data[, c("A4", "A5", "AA", "AB", "BA", "BB")] <- NULL

m_best <- "
bd ~ c(a, A) * B5
rdry ~ c(b, B) * bd
h3h1 ~ c(c, C) * bd + c(d, D) * rdry + c(e, E) * B5
idry ~ c(f, F) * bd + c(g, G) * rdry + c(h, H) * h3h1

#B4

#bd
B4_bd_dir_sun := 0
B4_bd_dir_shade := 0
B4_bd_dir_diff := 0

B4_bd_ind_sun := 0
B4_bd_ind_shade := 0
B4_bd_ind_diff := 0

#rdry
B4_rdry_dir_sun := 0
B4_rdry_dir_shade := 0
B4_rdry_dir_diff := 0

B4_rdry_ind_sun := 0
B4_rdry_ind_shade := 0
B4_rdry_ind_diff := 0

#h3h1
B4_h3h1_dir_sun := 0
B4_h3h1_dir_shade := 0
B4_h3h1_dir_diff := 0

B4_h3h1_ind_sun := 0
B4_h3h1_ind_shade := 0
B4_h3h1_ind_diff := 0

#idry
B4_idry_dir_sun := 0
B4_idry_dir_shade := 0
B4_idry_dir_diff := 0

B4_idry_ind_sun := 0
B4_idry_ind_shade := 0
B4_idry_ind_diff := 0

#B5

#bd
B5_bd_dir_sun := a
B5_bd_dir_shade := A
B5_bd_dir_diff := A - a

B5_bd_ind_sun := 0
B5_bd_ind_shade := 0
B5_bd_ind_diff := 0

#rdry
B5_rdry_dir_sun := 0
B5_rdry_dir_shade := 0
B5_rdry_dir_diff := 0

B5_rdry_ind_sun := a * b
B5_rdry_ind_shade := A * B
B5_rdry_ind_diff := (A * B) - (a * b)

#h3h1
B5_h3h1_dir_sun := e
B5_h3h1_dir_shade := E
B5_h3h1_dir_diff := E - e

B5_h3h1_ind_sun := (a * c) + (a * b * d)
B5_h3h1_ind_shade := (A * C) + (A * B * D)
B5_h3h1_ind_diff := ((A * C) + (A * B * D)) - ((a * c) + (a * b * d))

#idry
B5_idry_dir_sun := 0
B5_idry_dir_shade := 0
B5_idry_dir_diff := 0

B5_idry_ind_sun := (a * f) + (a * b * g) + (a * c * h) + (a * b * d * h) + (e * h)
B5_idry_ind_shade := (A * F) + (A * B * G) + (A * C * H) + (A * B * D * H) + (E * H)
B5_idry_ind_diff := ((A * F) + (A * B * G) + (A * C * H) + (A * B * D * H) + (E * H)) - ((a * f) + (a * b * g) + (a * c * h) + (a * b * d * h) + (e * h))
"

fit_best <- sem(m_best, data = cvi_data, group = "treatment")
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
path_parameters_eff$pop <- "cvi"

setwd("/Users/jkhta/Desktop/nam_cam_fixing/34 - pop_path_analysis/output/")
fwrite(path_parameters_eff, "cvi_path_eff_alt.csv", sep = ",", row.names = FALSE)


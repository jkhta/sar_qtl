#path analysis on bur group
library(data.table)
library(lavaan)
library(semPlot)

rm(list = ls())

setwd("/Users/jkhta/Desktop/nam_cam_fixing/34 - pop_path_analysis/output/")

jea_data <- fread("jea_col_pheno_qtl_combo_merge.csv",
                  sep = ",",
                  header = TRUE,
                  stringsAsFactors = FALSE)
setnames(jea_data, c("r_dry", "h3_h1", "i_dry"), c("rdry", "h3h1", "idry"))
jea_data$treatment <- factor(jea_data$treatment, levels = c("Sun", "Shade"))
head(jea_data)

jea_data[, c("A4", "A5", "AA", "AB", "BA", "BB")] <- NULL

m_best <- "
bd ~ c(a, A) * B4
rdry ~ c(b, B) * bd + c(c, C) * B4 + c(d, D) * B5
h3h1 ~ c(e, E) * bd + c(f, F) * rdry
idry ~ c(g, G) * bd + c(h, H) * rdry + c(i, I) * h3h1 + c(j, J) * B4 + c(k, K) * B5

#B4

#bd
B4_bd_dir_sun := a
B4_bd_dir_shade := A
B4_bd_dir_diff := A - a

B4_bd_ind_sun := 0
B4_bd_ind_shade := 0
B4_bd_ind_diff := 0

#rdry
B4_rdry_dir_sun := c
B4_rdry_dir_shade := C
B4_rdry_dir_diff := C - c

B4_rdry_ind_sun := a * b
B4_rdry_ind_shade := A * B
B4_rdry_ind_diff := (A * B) - (a * b)

#h3h1
B4_h3h1_dir_sun := 0
B4_h3h1_dir_shade := 0
B4_h3h1_dir_diff := 0

B4_h3h1_ind_sun := (a * e) + (a * b * f) + (c * f)
B4_h3h1_ind_shade := (A * E) + (A * B * F) + (C * F)
B4_h3h1_ind_diff := ((A * E) + (A * B * F) + (C * F)) - ((a * e) + (a * b * f) + (c * f))

#idry
B4_idry_dir_sun := j
B4_idry_dir_shade := J
B4_idry_dir_diff := J - j

B4_idry_ind_sun := (a * g) + (a * b * h) + (a * e * i) + (a * b * f * i) + (c * h) + (c * f * i)
B4_idry_ind_shade := (A * G) + (A * B * H) + (A * E * I) + (A * B * F * I) + (C * H) + (C * F * I)
B4_idry_ind_diff := ((A * G) + (A * B * H) + (A * E * I) + (A * B * F * I) + (C * H) + (C * F * I)) - ((a * g) + (a * b * h) + (a * e * i) + (a * b * f * i) + (c * h) + (c * f * i))

#B5

#bd
B5_bd_dir_sun := 0
B5_bd_dir_shade := 0
B5_bd_dir_diff := 0

B5_bd_ind_sun := 0
B5_bd_ind_shade := 0
B5_bd_ind_diff := 0

#rdry
B5_rdry_dir_sun := d
B5_rdry_dir_shade := D
B5_rdry_dir_diff := D - d

B5_rdry_ind_sun := 0
B5_rdry_ind_shade := 0
B5_rdry_ind_diff := 0

#h3h1
B5_h3h1_dir_sun := 0
B5_h3h1_dir_shade := 0
B5_h3h1_dir_diff := 0

B5_h3h1_ind_sun := d * f
B5_h3h1_ind_shade := D * F
B5_h3h1_ind_diff := (D * F) - (d * f)

#idry
B5_idry_dir_sun := k
B5_idry_dir_shade := K
B5_idry_dir_diff := K - k

B5_idry_ind_sun := (d * h) + (d * f * i)
B5_idry_ind_shade := (D * H) + (D * F * I)
B5_idry_ind_diff := ((D * H) + (D * F * I)) - ((d * h) + (d * f * i))
"

fit_best <- sem(m_best, data = jea_data, group = "treatment")
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
path_parameters_eff$pop <- "jea"

setwd("/Users/jkhta/Desktop/nam_cam_fixing/34 - pop_path_analysis/output/")
fwrite(path_parameters_eff, "jea_path_eff_alt.csv", sep = ",", row.names = FALSE)



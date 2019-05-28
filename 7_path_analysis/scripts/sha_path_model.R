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

m_best <- "
bd ~ c(a, A) * AB + c(b, B) * BA
rdry ~ c(c, C) * bd + c(d, D) * BA + c(e, E) * BB
h3h1 ~ c(f, F) * bd + c(g, G) * rdry + c(h, H) * AB + c(i, I) * BB
idry ~ c(j, J) * bd + c(k, K) * rdry + c(l, L) * h3h1

#AB
#bd 
AB_bd_dir_sun := a
AB_bd_dir_shade := A
AB_bd_dir_diff := A - a

AB_bd_ind_sun := 0
AB_bd_ind_shade := 0
AB_bd_ind_diff := 0

#r_dry
AB_rdry_dir_sun := 0
AB_rdry_dir_shade := 0
AB_rdry_dir_diff := 0

AB_rdry_ind_sun := a * c
AB_rdry_ind_shade := A * C
AB_rdry_ind_diff := (A * C) - (a * c)

#h3_h1
AB_h3h1_dir_sun := h
AB_h3h1_dir_shade := H
AB_h3h1_dir_diff := H - h

AB_h3h1_ind_sun := (a * f) + (a * c * g)
AB_h3h1_ind_shade := (A * F) + (A * C * G)
AB_h3h1_ind_diff := ((A * F) + (A * C * G)) - ((a * f) + (a * c * g))

#i_dry
AB_idry_dir_sun := 0
AB_idry_dir_shade := 0
AB_idry_dir_diff := 0

AB_idry_ind_sun := (a * j) + (a * c * k) + (a * f * l) + (a * c * g * l) + (h * l)
AB_idry_ind_shade := (A * J) + (A * C * K) + (A * F * L) + (A * C * G * L) + (H * L)
AB_idry_ind_diff := ((A * J) + (A * C * K) + (A * F * L) + (A * C * G * L) + (H * L)) - ((a * j) + (a * c * k) + (a * f * l) + (a * c * g * l) + (h * l))

#BA
#bd
BA_bd_dir_sun := b
BA_bd_dir_shade := B
BA_bd_dir_diff := B - b

BA_bd_ind_sun := 0
BA_bd_ind_shade := 0
BA_bd_ind_diff := 0

#r_dry
BA_rdry_dir_sun := d
BA_rdry_dir_shade := D
BA_rdry_dir_diff := D - d

BA_rdry_ind_sun := b * c
BA_rdry_ind_shade := B * C
BA_rdry_ind_diff := (B * C) - (b * c)

#h3_h1
BA_h3h1_dir_sun := 0
BA_h3h1_dir_shade := 0
BA_h3h1_dir_diff := 0

BA_h3h1_ind_sun := (b * f) + (b * c * g) + (d * g)
BA_h3h1_ind_shade := (B * F) + (B * C * G) + (D * G)
BA_h3h1_ind_diff := ((B * F) + (B * C * G) + (D * G)) - ((b * f) + (b * c * g) + (d * g))

#i_dry
BA_idry_dir_sun := 0
BA_idry_dir_shade := 0
BA_idry_dir_diff := 0

BA_idry_ind_sun := (b * j) + (b * c * k) + (b * f * l) + (b * c * g * l) + (d * k) + (d * g * l)
BA_idry_ind_shade := (B * K) + (B * C * K) + (B * F * L) + (B * C * G * L) + (D * K) + (D * G * L)
BA_idry_ind_diff := ((B * K) + (B * C * K) + (B * F * L) + (B * C * G * L) + (D * K) + (D * G * L)) - ((b * j) + (b * c * k) + (b * f * l) + (b * c * g * l) + (d * k) + (d * g * l))

#BB
#bd 
BB_bd_dir_sun := 0
BB_bd_dir_shade := 0
BB_bd_dir_diff := 0

BB_bd_ind_sun := 0
BB_bd_ind_shade := 0
BB_bd_ind_diff := 0

#r_dry
BB_rdry_dir_sun := e
BB_rdry_dir_shade := E
BB_rdry_dir_diff := E - e

BB_rdry_ind_sun := 0
BB_rdry_ind_shade := 0
BB_rdry_ind_diff := 0

#h3_h1
BB_h3h1_dir_sun := i
BB_h3h1_dir_shade := I
BB_h3h1_dir_diff := I - i

BB_h3h1_ind_sun := e * g
BB_h3h1_ind_shade := E * G
BB_h3h1_ind_diff := (E * G) - (e * g)

#i_dry
BB_idry_dir_sun := 0
BB_idry_dir_shade := 0
BB_idry_dir_diff := 0

BB_idry_ind_sun := (e * k) + (e * g * l) + (i * l)
BB_idry_ind_shade := (E * K) + (E * G * L) + (I * L)
BB_idry_ind_diff := ((E * K) + (E * G * L) + (I * L)) - ((e * k) + (e * g * l) + (i * l))
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
fwrite(path_parameters_eff, "sha_path_eff.csv", sep = ",", row.names = FALSE)





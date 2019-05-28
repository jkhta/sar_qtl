#path analysis on bur group
library(data.table)
library(lavaan)
library(semPlot)

rm(list = ls())

setwd("/Users/jkhta/Desktop/nam_cam_fixing/34 - pop_path_analysis/output/")

ita_data <- fread("ita_col_pheno_qtl_combo_merge.csv",
                  sep = ",",
                  header = TRUE,
                  stringsAsFactors = FALSE)
setnames(ita_data, c("r_dry", "h3_h1", "i_dry"), c("rdry", "h3h1", "idry"))
ita_data$treatment <- factor(ita_data$treatment, levels = c("Sun", "Shade"))
head(ita_data)

m_best <- "
bd ~ c(a, A) * AB + c(b, B) * BA
rdry ~ c(c, C) * bd + c(d, D) * AB
h3h1 ~ c(e, E) * bd + c(f, F) * rdry + c(g, G) * AB + c(h, H) * BA
idry ~ c(i, I) * bd + c(j, J) * rdry + c(k, K) * h3h1 + c(l, L) * BA

#AB
#bd
AB_bd_dir_sun := a
AB_bd_dir_shade := A
AB_bd_dir_diff := A - a

AB_bd_ind_sun := 0
AB_bd_ind_shade := 0
AB_bd_ind_diff := 0

#r_dry
AB_rdry_dir_sun := d
AB_rdry_dir_shade := D
AB_rdry_dir_diff := D - d

AB_rdry_ind_sun := a * c
AB_rdry_ind_shade := A * C
AB_rdry_ind_diff := (A * C) - (a * c)

#h3_h1
AB_h3h1_dir_sun := g
AB_h3h1_dir_shade := G
AB_h3h1_dir_diff := G - g

AB_h3h1_ind_sun := (a * e) + (a * c * f) + (d * f)
AB_h3h1_ind_shade := (A * E) + (A * C * F) + (D * F)
AB_h3h1_ind_diff := ((A * E) + (A * C * F) + (D * F)) - ((a * e) + (a * c * f) + (d * f))

#i_dry
AB_idry_dir_sun := 0
AB_idry_dir_shade := 0
AB_idry_dir_diff := 0

AB_idry_ind_sun := (a * i) + (a * c * j) + (a * e * k) + (a * c * f * k) + (d * j) + (d * f * k) + (g * k)
AB_idry_ind_shade := (A * I) + (A * C * J) + (A * E * K) + (A * C * F * K) + (D * J) + (D * F * K) + (G * K)
AB_idry_ind_diff := ((A * I) + (A * C * J) + (A * E * K) + (A * C * F * K) + (D * J) + (D * F * K) + (G * K)) - ((a * i) + (a * c * j) + (a * e * k) + (a * c * f * k) + (d * j) + (d * f * k) + (g * k))

#BA
#bd
BA_bd_dir_sun := b
BA_bd_dir_shade := B
BA_bd_dir_diff := B - b

BA_bd_ind_sun := 0
BA_bd_ind_shade := 0
BA_bd_ind_diff := 0

#r_dry
BA_rdry_dir_sun := 0
BA_rdry_dir_shade := 0
BA_rdry_dir_diff := 0

BA_rdry_ind_sun := b * c
BA_rdry_ind_shade := B * C
BA_rdry_ind_diff := (B * C) - (b * c)

#h3_h1
BA_h3h1_dir_sun := h
BA_h3h1_dir_shade := H
BA_h3h1_dir_diff := H - h

BA_h3h1_ind_sun := (b * e) + (b * c * f)
BA_h3h1_ind_shade := (B * E) + (B * C * F)
BA_h3h1_ind_diff := ((B * E) + (B * C * F)) - ((b * e) + (b * c * f))

#i_dry
BA_idry_dir_sun := l
BA_idry_dir_shade := L
BA_idry_dir_diff := L - l

BA_idry_ind_sun := (b * i) + (b * c * j) + (b * e * k) + (b * c * f * k) + (h * k)
BA_idry_ind_shade := (B * I) + (B * C * J) + (B * E * K) + (B * C * F * K) + (H * K)
BA_idry_ind_diff := ((B * I) + (B * C * J) + (B * E * K) + (B * C * F * K) + (H * K)) - ((b * i) + (b * c * j) + (b * e * k) + (b * c * f * k) + (h * k))

#BB
#bd
BB_bd_dir_sun := 0
BB_bd_dir_shade := 0
BB_bd_dir_diff := 0

BB_bd_ind_sun := 0
BB_bd_ind_shade := 0
BB_bd_ind_diff := 0

#r_dry
BB_rdry_dir_sun := 0
BB_rdry_dir_shade := 0
BB_rdry_dir_diff := 0

BB_rdry_ind_sun := 0
BB_rdry_ind_shade := 0
BB_rdry_ind_diff := 0

#h3_h1
BB_h3h1_dir_sun := 0
BB_h3h1_dir_shade := 0
BB_h3h1_dir_diff := 0

BB_h3h1_ind_sun := 0
BB_h3h1_ind_shade := 0
BB_h3h1_ind_diff := 0

#i_dry
BB_idry_dir_sun := 0
BB_idry_dir_shade := 0
BB_idry_dir_diff := 0

BB_idry_ind_sun := 0
BB_idry_ind_shade := 0
BB_idry_ind_diff := 0
"

fit_best <- sem(m_best, data = ita_data, group = "treatment")
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
path_parameters_eff$pop <- "ita"

setwd("/Users/jkhta/Desktop/nam_cam_fixing/34 - pop_path_analysis/output/")
fwrite(path_parameters_eff, "ita_path_eff.csv", sep = ",", row.names = FALSE)


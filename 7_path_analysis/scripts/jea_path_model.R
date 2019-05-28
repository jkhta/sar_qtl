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

m_best <- "
bd ~ c(a, A) * BA + c(b, B) * BB
rdry ~ c(c, C) * bd + c(d, D) * BA
h3h1 ~ c(e, E) * bd + c(f, F) * rdry
idry ~ c(g, G) * bd + c(h, H) * rdry + c(i, I) * h3h1 + c(j, J) * AB + c(k, K) * BA + c(l, L) * BB

#BA
#bd
BA_bd_dir_sun := a
BA_bd_dir_shade := A
BA_bd_dir_diff := A - a

BA_bd_ind_sun := 0
BA_bd_ind_shade := 0
BA_bd_ind_diff := 0

#r_dry
BA_rdry_dir_sun := d
BA_rdry_dir_shade := D
BA_rdry_dir_diff := D - d

BA_rdry_ind_sun := a * c
BA_rdry_ind_shade := A * C
BA_rdry_ind_diff := (A * C) - (a * c)

#h3_h1
BA_h3h1_dir_sun := 0
BA_h3h1_dir_shade := 0
BA_h3h1_dir_diff := 0

BA_h3h1_ind_sun := (a * e) + (a * c * f) + (d * f)
BA_h3h1_ind_shade := (A * E) + (A * C * F) + (D * F)
BA_h3h1_ind_diff := ((A * E) + (A * C * F) + (D * F)) - ((a * e) + (a * c * f) + (d * f))

#i_dry
BA_idry_dir_sun := k
BA_idry_dir_shade := K
BA_idry_dir_diff := K - k

BA_idry_ind_sun := (a * g) + (a * c * h) + (a * e * i) + (a * c * f * i) + (d * h) + (d * f * i)
BA_idry_ind_shade := (A * G) + (A * C * H) + (A * E * I) + (D * H) + (D * F * I)
BA_idry_ind_diff := ((A * G) + (A * C * H) + (A * E * I) + (D * H) + (D * F * I)) - ((a * g) + (a * c * h) + (a * e * i) + (d * h) + (d * f * i))

#BB
#bd
BB_bd_dir_sun := b
BB_bd_dir_shade := B
BB_bd_dir_diff := B - b

BB_bd_ind_sun := 0
BB_bd_ind_shade := 0
BB_bd_ind_diff := 0

#r_dry
BB_rdry_dir_sun := 0
BB_rdry_dir_shade := 0
BB_rdry_dir_diff := 0

BB_rdry_ind_sun := b * c
BB_rdry_ind_shade := B * C
BB_rdry_ind_diff := (B * C) - (b * c)

#h3_h1
BB_h3h1_dir_sun := 0
BB_h3h1_dir_shade := 0
BB_h3h1_dir_diff := 0

BB_h3h1_ind_sun := (b * e) + (b * c * f)
BB_h3h1_ind_shade := (B * E) + (B * C * F)
BB_h3h1_ind_diff := ((B * E) + (B * C * F)) - ((b * e) + (b * c * f))

#i_dry
BB_idry_dir_sun := l
BB_idry_dir_shade := L
BB_idry_dir_diff := L - l

BB_idry_ind_sun := (b * g) + (b * c * h) + (b * e * i) + (b * c * f * i) 
BB_idry_ind_shade := (B * G) + (B * C * H) + (B * E * I) + (B * C * F * I)
BB_idry_ind_diff := ((B * G) + (B * C * H) + (B * E * I) + (B * C * F * I)) - ((b * g) + (b * c * h) + (b * e * i) + (b * c * f * i))

#AB
#bd
AB_bd_dir_sun := 0
AB_bd_dir_shade := 0
AB_bd_dir_diff := 0

AB_bd_ind_sun := 0
AB_bd_ind_shade := 0
AB_bd_ind_diff := 0

#r_dry
AB_rdry_dir_sun := 0
AB_rdry_dir_shade := 0
AB_rdry_dir_diff := 0

AB_rdry_ind_sun := 0
AB_rdry_ind_shade := 0
AB_rdry_ind_diff := 0

#h3_h1
AB_h3h1_dir_sun := 0
AB_h3h1_dir_shade := 0
AB_h3h1_dir_diff := 0

AB_h3h1_ind_sun := 0
AB_h3h1_ind_shade := 0
AB_h3h1_ind_diff := 0

#i_dry
AB_idry_dir_sun := j
AB_idry_dir_shade := J
AB_idry_dir_diff := J - j

AB_idry_ind_sun := 0
AB_idry_ind_shade := 0
AB_idry_ind_diff := 0
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
fwrite(path_parameters_eff, "jea_path_eff.csv", sep = ",", row.names = FALSE)



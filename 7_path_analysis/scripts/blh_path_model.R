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

m_best <- "
bd ~ c(a, A) * BA + c(b, B) * BB
rdry ~ c(c, C) * bd + c(d, D) * AB + c(e, E) * BA + c(f, F) * BB
h3h1 ~ c(g, G) * AB
idry ~ c(h, H) * bd + c(i, I) * rdry + c(j, J) * h3h1 + c(k, K) * BA

#BA
#bd
BA_bd_dir_sun := a
BA_bd_dir_shade := A
BA_bd_dir_diff := A - a

BA_bd_ind_sun := 0
BA_bd_ind_shade := 0
BA_bd_ind_diff := 0

#r_dry
BA_rdry_dir_sun := e
BA_rdry_dir_shade := E
BA_rdry_dir_diff := E - e

BA_rdry_ind_sun := a * c
BA_rdry_ind_shade := A * C
BA_rdry_ind_diff := (A * C) - (a * c)

#h3_h1
BA_h3h1_dir_sun := 0
BA_h3h1_dir_shade := 0
BA_h3h1_dir_diff := 0

BA_h3h1_ind_sun := 0
BA_h3h1_ind_shade := 0
BA_h3h1_ind_diff := 0

#i_dry
BA_idry_dir_sun := k
BA_idry_dir_shade := K
BA_idry_dir_diff := K - k

BA_idry_ind_sun := (a * h) + (a * c * i) + (e * i)
BA_idry_ind_shade := (A * H) + (A * C * I) + (E * I)
BA_idry_ind_diff := ((A * H) + (A * C * I) + (E * I)) - ((a * h) + (a * c * i) + (e * i))

#BB
#bd
BB_bd_dir_sun := b
BB_bd_dir_shade := B
BB_bd_dir_diff := B - b

BB_bd_ind_sun := 0
BB_bd_ind_shade := 0
BB_bd_ind_diff := 0

#r_dry
BB_rdry_dir_sun := f
BB_rdry_dir_shade := F
BB_rdry_dir_diff := F - f

BB_rdry_ind_sun := b * c
BB_rdry_ind_shade := B * C
BB_rdry_ind_diff := (B * C) - (b * c)

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

BB_idry_ind_sun := (b * h) + (b * c * i) + (f * i)
BB_idry_ind_shade := (B * H) + (B * C * I) + (F * I)
BB_idry_ind_diff := ((B * H) + (B * C * I) + (F * I)) - ((B * H) + (B * C * I) + (F * I))

#AB
#bd
AB_bd_dir_sun := 0
AB_bd_dir_shade := 0
AB_bd_dir_diff := 0

AB_bd_ind_sun := 0
AB_bd_ind_shade := 0
AB_bd_ind_diff := 0

#r_dry
AB_rdry_dir_sun := d
AB_rdry_dir_shade := D
AB_rdry_dir_diff := D - d

AB_rdry_ind_sun := 0
AB_rdry_ind_shade := 0
AB_rdry_ind_diff := 0

#h3_h1
AB_h3h1_dir_sun := g
AB_h3h1_dir_shade := G
AB_h3h1_dir_diff := G - g

AB_h3h1_ind_sun := 0
AB_h3h1_ind_shade := 0
AB_h3h1_ind_diff := 0

#i_dry
AB_idry_dir_sun := 0
AB_idry_dir_shade := 0
AB_idry_dir_diff := 0

AB_idry_ind_sun := (d * i) + (g * j)
AB_idry_ind_shade := (D * I) + (G * J)
AB_idry_ind_diff := ((D * I) + (G * J)) - ((d * i) + (g * j))
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
fwrite(path_parameters_eff, "blh_path_eff.csv", sep = ",", row.names = FALSE)

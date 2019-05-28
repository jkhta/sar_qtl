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

m_best <- "
bd ~ AB + BA + BB
rdry ~ c(a, A) * bd + c(b, B) * AB + c(c, C) * BB
h3h1 ~ c(d, D) * bd + c(e, E) * rdry + c(f, F) * AB + c(g, G) * BB
idry ~ c(h, H) * bd + c(i, I) * rdry + c(j, J) * h3h1 + c(h, H) * BA

#AB
#bd 
AB_bd_dir_sun := 0
AB_bd_dir_shade := 0
AB_bd_dir_diff := 0

AB_bd_ind_sun := 0
AB_bd_ind_shade := 0
AB_bd_ind_diff := 0

#r_dry
AB_rdry_dir_sun := b
AB_rdry_dir_shade := B
AB_rdry_dir_diff := B - b

AB_rdry_ind_sun := 0
AB_rdry_ind_shade := 0
AB_rdry_ind_diff := 0

#h3_h1
AB_h3h1_dir_sun := f
AB_h3h1_dir_shade := F
AB_h3h1_dir_diff := F - f

AB_h3h1_ind_sun := b * e
AB_h3h1_ind_shade := B * E
AB_h3h1_ind_diff := (B * E) - (b * e)

#i_dry
AB_idry_dir_sun := 0
AB_idry_dir_shade := 0
AB_idry_dir_diff := 0

AB_idry_ind_sun := (b * i) + (b * e * j) + (f * j)
AB_idry_ind_shade := (B * I) + (B * E * J) + (F * J)
AB_idry_ind_diff := ((B * I) + (B * E * J) + (F * J)) - ((b * i) + (b * e * j) + (f * j))

#BB
#bd
BB_bd_dir_sun := 0
BB_bd_dir_shade := 0
BB_bd_dir_diff := 0

BB_bd_ind_sun := 0
BB_bd_ind_shade := 0
BB_bd_ind_diff := 0

#r_dry
BB_rdry_dir_sun := c
BB_rdry_dir_shade := C
BB_rdry_dir_diff := C - c

BB_rdry_ind_sun := 0
BB_rdry_ind_shade := 0
BB_rdry_ind_diff := 0

#h3_h1
BB_h3h1_dir_sun := g
BB_h3h1_dir_shade := G
BB_h3h1_dir_diff := G - g

BB_h3h1_ind_sun := c * e
BB_h3h1_ind_shade := C * E
BB_h3h1_ind_diff := (C * E) - (c * e)

#i_dry
BB_idry_dir_sun := 0
BB_idry_dir_shade := 0
BB_idry_dir_diff := 0

BB_idry_ind_sun := (c * i) + (c * e * j) + (g * j)
BB_idry_ind_shade := (C * I) + (C * E * J) + (G * J)
BB_idry_ind_diff := ((C * I) + (C * E * J) + (G * J)) - ((c * i) + (c * e * j) + (g * j))

#BA
#bd
BA_bd_dir_sun := 0
BA_bd_dir_shade := 0
BA_bd_dir_diff := 0

BA_bd_ind_sun := 0
BA_bd_ind_shade := 0
BA_bd_ind_diff := 0

#r_dry
BA_rdry_dir_sun := 0
BA_rdry_dir_shade := 0
BA_rdry_dir_diff := 0

BA_rdry_ind_sun := 0
BA_rdry_ind_shade := 0
BA_rdry_ind_diff := 0

#h3_h1
BA_h3h1_dir_sun := 0
BA_h3h1_dir_shade := 0
BA_h3h1_dir_diff := 0

BA_h3h1_ind_sun := 0
BA_h3h1_ind_shade := 0
BA_h3h1_ind_diff := 0

#i_dry
BA_idry_dir_sun := h
BA_idry_dir_shade := H
BA_idry_dir_diff := H - h

BA_idry_ind_sun := 0
BA_idry_ind_shade := 0
BA_idry_ind_diff := 0
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
fwrite(path_parameters_eff, "bur_path_eff.csv", sep = ",", row.names = FALSE)

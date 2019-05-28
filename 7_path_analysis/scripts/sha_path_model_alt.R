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

sha_data[, c("A4", "A5", "AA", "AB", "BA", "BB")] <- NULL

m_best <- "
bd ~ c(a, A) * B4 + c(b, B) * B5
rdry ~ c(c, C) * bd + c(d, D) * B4
h3h1 ~ c(e, E) * bd + c(f, F) * rdry + c(g, G) * B5
idry ~ c(h, H) * bd + c(i, I) * rdry + c(j, J) * h3h1

#B4

#bd
B4_bd_dir_sun := a
B4_bd_dir_shade := A
B4_bd_dir_diff := A - a

B4_bd_ind_sun := 0
B4_bd_ind_shade := 0
B4_bd_ind_diff := 0

#rdry
B4_rdry_dir_sun := d
B4_rdry_dir_shade := D
B4_rdry_dir_diff := D - d

B4_rdry_ind_sun := a * c
B4_rdry_ind_shade := A * C
B4_rdry_ind_diff := (A * C) - (a * c)

#h3h1
B4_h3h1_dir_sun := 0
B4_h3h1_dir_shade := 0
B4_h3h1_dir_diff := 0

B4_h3h1_ind_sun := (a * e) + (a * c * f) + (d * f)
B4_h3h1_ind_shade := (A * E) + (A * C * F) + (D * F)
B4_h3h1_ind_diff := ((A * E) + (A * C * F) + (D * F)) - ((a * e) + (a * c * f) + (d * f))

#idry
B4_idry_dir_sun := 0
B4_idry_dir_shade := 0
B4_idry_dir_diff := 0

B4_idry_ind_sun := (a * h) + (a * c * i) + (a * e * j) + (a * c * f * j) + (d * i) + (d * f * j)
B4_idry_ind_shade := (A * H) + (A * C * I) + (A * E * J) + (A * C * F * J) + (D * I) + (D * F * J)
B4_idry_ind_diff := ((A * H) + (A * C * I) + (A * E * J) + (A * C * F * J) + (D * I) + (D * F * J)) - ((a * h) + (a * c * i) + (a * e * j) + (a * c * f * j) + (d * i) + (d * f * j))

#B5

#bd
B5_bd_dir_sun := b
B5_bd_dir_shade := B
B5_bd_dir_diff := B - b

B5_bd_ind_sun := 0
B5_bd_ind_shade := 0
B5_bd_ind_diff := 0

#rdry
B5_rdry_dir_sun := 0
B5_rdry_dir_shade := 0
B5_rdry_dir_diff := 0

B5_rdry_ind_sun := b * c
B5_rdry_ind_shade := B * C
B5_rdry_ind_diff := (B * C) - (b * c)

#h3h1
B5_h3h1_dir_sun := g
B5_h3h1_dir_shade := G
B5_h3h1_dir_diff := G - g

B5_h3h1_ind_sun := (b * e) + (b * c * f) 
B5_h3h1_ind_shade := (B * E) + (B * C * F) 
B5_h3h1_ind_diff := ((B * E) + (B * C * F)) - ((b * e) + (b * c * f))

#idry
B5_idry_dir_sun := 0
B5_idry_dir_shade := 0
B5_idry_dir_diff := 0

B5_idry_ind_sun := (b * h) + (b * c * i) + (b * e * j) + (b * c * f * j) + (g * j)
B5_idry_ind_shade := (B * H) + (B * C * I) + (B * E * J) + (B * C * F * J) + (G * J)
B5_idry_ind_diff := ((B * H) + (B * C * I) + (B * E * J) + (B * C * F * J) + (G * J)) - ((b * h) + (b * c * i) + (b * e * j) + (b * c * f * j) + (g * j))
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
fwrite(path_parameters_eff, "sha_path_eff_alt.csv", sep = ",", row.names = FALSE)





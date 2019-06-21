#path analysis on bur group
library(data.table)
library(lavaan)
library(semPlot)
library(lettercase)

rm(list = ls())

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/7_path_analysis/output/")

ita_data <- fread("29RV_blups_qtl_combo_merge.csv",
                  sep = ",",
                  header = TRUE,
                  stringsAsFactors = FALSE)
setnames(ita_data, c("bd_geno", "h3_h1_geno", "i_dry_geno", "r_dry_geno"), c("bd", "h3h1", "idry", "rdry"))
head(ita_data)

ita_data[, c("A4", "A42", "A5", "AA", "AB", "BA", "BB")] <- NULL

m_best <- "
bd ~ B4 + B5
#rdry ~ bd + B42 + B5
"

fit_best <- sem(m_best, data = ita_data)
summary(fit_best, fit.measures = TRUE)


ita_data <- fread("ita_col_pheno_qtl_combo_merge.csv",
                  sep = ",",
                  header = TRUE,
                  stringsAsFactors = FALSE)
setnames(ita_data, c("r_dry", "h3_h1", "i_dry"), c("rdry", "h3h1", "idry"))
ita_data$treatment <- factor(ita_data$treatment, levels = c("Sun", "Shade"))
head(ita_data)

ita_data[, c("A4", "A42", "A5", "AA", "AB", "BA", "BB")] <- NULL

m_best <- "
bd ~ B4 + B5
"

fit_best <- sem(m_best, data = ita_data, group = "treatment")
summary(fit_best, fit.measures = TRUE)

m_best <- "
bd ~ c(A, a) * B4 + c(B, b) * B5
rdry ~ c(C, c) * bd + c(D, d) * B42 + c(E, e) * B5
h3h1 ~ c(F, f) * bd + c(G, g) * rdry + c(H, h) * B4 + c(I, i) * B5
idry ~ c(J, j) * bd + c(K, k) * rdry + c(L, l) * h3h1 + c(M, m) * B4

#B4

#bd
B4_bd_dir_sun := a
B4_bd_dir_shade := A
B4_bd_dir_diff := A - a

B4_bd_ind_sun := 0
B4_bd_ind_shade := 0
B4_bd_ind_diff := 0

#rdry
B4_rdry_dir_sun := 0
B4_rdry_dir_shade := 0
B4_rdry_dir_diff := 0

B4_rdry_ind_sun := a * c
B4_rdry_ind_shade := A * C
B4_rdry_ind_diff := (A * C) - (a * c)

#h3h1
B4_h3h1_dir_sun := h
B4_h3h1_dir_shade := H
B4_h3h1_dir_diff := H - h

B4_h3h1_ind_sun := (a * f) + (a * c  *g)
B4_h3h1_ind_shade := (A * F) + (A * C * G)
B4_h3h1_ind_diff := ((A * F) + (A * C * G)) - ((a * f) + (a * c  *g))

#idry
B4_idry_dir_sun := m
B4_idry_dir_shade := M
B4_idry_dir_diff := M - m

B4_idry_ind_sun := (a * j) + (a * c * k) + (a * f * l) + (a * c * g * l)
B4_idry_ind_shade := (A * J) + (A * C * K) + (A * F * L) + (A * C * G * L)
B4_idry_ind_diff := ((A * J) + (A * C * K) + (A * F * L) + (A * C * G * L)) - ((a * j) + (a * c * k) + (a * f * l) + (a * c * g * l))

#B42

#bd
B42_bd_dir_sun := 0
B42_bd_dir_shade := 0
B42_bd_dir_diff := 0

B42_bd_ind_sun := 0
B42_bd_ind_shade := 0
B42_bd_ind_diff := 0

#rdry
B42_rdry_dir_sun := d
B42_rdry_dir_shade := D
B42_rdry_dir_diff := D - d

B42_rdry_ind_sun := 0
B42_rdry_ind_shade := 0
B42_rdry_ind_diff := 0

#h3h1
B42_h3h1_dir_sun := 0
B42_h3h1_dir_shade := 0
B42_h3h1_dir_diff := 0

B42_h3h1_ind_sun := d * g
B42_h3h1_ind_shade := (D * G)
B42_h3h1_ind_diff := (D * G) - (d * g)

#idry
B42_idry_dir_sun := 0
B42_idry_dir_shade := 0
B42_idry_dir_diff := 0

B42_idry_ind_sun := (d * k) + (d * g * l)
B42_idry_ind_shade := (D * K) + (D * G * L)
B42_idry_ind_diff := ((D * K) + (D * G * L)) - ((d * k) + (d * g * l))

#B5

#bd
B5_bd_dir_sun := b
B5_bd_dir_shade := B
B5_bd_dir_diff := B - b

B5_bd_ind_sun := 0
B5_bd_ind_shade := 0
B5_bd_ind_diff := 0

#rdry
B5_rdry_dir_sun := e
B5_rdry_dir_shade := E
B5_rdry_dir_diff := E - e

B5_rdry_ind_sun := b * c
B5_rdry_ind_shade := B * C
B5_rdry_ind_diff := (B * C) - (b * c)

#h3h1
B5_h3h1_dir_sun := i
B5_h3h1_dir_shade := I
B5_h3h1_dir_diff := I - i

B5_h3h1_ind_sun := (b * f) + (b * c * g) 
B5_h3h1_ind_shade := (B * F) + (B * C * G)
B5_h3h1_ind_diff := ((B * F) + (B * C * G)) - ((b * f) + (b * c * g))

#idry
B5_idry_dir_sun := 0
B5_idry_dir_shade := 0
B5_idry_dir_diff := 0

B5_idry_ind_sun := (b * j) + (b * c * k) + (b * f * l) + (b * c * g * l) + (e * k) + (e * g * l) + (i * l)
B5_idry_ind_shade := (B * J) + (B * C * K) + (B * F * L) + (B * C * G * L) + (E * K) + (E * G * L) + (I * L)
B5_idry_ind_diff := ((B * J) + (B * C * K) + (B * F * L) + (B * C * G * L) + (E * K) + (E * G * L) + (I * L)) - ((b * j) + (b * c * k) + (b * f * l) + (b * c * g * l) + (e * k) + (e * g * l) + (i * l))
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

path_parameters_qtl <- subset(path_parameters, lhs %in% c("bd", "rdry", "h3h1", "idry") & rhs %in% c("B4", "B42", "B5") & !grepl("~~", op), select = c(lhs, op, rhs, group, est, label))
path_parameters_qtl$formula <- with(path_parameters_qtl, paste(lhs, op, rhs, sep = " "))
path_parameters_qtl_sun <- path_parameters_qtl[str_is(path_parameters_qtl$label, str_lower_case), ]
required_formula <- c("bd ~ B4", "bd ~ B42", "bd ~ B5", "rdry ~ B4", "rdry ~ B42", "rdry ~ B5", "h3h1 ~ B4", "h3h1 ~ B42", "h3h1 ~ B5", "idry ~ B4", "idry ~ B42", "idry ~ B5") 
path_parameters_qtl_effects_sun <- data.frame(matrix(ncol = length(required_formula)))

for (i in 1:length(required_formula)) {
  path_parameters_qtl_effects_sun[, i] <- tryCatch(subset(path_parameters_qtl_sun, formula == required_formula[i])$est, error = NA)
}
colnames(path_parameters_qtl_effects_sun) <- required_formula
path_parameters_qtl_effects_sun$env <- "Sun"
path_parameters_qtl_effects_sun$pop <- "ita"

path_parameters_qtl_shade <- path_parameters_qtl[str_is(path_parameters_qtl$label, str_upper_case), ]
path_parameters_qtl_effects_shade <- data.frame(matrix(ncol = length(required_formula)))

for (i in 1:length(required_formula)) {
  path_parameters_qtl_effects_shade[, i] <- tryCatch(subset(path_parameters_qtl_shade, formula == required_formula[i])$est, error = NA)
}
colnames(path_parameters_qtl_effects_shade) <- required_formula
path_parameters_qtl_effects_shade$env <- "Shade"
path_parameters_qtl_effects_shade$pop <- "ita"

path_parameters_trait <- subset(path_parameters, lhs %in% c("bd", "rdry", "h3h1", "idry") & rhs %in% c("bd", "rdry", "h3h1", "idry") & !grepl("~~", op), select = c(lhs, op, rhs, group, est, label))
path_parameters_trait$formula <- with(path_parameters_trait, paste(lhs, op, rhs, sep = " "))
path_parameters_trait_sun <- path_parameters_trait[str_is(path_parameters_trait$label, str_lower_case), ]
required_formula <- c("rdry ~ bd", "h3h1 ~ bd", "h3h1 ~ rdry", "idry ~ bd", "idry ~ rdry", "idry ~ h3h1")
path_trait_correlations_sun <- data.frame(matrix(ncol = length(required_formula)))

for (i in 1:length(required_formula)) {
  path_trait_correlations_sun[, i] <- tryCatch(subset(path_parameters_trait_sun, formula == required_formula[i])$est, error = NA)
}
colnames(path_trait_correlations_sun) <- required_formula
path_trait_correlations_sun$env <- "Sun"
path_trait_correlations_sun$pop <- "ita"

path_parameters_trait_shade <- path_parameters_trait[str_is(path_parameters_trait$label, str_upper_case), ]
path_trait_correlations_shade <- data.frame(matrix(ncol = length(required_formula)))

for (i in 1:length(required_formula)) {
  path_trait_correlations_shade[, i] <- tryCatch(subset(path_parameters_trait_shade, formula == required_formula[i])$est, error = NA)
}
colnames(path_trait_correlations_shade) <- required_formula
path_trait_correlations_shade$env <- "Shade"
path_trait_correlations_shade$pop <- "ita"

path_parameters_eff <- subset(path_parameters, grepl("sun|shade|diff", lhs))
path_parameters_eff$allele_comb <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[1])
path_parameters_eff$trait <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[2])
path_parameters_eff$effect_type <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[3])
path_parameters_eff$env <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[4])
nrow(path_parameters_eff)
path_parameters_eff$pop <- "ita"

setwd("/Users/James/Documents/GitHub/sar_qtl/7_path_analysis/output/")
fwrite(path_parameters_eff, "ita_path_eff_alt.csv", sep = ",", row.names = FALSE)

fwrite(path_trait_correlations_sun, "ita_trait_eff_alt_sun.csv", sep = ",", row.names = FALSE, na = "NA")
fwrite(path_trait_correlations_shade, "ita_trait_eff_alt_shade.csv", sep = ",", row.names = FALSE, na = "NA")

fwrite(path_parameters_qtl_effects_sun, "ita_qtl_eff_alt_sun.csv", sep = ",", row.names = FALSE, na = "NA")
fwrite(path_parameters_qtl_effects_shade, "ita_qtl_eff_alt_shade.csv", sep = ",", row.names = FALSE, na = "NA")

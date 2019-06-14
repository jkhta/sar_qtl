#path analysis on bur group
library(data.table)
library(lavaan)
library(semPlot)
library(lettercase)

rm(list = ls())

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/7_path_analysis/output/")

bur_data <- fread("bur_col_pheno_qtl_combo_merge.csv",
                  sep = ",",
                  header = TRUE,
                  stringsAsFactors = FALSE)
setnames(bur_data, c("r_dry", "h3_h1", "i_dry"), c("rdry", "h3h1", "idry"))
bur_data$treatment <- factor(bur_data$treatment, levels = c("Sun", "Shade"))
head(bur_data)

bur_data[, c("A4", "A42", "A5", "AA", "AB", "BA", "BB")] <- NULL

m_best <- "
bd ~ c(A, a) * B42
rdry ~ c(B, b) * bd + c(C, c) * B5
h3h1 ~ c(D, d) * bd + c(E, e) * rdry + c(F, f) * B5
idry ~ c(G, g)* bd + c(H, h) * rdry + c(I, i) * h3h1 + c(J, j) * B42 + c(K, k) * B5

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

#B42

#bd
B42_bd_dir_sun := a
B42_bd_dir_shade := A
B42_bd_dir_diff := A - a

B42_bd_ind_sun := 0
B42_bd_ind_shade := 0
B42_bd_ind_diff := 0

#rdry
B42_rdry_dir_sun := 0
B42_rdry_dir_shade := 0
B42_rdry_dir_diff := 0

B42_rdry_ind_sun := a * b
B42_rdry_ind_shade := A * B
B42_rdry_ind_diff := (A * B) - (a * b)

#h3h1
B42_h3h1_dir_sun := 0
B42_h3h1_dir_shade := 0
B42_h3h1_dir_diff := 0

B42_h3h1_ind_sun := 0
B42_h3h1_ind_shade := 0
B42_h3h1_ind_diff := 0

#idry
B42_idry_dir_sun := j
B42_idry_dir_shade := J
B42_idry_dir_diff := J - j

B42_idry_ind_sun := (a * g) + (a * b * h) + (a * d * i) + (a * b * e * i)
B42_idry_ind_shade := (A * G) + (A * B * H) + (A * D * I) + (A * B * E * I)
B42_idry_ind_diff := ((A * G) + (A * B * H) + (A * D * I) + (A * B * E * I)) - ((a * g) + (a * b * h) + (a * d * i) + (a * b * e * i))

#B5

#bd
B5_bd_dir_sun := 0
B5_bd_dir_shade := 0
B5_bd_dir_diff := 0

B5_bd_ind_sun := 0
B5_bd_ind_shade := 0
B5_bd_ind_diff := 0

#rdry
B5_rdry_dir_sun := c
B5_rdry_dir_shade := C
B5_rdry_dir_diff := C - c

B5_rdry_ind_sun := 0
B5_rdry_ind_shade := 0
B5_rdry_ind_diff := 0

#h3h1
B5_h3h1_dir_sun := f
B5_h3h1_dir_shade := F
B5_h3h1_dir_diff := F - f

B5_h3h1_ind_sun := c * e
B5_h3h1_ind_shade := C * E
B5_h3h1_ind_diff := (C * E) - (c * e)

#idry
B5_idry_dir_sun := k
B5_idry_dir_shade := K
B5_idry_dir_diff := K - k

B5_idry_ind_sun := (c * h) + (c * e * i) + (f * i)
B5_idry_ind_shade := (C * H) + (C * E * I) + (F * I)
B5_idry_ind_diff := ((C * H) + (C * E * I) + (F * I)) - ((c * h) + (c * e * i) + (f * i))
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
path_trait_correlations_sun$pop <- "bur"

path_parameters_trait_shade <- path_parameters_trait[str_is(path_parameters_trait$label, str_upper_case), ]
path_trait_correlations_shade <- data.frame(matrix(ncol = length(required_formula)))

for (i in 1:length(required_formula)) {
  path_trait_correlations_shade[, i] <- tryCatch(subset(path_parameters_trait_shade, formula == required_formula[i])$est, error = NA)
}
colnames(path_trait_correlations_shade) <- required_formula
path_trait_correlations_shade$env <- "Shade"
path_trait_correlations_shade$pop <- "bur"

path_parameters_eff <- subset(path_parameters, grepl("sun|shade|diff", lhs))
path_parameters_eff$allele_comb <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[1])
path_parameters_eff$trait <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[2])
path_parameters_eff$effect_type <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[3])
path_parameters_eff$env <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[4])
nrow(path_parameters_eff)
path_parameters_eff$pop <- "bur"

setwd("/Users/James/Documents/GitHub/sar_qtl/7_path_analysis/output/")
fwrite(path_parameters_eff, "bur_path_eff_alt.csv", sep = ",", row.names = FALSE)

fwrite(path_trait_correlations_sun, "bur_trait_eff_alt_sun.csv", sep = ",", row.names = FALSE, na = "NA")
fwrite(path_trait_correlations_shade, "bur_trait_eff_alt_shade.csv", sep = ",", row.names = FALSE, na = "NA")

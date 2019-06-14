#path analysis on bur group
library(data.table)
library(lavaan)
library(semPlot)
library(lettercase)

rm(list = ls())

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/7_path_analysis/output/")

cvi_data <- fread("cvi_col_pheno_qtl_combo_merge.csv",
                  sep = ",",
                  header = TRUE,
                  stringsAsFactors = FALSE)
setnames(cvi_data, c("r_dry", "h3_h1", "i_dry"), c("rdry", "h3h1", "idry"))
cvi_data$treatment <- factor(cvi_data$treatment, levels = c("Sun", "Shade"))
head(cvi_data)

cvi_data[, c("A4", "A42", "A5", "AA", "AB", "BA", "BB")] <- NULL

m_best <- "
bd ~ c(a, A) * B42 + c(b, B) * B5
rdry ~ c(c, C) * bd
h3h1 ~ c(d, D) * bd + c(e, E) * rdry + c(f, F) * B42 + c(g, G) * B5
idry ~ c(h, H) * bd + c(i, I) * rdry + c(j, J) * h3h1

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

B42_rdry_ind_sun := a * c
B42_rdry_ind_shade := A * C
B42_rdry_ind_diff := (A * C) - (a * c)

#h3h1
B42_h3h1_dir_sun := f
B42_h3h1_dir_shade := F
B42_h3h1_dir_diff := F - f

B42_h3h1_ind_sun := (a * d) + (a * c * e)
B42_h3h1_ind_shade := (A * D) + (A * C * E)
B42_h3h1_ind_diff := ((A * D) + (A * C * E)) - ((a * d) + (a * c * e))

#idry
B42_idry_dir_sun := 0
B42_idry_dir_shade := 0
B42_idry_dir_diff := 0

B42_idry_ind_sun := (a * h) + (a * c * i) + (a * d * j) + (a * c * e * j) + (f * j)
B42_idry_ind_shade := (A * H) + (A * C * I) + (A * D * J) + (A * C * E * J) + (F * J)
B42_idry_ind_diff := ((A * H) + (A * C * I) + (A * D * J) + (A * C * E * J) + (F * J)) - ((a * h) + (a * c * i) + (a * d * j) + (a * c * e * j) + (f * j))

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

B5_h3h1_ind_sun := (b * d) + (b * c * e)
B5_h3h1_ind_shade := (B * D) + (B * C * E)
B5_h3h1_ind_diff := ((B * D) + (B * C * E)) - ((b * d) + (b * c * e))

#idry
B5_idry_dir_sun := 0
B5_idry_dir_shade := 0
B5_idry_dir_diff := 0

B5_idry_ind_sun := (b * h) + (b * c * i) + (b * d * j) + (b * c * e * j) + (g * j)
B5_idry_ind_shade := (B * H) + (B * C * I) + (B * D * J) + (B * C * E * J) + (G * J)
B5_idry_ind_diff := ((B * H) + (B * C * I) + (B * D * J) + (B * C * E * J) + (G * J)) - ((b * h) + (b * c * i) + (b * d * j) + (b * c * e * j) + (g * j))
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
path_trait_correlations_sun$pop <- "cvi"

path_parameters_trait_shade <- path_parameters_trait[str_is(path_parameters_trait$label, str_upper_case), ]
path_trait_correlations_shade <- data.frame(matrix(ncol = length(required_formula)))

for (i in 1:length(required_formula)) {
    path_trait_correlations_shade[, i] <- tryCatch(subset(path_parameters_trait_shade, formula == required_formula[i])$est, error = NA)
}
colnames(path_trait_correlations_shade) <- required_formula
path_trait_correlations_shade$env <- "Shade"
path_trait_correlations_shade$pop <- "cvi"

path_parameters_eff <- subset(path_parameters, grepl("sun|shade|diff", lhs))
path_parameters_eff$allele_comb <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[1])
path_parameters_eff$trait <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[2])
path_parameters_eff$effect_type <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[3])
path_parameters_eff$env <- sapply(strsplit(path_parameters_eff$label, split = "_"), function(x) x[4])
nrow(path_parameters_eff)
path_parameters_eff$pop <- "cvi"

setwd("/Users/James/Documents/GitHub/sar_qtl/7_path_analysis/output/")
fwrite(path_parameters_eff, "cvi_path_eff_alt.csv", sep = ",", row.names = FALSE)

fwrite(path_trait_correlations_sun, "cvi_trait_eff_alt_sun.csv", sep = ",", row.names = FALSE, na = "NA")
fwrite(path_trait_correlations_shade, "cvi_trait_eff_alt_shade.csv", sep = ",", row.names = FALSE, na = "NA")


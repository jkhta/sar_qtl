# sar_qtl

# About
The code, figures, and supplementary data for "Shade alters indirect QTL effect sizes across vegetative and reproductive development in Arabidopsis thaliana" manuscript. The repository is split by the different steps of the QTL mapping pipeline, as well as the figures.
1_phenotype_standardization: Transforms and standardizes the phenotype data using the Box-Cox procedure.
2_lmm_models: Fits Bayesian mixed-models to the data.
3_lmm_blups: Extracts best-linear unbiased predictors (BLUPs) for each trait from the models.
4_geno_prob_array: Estimates genotype probabilities for each parent from genotype data.
5_nam_marker_positions: Creates a file
6_push_button: Creates the files for QTL mapping, generates permutations for QTL mapping, maps QTL using a stepwise approach, and calculates confidence intervals for each QTL.
7_path_analysis: Quantifies direct and indirect effects for each QTL, for each trait and within each population.
figures: Containts figures and supplemental data for the manuscript.


# Code

The R code for the QTL mapping pipeline was run on a cluster using the Slurm Workload Manager.
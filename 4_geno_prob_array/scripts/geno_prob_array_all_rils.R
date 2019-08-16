#this will generate a genotype probability matrix for the markers for Haley-Knott regression
library(qtl)
library(data.table)
library(lme4qtl)
library(lmerTest)
library(Matrix)
library(pedigreemm)
library(ggsci)
library(stringr)
library(parallel)
library(plyr)

rm(list = ls())

#reading in the new marker states for rqtl; then going to rename the genotypes
#so that they're easier to access
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/4_geno_prob_array/input/")
nam_pheno_geno <- fread("merged_NAM_lines_all_chr_11_NRP_joint_linkage_map_CSVR.csv",
                        sep = ",",
                        header = TRUE,
                        stringsAsFactors = FALSE)

#grabbing the genotype information
nam_geno <- data.frame(rils = colnames(nam_pheno_geno)[grepl("X", colnames(nam_pheno_geno))])

#removing the X in the genotype name
nam_geno$rils <- gsub("X", "", nam_geno$rils)

#getting the split between the rils in order to fill the rest with 0's
nam_geno_split <- strsplit(nam_geno$rils, split = "RV")
nam_geno$pop_no <- sapply(nam_geno_split, function(x) x[1])
nam_geno$ril_no <- sapply(nam_geno_split, function(x) x[2])
nam_geno$ril_no <- str_pad(nam_geno$ril_no, 3, side = "left", pad = 0)
nam_geno$geno <- with(nam_geno, paste(pop_no, "RV", ril_no, sep = ""))
nam_geno$geno <- factor(nam_geno$geno)
nam_geno$pop <- paste(nam_geno$pop_no, "RV", sep = "")
nam_geno$geno_pop <- with(nam_geno, paste(pop, geno, sep = "_"))

#reading in the marker data
cross <- read.cross('csvr','.','merged_NAM_lines_all_chr_11_NRP_joint_linkage_map_CSVR.csv',
                    genotype = c('AA','BB'),
                    na.strings = '-', 
                    crosstype = 'riself')

#grabbing the line names
cross_line_names <- as.character(cross$pheno$marker)
cross_line_names <- gsub("X", "", cross_line_names)

#adding an RV to the last numbers to match with the genotype data
cross_line_families <- sapply(strsplit(cross_line_names, "RV"), function(x) x[1])
cross_line_families <- paste(cross_line_families, "RV", sep = "")

#generating different genotype numbers to match the phenotype ril names and the marker
#ril names
cross_line_numbers <- sapply(strsplit(cross_line_names, "RV"), function(x) x[2])
cross_line_numbers <- str_pad(cross_line_numbers, 3, pad = "0")
cross_line_genotypes <- paste(cross_line_families, cross_line_numbers, sep = "")
cross_line_geno_pop <- paste(cross_line_families, cross_line_genotypes, sep = "_")

#checking to see if all lines are represented
all(cross_line_geno_pop %in% nam_geno$geno_pop)

#checking to see which lines aren't part of the rqtl markers
nam_geno_not_overlapping <- nam_geno$geno_pop[!(nam_geno$geno_pop %in% cross_line_geno_pop)]

#first calculating the genotype probabilities and then extract the probabilities from
#the cross object
cross <- calc.genoprob(cross, step = 0)
cross <- argmax.geno(cross)

#getting the geno argmax of the cross file in order to get information on the markers 
geno_argmax <- pull.argmaxgeno(cross, include.pos.info = T)
geno_argmax_info <- geno_argmax[, 1:3]
geno_argmax_geno_info <- geno_argmax[, 4:ncol(geno_argmax)]
colnames(geno_argmax_geno_info) <- paste(cross_line_families, cross_line_genotypes, sep = "_")

#seems like you can calculate both genotype probabilities and max arg for the cross
#object, and then pull them both out
geno_prob <- pull.genoprob(cross, omit.first.prob = TRUE)

#grabbing the chromosome, position and etc data
geno_prob <- data.table(geno_prob, keep.rownames = TRUE)

#renaming the column
colnames(geno_prob)[1] <- "geno"
geno_prob$geno <- paste(cross_line_families, cross_line_genotypes, sep = "_")

#changing the column names to be the marker names
#here i need the pulled argmax genotypes in order to label the marker names for
#the genotype probabilities
colnames(geno_prob)[2:ncol(geno_prob)] <- rownames(geno_argmax)

#grabbing the ril names from the column names from the argmax data frame
ril_names <- colnames(geno_argmax_geno_info)
rownames(geno_prob) <- ril_names
geno_prob <- data.table(geno_prob, keep.rownames = TRUE)
colnames(geno_prob)[1] <- "geno"

#need to transform the geno probabilities matrix into an array where each matrix is the probability
#that an individual is a particular genotype 
cross_line_code <- data.frame(family = unique(cross_line_families),
                              code = LETTERS[2:(1+length(unique(cross_line_families)))],
                              stringsAsFactors = FALSE)
cross_line_code$number <- 2:(nrow(cross_line_code) + 1)

#trying out a single pop vector
geno_pop <- sapply(strsplit(geno_prob$geno, split = "_"), function(x) x[1])

#converting pop to number based on code
geno_pop_to_no <- mapvalues(geno_pop, from = cross_line_code$family, to = cross_line_code$number)

geno_line_names <- geno_prob$geno

#changing genotype to number
geno_prob$geno <- as.numeric(geno_pop_to_no)

#this creates a genotype probability vector and then fills the vector by the 
#appropriate genotype probabilities in the correct vector slot
#so a genotype from family 2 will have the first and second slot filled out
geno_probability_vector_returner <- function(marker_prob) {
  empty_vector <- rep(0, 8)
  empty_vector[c(1, as.numeric(marker_prob[1]))] <- c(1 - marker_prob[2], marker_prob[2])
  empty_vector <- as.data.frame(t(as.data.frame(unlist(empty_vector))))
  return(empty_vector)
}

geno_prob <- as.data.frame(geno_prob)

#this generates a vector of 8 numbers for each line at the marker
marker_probability_dt_generator <- function(geno_cols) {
  marker_probability_dt <- rbindlist(lapply(1:nrow(geno_cols), function(x) geno_probability_vector_returner(geno_cols[x, ])))
  return(marker_probability_dt)
}

#making a genotype probability array variable from the genotype probabilities of the vectors
#this takes a while; it takes the first column and a genotype probability column
marker_probability_dt <- mclapply(2:ncol(geno_prob), 
                                  function(x) marker_probability_dt_generator(subset(geno_prob, select = c(1, x))),
                                  mc.cores = detectCores())

#renaming the elements of each list (a matrix of genotype probabilities) by the marker name
names(marker_probability_dt) <- colnames(geno_prob)[2:ncol(geno_prob)]
marker_probability_df <- lapply(marker_probability_dt, function(x) as.data.frame(x))

for (i in 1:length(marker_probability_df)) {
  rownames(marker_probability_df[[i]]) <- geno_line_names
}

#grabbing the genotypes that are just part of the kinship matrix
#don't need to run after saving the probability array
Genotypes <- marker_probability_df

#using the written rds file to read in the genotype probabilities array
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/4_geno_prob_array/output/")
saveRDS(Genotypes, "nam_rqtl_geno_prob_array_comp_11_all_rils_NEW.RDS")

#this script will verify the genotypes from published marker data to Marc's marker data
#populations are 21RV, 20RV, 8RV, 29RV, 28RV, 27RV, 13RV

rm(list = ls())

#first reading in the genotype probability file
setwd("/Users/James/Documents/GitHub/sar_qtl/6_push_button/")
geno_array <- readRDS("nam_rqtl_geno_prob_array_final_0.99_no_kinship.rds")
marker_anno <- fread("nam_marker_info_final.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
marker_anno$snp <- paste("m", marker_anno$snp, sep = "_")

#reading in the qtl data to extract the proper markers
setwd("/Users/James/Documents/GitHub/sar_qtl/figures/qtl_table/")
bd_qtl <- fread("bd_geno_qtl_ci.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
fri_flc_qtl <- c("m_4_407208", "m_5_3799350")

#subsetting from the genotype probability array FRI and FLC
names(geno_array) <- paste("m", names(geno_array), sep = "_")
geno_array_FRI_FLC <- geno_array[fri_flc_qtl]

library(data.table)

#grabbing genotypes for FRI and FLC individually to compare against the past genotype data
geno_array_FRI_FLC_col <- do.call(cbind, lapply(geno_array_FRI_FLC, function(x) ifelse(x[, 1] >= 0.5, "A", "B")))
geno_array_FRI <- data.frame(geno = rownames(geno_array_FRI_FLC_col), genotypes = geno_array_FRI_FLC_col[, 1])
geno_array_FLC <- data.frame(geno = rownames(geno_array_FRI_FLC_col), genotypes = geno_array_FRI_FLC_col[, 2])

#combining the QTL information with the cM information
marker_anno_FRI_FLC <- subset(marker_anno, snp %in% fri_flc_qtl)
marker_anno_FRI_FLC$Start_Mb <- as.numeric(sapply(strsplit(marker_anno_FRI_FLC$snp, split = "_"), function(x) x[3]))

#reading in the past genotypes data
setwd("/Users/James/Documents/GitHub/sar_qtl/genotype_validation/genotype/")
past_genotypes <- lapply(list.files(pattern = "RV_genotypes.csv"), function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE))
genotype_names <- sapply(strsplit(list.files(pattern = "RV_genotypes.csv"), split = "_"), function(x) x[1])

#changing C's and D's to NAs, and then adding genotype names to the columns
for (i in 1:length(past_genotypes)) {
    past_genotypes[[i]][past_genotypes[[i]] == "C"] <- "NA"
    past_genotypes[[i]][past_genotypes[[i]] == "D"] <- "NA"
    colnames(past_genotypes[[i]])[5:ncol(past_genotypes[[i]])] <- paste(paste(genotype_names[i], genotype_names[i], sep = "_"), colnames(past_genotypes[[i]])[5:ncol(past_genotypes[[i]])], sep = "")
}

#grabbing only the markers on chromosome 4
past_genotypes_chr_4 <- lapply(past_genotypes, function(x) subset(x, Chr == 4))

#calculating the differences between Marc's genotype data and the publicly available genotype data
for (i in 1:length(past_genotypes_chr_4)) {
    past_genotypes_chr_4[[i]]$cM_diff <- abs(past_genotypes_chr_4[[i]]$CR_cM - marker_anno_FRI_FLC$cM[1])
    past_genotypes_chr_4[[i]]$Mb_diff <- abs(past_genotypes_chr_4[[i]]$Start_Mb * 1000000 - marker_anno_FRI_FLC$Start_Mb[1])
}

past_genotypes_chr_4_cM_min <- lapply(past_genotypes_chr_4, function(x) subset(x, cM_diff == min(cM_diff))$cM_diff)
past_genotypes_chr_4_Mb_min <- lapply(past_genotypes_chr_4, function(x) subset(x, Mb_diff == min(Mb_diff))$Mb_diff)

#two different metrics of closeness: cM (genetic distance) or Mb (physical distance); I think physical distance would be a lot more consistent
past_genotypes_FRI <- lapply(past_genotypes_chr_4, function(x) t(subset(x, cM_diff == min(cM_diff))))
past_genotypes_FRI <- lapply(past_genotypes_chr_4, function(x) t(subset(x, Mb_diff == min(Mb_diff))))

#genotype plots
FRI_correlation <- c()

for (i in 1:length(past_genotypes_FRI)) {
    past_genotypes_FRI[[i]] <- data.frame(geno = rownames(past_genotypes_FRI[[i]]), genotypes = past_genotypes_FRI[[i]][, 1])
    FRI_merged <- merge(geno_array_FRI, past_genotypes_FRI[[i]], by = "geno")
    FRI_cor <- with(FRI_merged, cor(as.numeric(genotypes.x), as.numeric(genotypes.y), use = "pairwise.complete.obs"))
    single_pop_FRI_cor <- data.frame(pop = genotype_names[i], cor = FRI_cor, num_geno = nrow(FRI_merged))
    FRI_correlation <- rbindlist(list(FRI_correlation, single_pop_FRI_cor))
}

FRI_correlation

past_genotypes_chr_5 <- lapply(past_genotypes, function(x) subset(x, Chr == 5))

for (i in 1:length(past_genotypes_chr_5)) {
    past_genotypes_chr_5[[i]]$cM_diff <- abs(past_genotypes_chr_5[[i]]$CR_cM - marker_anno_FRI_FLC$cM[2])
    past_genotypes_chr_5[[i]]$Mb_diff <- abs(past_genotypes_chr_5[[i]]$Start_Mb * 1000000 - marker_anno_FRI_FLC$Start_Mb[2])
}

past_genotypes_chr_5_cM_min <- lapply(past_genotypes_chr_5, function(x) subset(x, cM_diff == min(cM_diff))$cM_diff)
past_genotypes_chr_5_Mb_min <- lapply(past_genotypes_chr_5, function(x) subset(x, Mb_diff == min(Mb_diff))$Mb_diff)

past_genotypes_FLC <- lapply(past_genotypes_chr_5, function(x) t(subset(x, cM_diff == min(cM_diff))))
past_genotypes_FLC <- lapply(past_genotypes_chr_5, function(x) t(subset(x, Mb_diff == min(Mb_diff))))

#genotype plots
FLC_correlation <- c()

for (i in 1:length(past_genotypes_FLC)) {
    past_genotypes_FLC[[i]] <- data.frame(geno = rownames(past_genotypes_FLC[[i]]), genotypes = past_genotypes_FLC[[i]][, 1])
    FLC_merged <- merge(geno_array_FRI, past_genotypes_FLC[[i]], by = "geno")
    FLC_cor <- with(FLC_merged, cor(as.numeric(genotypes.x), as.numeric(genotypes.y), use = "pairwise.complete.obs"))
    single_pop_FLC_cor <- data.frame(pop = genotype_names[i], cor = FLC_cor, num_geno = nrow(FLC_merged))
    FLC_correlation <- rbindlist(list(FLC_correlation, single_pop_FLC_cor))
}

FLC_correlation


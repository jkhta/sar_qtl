#this script will see if individual plant shade responses are similar across replications
library(data.table)

rm(list = ls())

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/2_lmm_models/input/")

nam_std_data <- fread("nam_cam_data_combined_tformed_std.csv",
                      sep = ",",
                      header = TRUE,
                      stringsAsFactors = FALSE)
nam_std_data_subset <- subset(nam_std_data, rep != 0)

nam_rep_list <- lapply(unique(nam_std_data_subset$rep), function(x) subset(nam_std_data_subset, rep == x))

nam_shade_response_by_rep <- c()

for (i in 1:length(nam_rep_list)) {
  nam_rep_single <- nam_rep_list[[i]]
  for (j in unique(nam_rep_single$geno)) {
    nam_rep_geno_sun <- subset(nam_rep_single, geno == j & treatment == "Sun")
    nam_rep_geno_shade <- subset(nam_rep_single, geno == j & treatment == "Shade")
    
    if (nrow(nam_rep_geno_sun) == 0  | nrow(nam_rep_geno_shade) == 0) {
      next
    }
    
    nam_rep_geno_sun_mean <- mean(nam_rep_geno_sun$bd)
    nam_rep_geno_shade_mean <- mean(nam_rep_geno_shade$bd)
    nam_shade_response <- data.table(geno = unique(nam_rep_geno_sun$geno), 
                                     rep = unique(nam_rep_geno_sun$rep), 
                                     bd_shade_resp = nam_rep_geno_shade_mean - nam_rep_geno_sun_mean)
    nam_shade_response_by_rep <- rbindlist(list(nam_shade_response_by_rep, nam_shade_response))
  }
}

nam_bd_by_rep_list <- lapply(unique(nam_shade_response_by_rep$rep), function(x) subset(nam_shade_response_by_rep, rep == x))

for (i in 1:length(nam_bd_by_rep_list)) {
  colnames(nam_bd_by_rep_list[[i]])[3] <- paste(colnames(nam_bd_by_rep_list[[i]])[3], i, sep = "_")
}

test <- Reduce(function(x, y) merge(x, y, by = "geno"), nam_bd_by_rep_list)
test_sar <- subset(test, select = grepl("shade_resp", colnames(test)))
cor(test_sar)
pairs(test_sar)
test_1_2 <- merge(nam_bd_by_rep_list[[1]], nam_bd_by_rep_list[[2]], by = "geno", all = FALSE)
plot(test_1_2$bd_shade_resp_1, test_1_2$bd_shade_resp_2)
test_1_3 <- merge(nam_bd_by_rep_list[[1]], nam_bd_by_rep_list[[3]], by = "geno", all = FALSE)
test_1_4 <- merge(nam_bd_by_rep_list[[1]], nam_bd_by_rep_list[[4]], by = "geno", all = FALSE)
test_2_3 <- merge(nam_bd_by_rep_list[[2]], nam_bd_by_rep_list[[3]], by = "geno", all = FALSE)
test_2_4 <- merge(nam_bd_by_rep_list[[2]], nam_bd_by_rep_list[[4]], by = "geno", all = FALSE)
test_3_4 <- merge(nam_bd_by_rep_list[[3]], nam_bd_by_rep_list[[4]], by = "geno", all = FALSE)
